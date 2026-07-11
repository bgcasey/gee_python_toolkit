# ---
# title:   Extract Image Collection Values to Points
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Reduce each image in an ee.ImageCollection over point
#   locations (optionally buffered), copy image properties
#   onto the resulting features, rename extracted properties
#   with a reducer/buffer suffix, and export to Google Drive
#   as a CSV.
# ---

import ee


def image_collection_to_points(
    buffer_size,
    reducer,
    xy_points,
    aoi,
    image_collection,
    crs,
    scale,
    tile_scale,
    file_name,
):
    """Reduce an image collection over points and export.

    Applies a reducer to buffered point locations (or direct
    points when ``buffer_size`` is 0) for every image in the
    collection, renames the extracted image properties with a
    suffix, and exports the results to Google Drive as a CSV.

    Parameters
    ----------
    buffer_size : float
        Buffer size in meters around points. Use 0 to reduce
        directly at the point locations.
    reducer : ee.Reducer
        Reducer function to apply.
    xy_points : ee.FeatureCollection
        Collection of points.
    aoi : ee.Geometry
        Area of interest used to filter points.
    image_collection : ee.ImageCollection
        Image collection to sample.
    crs : str
        Coordinate reference system (CRS) to use.
    scale : float
        Scale in meters for reduction.
    tile_scale : float
        Tile scale for parallel processing.
    file_name : str
        Prefix for the exported file.

    Returns
    -------
    ee.FeatureCollection
        Feature collection with extracted image properties and
        renamed attributes.
    """
    # Step 1: Create suffix using buffer size and reducer type
    buffer_str = str(buffer_size)
    reducer_type = reducer.getInfo()["type"].split(".").pop()
    suffix = ee.String(reducer_type).cat("_").cat(buffer_str)

    # Step 2: Apply buffer to points or use them directly
    processed_points = xy_points.filterBounds(aoi).map(
        lambda pt: pt if buffer_size == 0 else pt.buffer(buffer_size)
    )

    # Step 3: Retrieve property names from points and images
    xy_properties = ee.Feature(xy_points.first()).propertyNames()
    img_properties = ee.Feature(
        image_collection.first()
    ).propertyNames()
    combined_properties = xy_properties.cat(img_properties)

    # Step 4: Reduce regions using the provided reducer,
    # and copy image properties onto each feature
    def reduce_image(img):
        return img.reduceRegions(
            collection=processed_points,
            crs=crs,
            reducer=reducer,
            scale=scale,
            tileScale=tile_scale,
        ).map(
            # Copy image properties (e.g., system:index) to the
            # feature
            lambda feature_with_reduction:
                feature_with_reduction.copyProperties(img)
        )

    reduced_region = image_collection.map(reduce_image).flatten()

    # Step 5: Rename extracted properties with a suffix
    def rename_properties(feature):
        new_properties = ee.Dictionary(
            feature.propertyNames()
            .map(
                lambda name: [
                    ee.Algorithms.If(
                        combined_properties.contains(name),
                        name,
                        ee.String(name).cat("_").cat(suffix),
                    ),
                    feature.get(name),
                ]
            )
            .flatten()
        )
        return ee.Feature(feature.geometry(), new_properties)

    renamed_feature_collection = reduced_region.map(
        rename_properties
    )

    # Step 6: Export the final feature collection to Google Drive
    task = ee.batch.Export.table.toDrive(
        collection=renamed_feature_collection,
        description=file_name,
        folder="gee_exports",
        fileNamePrefix=file_name,
        fileFormat="CSV",
    )
    task.start()

    return renamed_feature_collection
