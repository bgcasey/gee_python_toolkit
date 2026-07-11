# ---
# title:   Extract Image Collection Values to Features
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Reduce each image in an ee.ImageCollection over feature
#   regions (e.g., polygons or buffered points), copy image
#   properties onto the resulting features, rename extracted
#   properties with a reducer suffix, and export to Google
#   Drive as a CSV.
# ---

import ee


def image_collection_to_features(
    reducer,
    features,
    aoi,
    image_collection,
    crs,
    scale,
    tile_scale,
    file_name,
):
    """Reduce an image collection over features and export.

    Applies a reducer to feature regions (e.g., polygons or
    buffered points) for every image in the collection,
    renames the extracted image properties with a suffix, and
    exports the results to Google Drive as a CSV.

    Parameters
    ----------
    reducer : ee.Reducer
        Reducer function to apply.
    features : ee.FeatureCollection
        Collection of features (e.g., polygons, buffered
        points).
    aoi : ee.Geometry
        Area of interest used to filter features.
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
    # Step 1: Create suffix using reducer type
    reducer_type = reducer.getInfo()["type"].split(".").pop()
    suffix = ee.String(reducer_type)

    # Step 2: Filter features by AOI
    processed_features = features.filterBounds(aoi)

    # Step 3: Retrieve property names from features and images
    feature_properties = ee.Feature(
        features.first()
    ).propertyNames()
    img_properties = ee.Feature(
        image_collection.first()
    ).propertyNames()
    combined_properties = feature_properties.cat(img_properties)

    # Step 4: Reduce regions using the provided reducer,
    # and copy image properties onto each feature
    def reduce_image(img):
        return img.reduceRegions(
            collection=processed_features,
            crs=crs,
            reducer=reducer,
            scale=scale,
            tileScale=tile_scale,
        ).map(
            # Copy image properties (e.g., system:index)
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
