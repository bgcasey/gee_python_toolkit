# ---
# title:   Extract Image Values to Points
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Reduce a single ee.Image over point locations (optionally
#   buffered) and export the results to Google Drive as a CSV.
#   Extracted properties are renamed with a suffix built from
#   the reducer type and buffer size. Non-XY properties are
#   suffixed; original point properties are preserved.
# ---

import ee


def image_to_points(
    buffer_size,
    reducer,
    xy_points,
    aoi,
    image,
    crs,
    scale,
    tile_scale,
    file_name,
):
    """Reduce an image over points and export to Drive.

    Optionally buffers points, applies a reducer to the
    buffered regions (or directly to the points when
    ``buffer_size`` is 0), renames the resulting properties
    with a suffix, and exports the collection to Google Drive
    as a CSV.

    Parameters
    ----------
    buffer_size : float
        Buffer size (in meters) to apply to points. Use 0 to
        reduce directly at the point locations.
    reducer : ee.Reducer
        The reducer to apply to the buffered regions.
    xy_points : ee.FeatureCollection
        The collection of points for analysis.
    aoi : ee.FeatureCollection
        Area of interest used to filter points.
    image : ee.Image
        The image to sample at each point.
    crs : str
        Coordinate reference system (e.g., 'EPSG:4326').
    scale : float
        Nominal scale in meters of the projection to work at.
    tile_scale : float
        Scaling factor for large parallel computations.
    file_name : str
        Prefix for the exported file in Google Drive.

    Returns
    -------
    ee.FeatureCollection
        The collection of features with renamed properties.
    """
    # Step 1: Prepare a string suffix based on buffer size and
    # reducer type.
    buffer_str = str(buffer_size)
    reducer_info = reducer.getInfo()
    reducer_type = reducer_info["type"].split(".").pop()
    suffix = ee.String(reducer_type).cat("_").cat(buffer_str)

    # Step 2: Filter the input points by the AOI and buffer them
    # if buffer_size is not 0.
    processed_points = xy_points.filterBounds(aoi).map(
        lambda pt: pt if buffer_size == 0 else pt.buffer(buffer_size)
    )

    # Step 3: Reduce the image over the processed points.
    reduced_regions = image.reduceRegions(
        collection=processed_points,
        reducer=reducer,
        crs=crs,
        scale=scale,
        tileScale=tile_scale,
    )

    # Step 4: Rename properties of each feature to include the
    # suffix for non-XY properties.
    xy_properties = ee.Feature(xy_points.first()).propertyNames()

    def rename_properties(feature):
        new_properties = ee.Dictionary(
            feature.propertyNames()
            .map(
                lambda name: [
                    ee.Algorithms.If(
                        xy_properties.contains(name),
                        name,
                        ee.String(name).cat("_").cat(suffix),
                    ),
                    feature.get(name),
                ]
            )
            .flatten()
        )
        return ee.Feature(feature.geometry(), new_properties)

    renamed_feature_collection = reduced_regions.map(
        rename_properties
    )

    # Step 5: Export the resulting collection to Google Drive and
    # return it.
    task = ee.batch.Export.table.toDrive(
        collection=renamed_feature_collection,
        description=file_name,
        folder="gee_exports",
        fileNamePrefix=file_name,
        fileFormat="CSV",
    )
    task.start()

    return renamed_feature_collection
