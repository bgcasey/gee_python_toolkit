# ---
# title:   GEE Helper Functions
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   General-purpose helpers for working with Earth Engine
#   imagery: converting degrees to radians, combining and
#   normalizing images, building date lists, tiling an AOI,
#   filtering collections by band, reducing images and
#   collections to buffered points, exporting bands or
#   images to Drive, computing image statistics, and focal
#   (neighbourhood) statistics.
#
#   Several functions are adapted from the geeTools library
#   (https://github.com/aazuspan/geeTools).
# ---

import math

import ee


def deg2rad(deg):
    """Convert degrees to radians.

    Parameters
    ----------
    deg : ee.Number or ee.Image
        An angle in degrees.

    Returns
    -------
    ee.Number or ee.Image
        The angle in radians.
    """
    coeff = 180 / math.pi
    return deg.divide(coeff)


def combine_images(img_list, optional_parameters=None):
    """Combine a list of images into one multi-band image.

    Convenience wrapper over repeatedly calling addBands.

    Parameters
    ----------
    img_list : ee.List
        A list of images to combine. Images can be single-
        or multi-band.
    optional_parameters : dict, optional
        Optional parameters to override defaults:
        - prefix (bool): If True (default), band names are
          prefixed with the list index of the source image,
          allowing images with identical band names to be
          combined. If False, original band names are kept
          (an error is raised on duplicate names).
        - props (ee.Dictionary): Properties to store in the
          combined image. Defaults to the properties of the
          first image in img_list.

    Returns
    -------
    ee.Image
        An image with the bands of all images in img_list.
    """
    first = ee.Image(ee.List(img_list).get(0))

    # Start from defaults, then override with any provided
    # optional parameters.
    params = {
        "prefix": True,
        "props": first.toDictionary(first.propertyNames()),
    }
    if optional_parameters:
        params.update(optional_parameters)

    # Convert the list to a collection and collapse it into
    # a multiband image, storing the chosen properties.
    combined = (
        ee.ImageCollection.fromImages(img_list)
        .toBands()
        .set(params["props"])
    )

    if params["prefix"] is False:
        band_names = ee.List(
            img_list.map(lambda img: img.bandNames())
        ).flatten()
        combined = combined.rename(band_names)

    return combined


def normalize_image(img, optional_parameters=None):
    """Normalize each band of an image to the range 0-1.

    Parameters
    ----------
    img : ee.Image
        The image to normalize.
    optional_parameters : dict, optional
        Optional parameters to override defaults:
        - region (ee.Geometry): Area over which to compute
          image statistics. Defaults to None.
        - scale (float): Scale, in image units, at which to
          compute statistics. Defaults to None.
        - max_pixels (float): Maximum number of pixels to
          sample. Defaults to 1e13.

    Returns
    -------
    ee.Image
        The input image with all bands rescaled 0-1.
    """
    # Start from defaults, then override with any provided
    # optional parameters.
    params = {
        "region": None,
        "scale": None,
        "max_pixels": 1e13,
    }
    if optional_parameters:
        params.update(optional_parameters)

    min_img = img.reduceRegion(
        reducer=ee.Reducer.min(),
        geometry=params["region"],
        scale=params["scale"],
        maxPixels=params["max_pixels"],
    ).toImage(img.bandNames())

    max_img = img.reduceRegion(
        reducer=ee.Reducer.max(),
        geometry=params["region"],
        scale=params["scale"],
        maxPixels=params["max_pixels"],
    ).toImage(img.bandNames())

    return img.subtract(min_img).divide(
        max_img.subtract(min_img)
    )


def create_date_list(date_start, date_end, interval,
                     interval_type):
    """Generate a list of dates for time series analysis.

    Dates start from the given start date and advance by a
    fixed interval until the end date.

    Parameters
    ----------
    date_start : ee.Date
        Start date of the time series.
    date_end : ee.Date
        End date of the time series.
    interval : int
        Units to skip between dates in the series.
    interval_type : str
        Type of interval ('months', 'weeks', 'days',
        'years').

    Returns
    -------
    ee.List
        A list of dates for the time series.
    """
    # Total intervals between the start and end dates.
    n_intervals = date_end.difference(
        date_start, interval_type
    ).round()

    # Sequence from 0 to n_intervals, stepping by interval.
    dates = ee.List.sequence(0, n_intervals, interval)

    # Advance the start date by n intervals.
    def make_datelist(n):
        return date_start.advance(n, interval_type)

    dates = dates.map(make_datelist)

    return dates


def get_vis_params(image, band, aoi, scale):
    """Get min/max visualization parameters for a band.

    Parameters
    ----------
    image : ee.Image
        Image containing the band.
    band : str
        Name of the band.
    aoi : ee.Geometry
        Area of interest.
    scale : float
        Scale for the reduceRegion operation.

    Returns
    -------
    dict
        Visualization parameters with 'min', 'max', and a
        red-yellow-green 'palette'.
    """
    stats = image.select(band).reduceRegion(
        reducer=ee.Reducer.minMax(),
        geometry=aoi,
        scale=scale,
        bestEffort=True,
        tileScale=8,
    )

    min_val = stats.get(band + "_min").getInfo()
    max_val = stats.get(band + "_max").getInfo()

    return {
        "min": min_val,
        "max": max_val,
        "palette": ["red", "yellow", "green"],
    }


def split_aoi_into_tiles(aoi, tile_size):
    """Split an AOI into square tiles.

    Parameters
    ----------
    aoi : ee.Geometry
        Area of interest.
    tile_size : float
        Side length of each tile, in meters.

    Returns
    -------
    ee.List
        Tiles as ee.Geometry.Rectangle objects.
    """
    aoi = ee.Geometry(aoi)

    # Bounds of the AOI and its corner coordinates.
    bounds = aoi.bounds()
    coords = ee.List(bounds.coordinates().get(0))

    bottom_left = ee.Geometry.Point(coords.get(0))
    top_right = ee.Geometry.Point(coords.get(2))

    # Width and height of the AOI in meters.
    aoi_width = bottom_left.distance(
        ee.Geometry.Point([
            top_right.coordinates().get(0),
            bottom_left.coordinates().get(1),
        ])
    )
    aoi_height = bottom_left.distance(
        ee.Geometry.Point([
            bottom_left.coordinates().get(0),
            top_right.coordinates().get(1),
        ])
    )

    # Number of tiles needed horizontally and vertically.
    num_tiles_x = aoi_width.divide(tile_size).ceil()
    num_tiles_y = aoi_height.divide(tile_size).ceil()

    def make_column(i):
        def make_tile(j):
            x = ee.Number(
                bottom_left.coordinates().get(0)
            ).add(ee.Number(i).multiply(tile_size))
            y = ee.Number(
                bottom_left.coordinates().get(1)
            ).add(ee.Number(j).multiply(tile_size))
            return ee.Geometry.Rectangle([
                x, y, x.add(tile_size), y.add(tile_size)
            ])

        return ee.List.sequence(
            0, num_tiles_y.subtract(1)
        ).map(make_tile)

    tiles = ee.List.sequence(
        0, num_tiles_x.subtract(1)
    ).map(make_column).flatten()

    return tiles


def filter_collection_by_bands(collection, required_bands):
    """Filter a collection to images with all given bands.

    Parameters
    ----------
    collection : ee.ImageCollection
        The image collection to filter.
    required_bands : list of str
        The required band names.

    Returns
    -------
    ee.ImageCollection
        The filtered image collection.
    """
    def has_all_required_bands(image):
        band_names = image.bandNames()
        # min ensures every required band is present.
        has_all_bands = ee.List(required_bands).map(
            lambda band: band_names.contains(band)
        ).reduce(ee.Reducer.min())
        return image.set("hasAllBands", has_all_bands)

    collection_with_check = collection.map(
        has_all_required_bands
    )
    filtered = collection_with_check.filterMetadata(
        "hasAllBands", "equals", True
    )

    return filtered


def image_to_points(buffer_size, reducer, xy_points, aoi,
                    image, crs, scale, tile_scale,
                    file_name):
    """Reduce an image to buffered points and export as CSV.

    Buffers points (unless buffer_size is 0), applies the
    reducer to the resulting regions, renames the reduced
    properties with a reducer/buffer suffix, and exports the
    result to Google Drive.

    Parameters
    ----------
    buffer_size : float
        Buffer size to apply to points, in meters. Use 0 to
        reduce at the points directly.
    reducer : ee.Reducer
        Reducer to apply to the regions.
    xy_points : ee.FeatureCollection
        Points for analysis.
    aoi : ee.Geometry
        Area of interest used to filter points.
    image : ee.Image
        Image to reduce.
    crs : str
        Coordinate reference system to use.
    scale : float
        Scale, in meters, for the reduction.
    tile_scale : float
        Tile scale for parallel processing.
    file_name : str
        Prefix for the exported file.

    Returns
    -------
    ee.FeatureCollection
        The reduced collection with renamed properties.
    """
    buffer_str = str(buffer_size)

    # Build a suffix from the reducer type and buffer size.
    reducer_info = reducer.getInfo()
    reducer_type = reducer_info["type"].split(".")[-1]
    suffix = ee.String(reducer_type).cat("_").cat(buffer_str)

    # Buffer points unless buffer_size is 0.
    def buffer_point(pt):
        return pt if buffer_size == 0 else pt.buffer(
            buffer_size
        )

    processed_points = xy_points.filterBounds(aoi).map(
        buffer_point
    )

    # Property names to preserve from the input points.
    xy_properties = ee.Feature(
        xy_points.first()
    ).propertyNames()

    reduced_regions = image.reduceRegions(
        collection=processed_points,
        reducer=reducer,
        crs=crs,
        scale=scale,
        tileScale=tile_scale,
    )

    # Rename reduced properties, keeping the point ones.
    def rename_properties(feature):
        def rename_one(name):
            new_name = ee.Algorithms.If(
                xy_properties.contains(name),
                name,
                ee.String(name).cat("_").cat(suffix),
            )
            return [new_name, feature.get(name)]

        new_properties = ee.Dictionary(
            feature.propertyNames().map(rename_one).flatten()
        )
        return ee.Feature(feature.geometry(), new_properties)

    renamed_fc = reduced_regions.map(rename_properties)

    # Export the result to Google Drive.
    task = ee.batch.Export.table.toDrive(
        collection=renamed_fc,
        description=file_name,
        folder="gee_exports",
        fileNamePrefix=file_name,
        fileFormat="CSV",
    )
    task.start()

    return renamed_fc


def image_collection_to_points(buffer_size, reducer,
                               xy_points, aoi,
                               image_collection, crs, scale,
                               tile_scale, file_name):
    """Reduce a collection to buffered points, export CSV.

    Buffers points (unless buffer_size is 0), applies the
    reducer to each image in the collection, renames the
    reduced properties with a reducer/buffer suffix, and
    exports the result to Google Drive.

    Parameters
    ----------
    buffer_size : float
        Buffer size to apply to points, in meters. Use 0 to
        reduce at the points directly.
    reducer : ee.Reducer
        Reducer to apply to the regions.
    xy_points : ee.FeatureCollection
        Points for analysis.
    aoi : ee.Geometry
        Area of interest used to filter points.
    image_collection : ee.ImageCollection
        Image collection to reduce.
    crs : str
        Coordinate reference system to use.
    scale : float
        Scale, in meters, for the reduction.
    tile_scale : float
        Tile scale for parallel processing.
    file_name : str
        Prefix for the exported file.

    Returns
    -------
    ee.FeatureCollection
        The reduced collection with renamed properties.
    """
    buffer_str = str(buffer_size)

    # Build a suffix from the reducer type and buffer size.
    reducer_info = reducer.getInfo()
    reducer_type = reducer_info["type"].split(".")[-1]
    suffix = ee.String(reducer_type).cat("_").cat(buffer_str)

    image_collection = ee.ImageCollection(image_collection)

    # Buffer points unless buffer_size is 0.
    def buffer_point(pt):
        return pt if buffer_size == 0 else pt.buffer(
            buffer_size
        )

    processed_points = xy_points.filterBounds(aoi).map(
        buffer_point
    )

    # Property names to preserve from points and images.
    xy_properties = ee.Feature(
        xy_points.first()
    ).propertyNames()
    img_properties = ee.Feature(
        image_collection.first()
    ).propertyNames()
    combined_properties = xy_properties.cat(img_properties)

    # Reduce each image and copy its properties.
    def reduce_image(img):
        return img.reduceRegions(
            collection=processed_points,
            crs=crs,
            reducer=reducer,
            scale=scale,
            tileScale=tile_scale,
        ).map(lambda f: f.copyProperties(img))

    reduced_region = image_collection.map(
        reduce_image
    ).flatten()

    # Rename reduced properties, keeping the known ones.
    def rename_properties(feature):
        def rename_one(name):
            new_name = ee.Algorithms.If(
                combined_properties.contains(name),
                name,
                ee.String(name).cat("_").cat(suffix),
            )
            return [new_name, feature.get(name)]

        new_properties = ee.Dictionary(
            feature.propertyNames().map(rename_one).flatten()
        )
        return ee.Feature(feature.geometry(), new_properties)

    renamed_fc = reduced_region.map(rename_properties)

    # Export the result to Google Drive.
    task = ee.batch.Export.table.toDrive(
        collection=renamed_fc,
        description=file_name,
        folder="gee_exports",
        fileNamePrefix=file_name,
        fileFormat="CSV",
    )
    task.start()

    return renamed_fc


def export_bands_by_year(collection, aoi, folder, scale,
                         crs, file_name_fn):
    """Export each band of each image to Google Drive.

    Iterates over the collection client-side, clips each
    band to the AOI, and starts a Drive export task per
    band.

    Parameters
    ----------
    collection : ee.ImageCollection
        The image collection to process.
    aoi : ee.Geometry
        Area of interest to clip the images.
    folder : str
        Google Drive folder name for exports.
    scale : float
        Export scale in meters.
    crs : str
        Output coordinate reference system.
    file_name_fn : callable
        Function that returns a file name given an image
        and a band name.
    """
    col_list = collection.toList(collection.size())
    size = collection.size().getInfo()

    for i in range(size):
        try:
            img = ee.Image(col_list.get(i))
            bands = img.bandNames().getInfo()

            for band in bands:
                band_image = img.select(band).clip(aoi)
                file_name = file_name_fn(img, band)

                if not file_name or not isinstance(
                    file_name, str
                ):
                    raise ValueError(
                        "Invalid file name generated."
                    )

                task = ee.batch.Export.image.toDrive(
                    image=band_image,
                    description=file_name,
                    folder=folder,
                    fileNamePrefix=file_name,
                    region=aoi,
                    scale=scale,
                    crs=crs,
                    maxPixels=1e13,
                )
                task.start()
        except Exception as err:
            print(f"Error processing image: {err}")
            continue


def export_image_collection(collection, aoi, folder, scale,
                            crs, file_name_fn):
    """Export each image as a multi-band GeoTIFF to Drive.

    Iterates over the collection client-side, clips each
    image to the AOI, and starts a Drive export task per
    image.

    Parameters
    ----------
    collection : ee.ImageCollection
        The image collection to process.
    aoi : ee.Geometry
        Area of interest to clip the images.
    folder : str
        Google Drive folder name for exports.
    scale : float
        Export scale in meters.
    crs : str
        Output coordinate reference system.
    file_name_fn : callable
        Function that returns a file name given an image.
    """
    col_list = collection.toList(collection.size())
    size = collection.size().getInfo()

    for i in range(size):
        try:
            img = ee.Image(col_list.get(i))
            file_name = file_name_fn(img)

            if not file_name or not isinstance(
                file_name, str
            ):
                raise ValueError(
                    "Invalid file name generated."
                )

            clipped_image = img.clip(aoi)

            task = ee.batch.Export.image.toDrive(
                image=clipped_image,
                description=file_name,
                folder=folder,
                fileNamePrefix=file_name,
                region=aoi,
                scale=scale,
                crs=crs,
                maxPixels=1e13,
            )
            task.start()
        except Exception as err:
            print(f"Error processing image: {err}")
            continue


def calculate_image_stats(image, geometry, scale,
                          max_pixels, reducer):
    """Calculate statistics for a single image.

    Parameters
    ----------
    image : ee.Image
        Image for which statistics are calculated.
    geometry : ee.Geometry
        Geometry defining the area of interest.
    scale : float
        Scale (in meters) for the reducer.
    max_pixels : float
        Maximum number of pixels to process.
    reducer : ee.Reducer
        Reducer used to calculate statistics.

    Returns
    -------
    ee.Dictionary
        The calculated statistics.
    """
    return image.reduceRegion(
        reducer=reducer,
        geometry=geometry,
        scale=scale,
        bestEffort=True,
        maxPixels=max_pixels,
    )


def calculate_image_collection_stats(collection, geometry,
                                     scale, max_pixels,
                                     reducer):
    """Append statistics to each image in a collection.

    Parameters
    ----------
    collection : ee.ImageCollection
        The collection of images.
    geometry : ee.Geometry
        Geometry defining the area of interest.
    scale : float
        Scale (in meters) for the reducer.
    max_pixels : float
        Maximum number of pixels to process.
    reducer : ee.Reducer
        Reducer used to calculate statistics.

    Returns
    -------
    ee.ImageCollection
        The image collection with appended statistics.
    """
    def add_stats(image):
        stats = calculate_image_stats(
            image, geometry, scale, max_pixels, reducer
        )
        return image.set(stats)

    return collection.map(add_stats)


def export_stats_to_csv(stats_collection, file_name):
    """Export image statistics to a CSV on Google Drive.

    Parameters
    ----------
    stats_collection : ee.FeatureCollection
        Feature collection containing the statistics.
    file_name : str
        Name of the CSV file.
    """
    task = ee.batch.Export.table.toDrive(
        collection=stats_collection,
        description=file_name,
        folder="gee_tables",
        fileFormat="CSV",
    )
    task.start()


def focal_stats(image, kernel_size, shape,
                properties_to_copy=None):
    """Apply a focal mean and rename bands by kernel size.

    Parameters
    ----------
    image : ee.Image
        Input image to which the focal mean is applied.
    kernel_size : float
        Radius of the kernel, in meters.
    shape : str
        Kernel shape (e.g., 'circle' or 'square').
    properties_to_copy : list of str, optional
        Property names to copy from the original image.
        Defaults to ['system:time_start'].

    Returns
    -------
    ee.Image
        Focal-mean image with bands renamed to include the
        kernel size, plus the copied properties.
    """
    if not properties_to_copy:
        properties_to_copy = ["system:time_start"]

    # Focal mean using the requested kernel shape and size.
    focal = image.reduceNeighborhood(
        reducer=ee.Reducer.mean(),
        kernel=getattr(ee.Kernel, shape)(
            kernel_size, "meters"
        ),
    )

    band_names = focal.bandNames()

    # Append the kernel size to each band name.
    def append_size(band_name):
        return ee.String(band_name).cat("_").cat(
            ee.Number(kernel_size).format()
        )

    new_names = band_names.map(append_size)
    renamed = focal.rename(new_names)

    # Copy the requested properties from the source image.
    result = renamed.copyProperties(
        image, properties_to_copy
    )

    return result
