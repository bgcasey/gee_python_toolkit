# ---
# title:   Gap Filling Functions
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Fill gaps in Earth Engine imagery using Inverse
#   Distance Weighting (IDW) interpolation, applied band by
#   band. Each band is sampled to points, interpolated, and
#   recombined into a single multi-band image.
# ---

import ee


def apply_idw_interpolation(image, aoi, range_m, gamma,
                            num_pixels):
    """Fill image gaps with IDW interpolation.

    Parameters
    ----------
    image : ee.Image
        Input image with gaps to fill.
    aoi : ee.Geometry
        Area of interest for interpolation.
    range_m : float
        Maximum distance (in meters) to search for values.
    gamma : float
        Decay factor for the inverse distance weighting.
    num_pixels : int
        Number of pixels to sample for interpolation.

    Returns
    -------
    ee.Image
        The image with gaps filled by interpolation.
    """
    band_names = image.bandNames()

    # Interpolate a single band by name.
    def interpolate_band(band_name):
        band_name = ee.String(band_name)

        # Turn each sampled pixel into a point feature.
        def to_point(sample):
            lat = sample.get("latitude")
            lon = sample.get("longitude")
            value = sample.get(band_name)
            return ee.Feature(
                ee.Geometry.Point([lon, lat])
            ).set(band_name, value)

        # Sample the band to get known values.
        samples = (
            image.select([band_name])
            .addBands(ee.Image.pixelLonLat())
            .sample(
                region=aoi,
                numPixels=num_pixels,
                scale=30,
                projection="EPSG:4326",
            )
            .map(to_point)
        )

        # Global mean and standard deviation of samples.
        stats = samples.reduceColumns(
            reducer=ee.Reducer.mean().combine(
                reducer2=ee.Reducer.stdDev(),
                sharedInputs=True,
            ),
            selectors=[band_name],
        )

        # Apply IDW interpolation.
        interpolated = samples.inverseDistance(
            range=range_m,
            propertyName=band_name,
            mean=stats.get("mean"),
            stdDev=stats.get("stdDev"),
            gamma=gamma,
        )

        return interpolated.rename(band_name)

    # Interpolate every band and combine into one image.
    interpolated_bands = band_names.map(interpolate_band)

    interpolated_image = (
        ee.ImageCollection(interpolated_bands)
        .toBands()
        .clip(aoi)
    )

    return interpolated_image
