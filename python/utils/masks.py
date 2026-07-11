# ---
# title:   Image Masking Functions
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Masking helpers for Earth Engine imagery: mask by
#   forest age, Sentinel-2 cloud / snow / vegetation /
#   water masks, land cover masking, and a Dynamic World
#   tree mask. Each function takes an ee.Image and returns
#   the masked ee.Image.
# ---

import ee


def mask_by_forest_age(image):
    """Mask an image to forest older than a threshold.

    Uses the CA_forest_age_2019 dataset and keeps pixels
    with a forest age greater than 60 years.

    Parameters
    ----------
    image : ee.Image
        Image to mask.

    Returns
    -------
    ee.Image
        The masked image.
    """
    age = ee.Image(
        "projects/sat-io/open-datasets/CA_FOREST/"
        "CA_forest_age_2019"
    )

    # Keep pixels older than the age threshold.
    age_threshold = 60
    mask = age.gt(age_threshold)

    return image.updateMask(mask)


def mask_s2_clouds(image):
    """Mask clouds and cirrus in Sentinel-2 using QA60.

    Parameters
    ----------
    image : ee.Image
        Sentinel-2 image with a QA60 band.

    Returns
    -------
    ee.Image
        Cloud-masked image scaled to reflectance (0-1),
        preserving 'system:time_start'.
    """
    qa = image.select("QA60")

    # Bits 10 and 11 are clouds and cirrus, respectively.
    cloud_bit_mask = 1 << 10
    cirrus_bit_mask = 1 << 11

    # Both flags must be zero for clear conditions.
    mask = qa.bitwiseAnd(cloud_bit_mask).eq(0).And(
        qa.bitwiseAnd(cirrus_bit_mask).eq(0)
    )

    return (
        image.updateMask(mask)
        .divide(10000)
        .copyProperties(image, ["system:time_start"])
    )


def mask_s2_snow(image):
    """Mask snow pixels in Sentinel-2 using the SCL band.

    Parameters
    ----------
    image : ee.Image
        Sentinel-2 image with an SCL band.

    Returns
    -------
    ee.Image
        Image with snow pixels (SCL class 11) masked out.
    """
    scl = image.select("SCL")
    non_snow_mask = scl.neq(11)

    return image.updateMask(non_snow_mask)


def mask_s2_vegetation(image):
    """Mask a Sentinel-2 image to vegetation pixels.

    Parameters
    ----------
    image : ee.Image
        Sentinel-2 image with an SCL band.

    Returns
    -------
    ee.Image
        Image masked to vegetation pixels (SCL class 4).
    """
    scl = image.select("SCL")
    vegetation_mask = scl.eq(4)

    return image.updateMask(vegetation_mask)


def mask_s2_water(image):
    """Mask out water pixels using the SCL band.

    Parameters
    ----------
    image : ee.Image
        Sentinel-2 image with an SCL band.

    Returns
    -------
    ee.Image
        Image with water pixels (SCL class 6) masked out.
    """
    scl = image.select("SCL")
    wanted_pixels = scl.neq(6)

    return image.updateMask(wanted_pixels)


def mask_by_landcover(image):
    """Mask an image to a target land cover class.

    Uses the CA_FOREST_LC_VLCE2 land cover product (30 m)
    for 2019 and keeps pixels of the target class (210).

    Parameters
    ----------
    image : ee.Image
        Image to mask.

    Returns
    -------
    ee.Image
        The masked image.
    """
    mask_year = 2019
    mask_collection = ee.ImageCollection(
        "projects/sat-io/open-datasets/CA_FOREST_LC_VLCE2"
    ).filter(
        ee.Filter.calendarRange(mask_year, mask_year, "year")
    )

    # Land cover class to keep.
    landcover_class = 210
    mask_image = mask_collection.first()
    mask = mask_image.eq(landcover_class)

    return image.updateMask(mask)


def dynamic_world(image):
    """Mask an image to tree pixels using Dynamic World.

    Builds a mean probability composite from Dynamic World
    for the summer (June-September) of the image's year and
    keeps pixels whose most probable class is 'trees'.

    Parameters
    ----------
    image : ee.Image
        Image to mask. Must carry a 'date' property.

    Returns
    -------
    ee.Image
        Image masked to tree pixels.
    """
    probability_bands = [
        "water",
        "trees",
        "grass",
        "flooded_vegetation",
        "crops",
        "shrub_and_scrub",
        "built",
        "bare",
        "snow_and_ice",
    ]

    # One-year window starting at the image's date.
    start = ee.Date(image.get("date"))
    end = ee.Date(image.get("date")).advance(1, "year")
    date_range = ee.DateRange(start, end)

    dw = (
        ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
        .filterDate(date_range)
        .filter(ee.Filter.calendarRange(6, 9, "month"))
    )

    # Mean probability per class over the time period.
    dw_time_series = dw.select(probability_bands)
    mean_probability = dw_time_series.reduce(
        reducer=ee.Reducer.mean(),
        parallelScale=10,
    )

    # Class with the top mean probability per pixel.
    top_probability = (
        mean_probability.toArray()
        .arrayArgmax()
        .arrayGet(0)
        .rename("label")
    )

    # Class index 1 corresponds to 'trees'.
    mask = top_probability.eq(1)

    return image.updateMask(mask)
