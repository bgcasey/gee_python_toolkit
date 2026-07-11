# ---
# title:   Landsat Indices and Masks Functions
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Functions to calculate spectral indices and apply masks
#   to a time series of Landsat images. Indices include
#   vegetation, moisture, and stress-related measures. Masks
#   handle cloud, snow, fill, saturation, and QA filtering.
# ---

import ee

from utils import annual_forest_land_cover as forest_lc
from utils import masks

# Band name and threshold used to classify stressed forest
# pixels from the NDRS index.
BAND_NAME = "NDRS"
THRESHOLD = 0.5


def add_bsi(image):
    """Add a Bare Soil Index (BSI) band to an image.

    BSI = ((Red + SWIR) - (NIR + Blue)) /
          ((Red + SWIR) + (NIR + Blue))

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the BSI band added.
    """
    bsi = image.expression(
        "((Red + SWIR) - (NIR + Blue)) / "
        "((Red + SWIR) + (NIR + Blue))",
        {
            "NIR": image.select("SR_B4"),
            "Red": image.select("SR_B3"),
            "Blue": image.select("SR_B1"),
            "SWIR": image.select("SR_B5"),
        },
    ).rename("BSI")
    return image.addBands([bsi])


def add_dswi(image):
    """Add a Disease Stress Water Index (DSWI) band.

    DSWI = (NIR + Green) / (SWIR + Red), clamped to [0, 3].

    Galvao, L. S., Formaggio, A. R., and Tisot, D. A.
    (2005). Discrimination of sugarcane varieties in
    Southeastern Brazil with EO-1 Hyperion data. Remote
    Sens. Environ. 94, 523-534.

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the DSWI band added.
    """
    dswi = image.expression(
        "(NIR + Green) / (SWIR + Red)",
        {
            "NIR": image.select("SR_B4"),
            "Green": image.select("SR_B2"),
            "SWIR": image.select("SR_B5"),
            "Red": image.select("SR_B3"),
        },
    ).rename("DSWI")
    dswi_clamped = dswi.clamp(0, 3)
    return image.addBands(dswi_clamped)


def add_drs(image):
    """Add a Distance Red & SWIR (DRS) band to an image.

    DRS = sqrt((RED^2) + (SWIR^2))

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the DRS band added.
    """
    drs = image.expression(
        "sqrt(((RED) * (RED)) + ((SWIR) * (SWIR)))",
        {
            "SWIR": image.select("SR_B5"),
            "RED": image.select("SR_B3"),
        },
    ).rename("DRS")
    return image.addBands([drs])


def add_evi(image):
    """Add an Enhanced Vegetation Index (EVI) band.

    EVI = 2.5 * ((NIR - RED) /
          (NIR + 6 * RED - 7.5 * BLUE + 1)), clamped to
    [-2, 2].

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the EVI band added.
    """
    evi = image.expression(
        "2.5 * ((NIR - RED) / "
        "(NIR + 6 * RED - 7.5 * BLUE + 1))",
        {
            "NIR": image.select("SR_B4"),
            "RED": image.select("SR_B3"),
            "BLUE": image.select("SR_B1"),
        },
    ).rename("EVI")
    evi_clamped = evi.clamp(-2, 2)
    return image.addBands(evi_clamped)


def add_gndvi(image):
    """Add a Green NDVI (GNDVI) band to an image.

    GNDVI = (NIR - Green) / (NIR + Green).
    Gitelson and Merzlyak (1998).

    Parameters
    ----------
    image : ee.Image
        The input image.

    Returns
    -------
    ee.Image
        The image with the GNDVI band added.
    """
    gndvi = image.expression(
        "(NIR - Green) / (NIR + Green)",
        {
            "NIR": image.select("SR_B4"),
            "Green": image.select("SR_B2"),
        },
    ).rename("GNDVI")
    return image.addBands([gndvi])


def add_lai(image):
    """Add a Leaf Area Index (LAI) band to an image.

    LAI = 3.618 * EVI - 0.118, clamped to [0, 10].

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the LAI band added.
    """
    # Small constant to avoid division by zero.
    epsilon = 1e-10
    lai = image.expression(
        "3.618 * (EVI) - 0.118",
        {
            "EVI": image.expression(
                "2.5 * ((NIR - RED) / "
                "(NIR + 6 * RED - 7.5 * BLUE + 1 + "
                "epsilon))",
                {
                    "NIR": image.select("SR_B4"),
                    "RED": image.select("SR_B3"),
                    "BLUE": image.select("SR_B1"),
                    "epsilon": epsilon,
                },
            )
        },
    ).rename("LAI")
    lai = lai.clamp(0, 10)
    return image.addBands([lai])


def add_nbr(image):
    """Add a Normalized Burn Ratio (NBR) band to an image.

    NBR = (NIR - SWIR2) / (NIR + SWIR2)

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the NBR band added.
    """
    nbr = image.expression(
        "(NIR - SWIR2) / (NIR + SWIR2)",
        {
            "NIR": image.select("SR_B4"),
            "SWIR2": image.select("SR_B7"),
        },
    ).rename("NBR")
    return image.addBands([nbr])


def add_ndmi(image):
    """Add a Normalized Difference Moisture Index (NDMI).

    NDMI = (NIR - SWIR1) / (NIR + SWIR1)

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the NDMI band added.
    """
    ndmi = image.expression(
        "(NIR - SWIR1) / (NIR + SWIR1)",
        {
            "NIR": image.select("SR_B4"),
            "SWIR1": image.select("SR_B5"),
        },
    ).rename("NDMI")
    return image.addBands([ndmi])


def add_ndsi(image):
    """Add a Normalized Difference Snow Index (NDSI).

    NDSI = (Green - SWIR) / (Green + SWIR)

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the NDSI band added.
    """
    ndsi = image.expression(
        "(GREEN - SWIR) / (GREEN + SWIR)",
        {
            "GREEN": image.select("SR_B2"),
            "SWIR": image.select("SR_B5"),
        },
    ).rename("NDSI")
    return image.addBands([ndsi])


def add_ndvi(image):
    """Add a Normalized Difference Vegetation Index (NDVI).

    NDVI = (NIR - Red) / (NIR + Red)

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the NDVI band added.
    """
    ndvi = image.expression(
        "(NIR - Red) / (NIR + Red)",
        {
            "NIR": image.select("SR_B4"),
            "Red": image.select("SR_B3"),
        },
    ).rename("NDVI")
    return image.addBands([ndvi])


def add_ndwi(image):
    """Add a Normalized Difference Water Index (NDWI).

    NDWI = (Green - NIR) / (Green + NIR)

    Parameters
    ----------
    image : ee.Image
        The input image.

    Returns
    -------
    ee.Image
        The image with the NDWI band added.
    """
    ndwi = image.expression(
        "(Green - NIR) / (Green + NIR)",
        {
            "NIR": image.select("SR_B4"),
            "Green": image.select("SR_B2"),
        },
    ).rename("NDWI")
    return image.addBands([ndwi])


def add_ndrs(image, forest_types=None):
    """Add a Normalized Distance Red & SWIR (NDRS) band.

    Normalizes the DRS band within forest pixels and renames
    the band with a suffix based on forest class codes:

    - 210 : Coniferous  (_coni)
    - 220 : Broadleaf   (_deci)
    - 230 : Mixedwood   (_mixed)

    Forest data are sourced from
    https://gee-community-catalog.org/projects/ca_lc/.

    Parameters
    ----------
    image : ee.Image
        The image to process (must contain a 'DRS' band).
    forest_types : list of int, optional
        Forest type codes to include. Defaults to
        [210, 220, 230].

    Returns
    -------
    ee.Image
        The image with the renamed NDRS band added.
    """
    # Area of interest from the image geometry.
    aoi = image.geometry()

    # Extract the year from the image properties.
    year = ee.Number.parse(image.get("year"))

    # Define start and end dates based on the year.
    start_date = ee.Algorithms.If(
        year.gte(2019),
        ee.Date("2019-01-01"),
        ee.Date(year.format().cat("-01-01")),
    )
    end_date = ee.Algorithms.If(
        year.gte(2019),
        ee.Date("2019-12-31"),
        ee.Date(year.format().cat("-12-31")),
    )

    # Load landcover data for the specified period.
    lc_collection = forest_lc.lc_fn(start_date, end_date, aoi)
    landcover_image = ee.Image(lc_collection.first()).select(
        "forest_lc_class"
    )

    # Default to all three forest types.
    if forest_types is None:
        forest_types = [210, 220, 230]

    # Create a mask for the specified forest types.
    forest_mask = landcover_image.remap(
        forest_types,
        ee.List.repeat(1, len(forest_types)),
        0,
    )

    # Apply the forest mask to the DRS band.
    drs = image.select("DRS")
    masked_drs = drs.updateMask(forest_mask)

    # Calculate min and max of DRS for forest pixels.
    min_max = masked_drs.reduceRegion(
        reducer=ee.Reducer.minMax(),
        geometry=aoi.bounds(),
        scale=1000,
        maxPixels=1e10,
        bestEffort=True,
        tileScale=8,
    )
    drs_min = ee.Number(min_max.get("DRS_min"))
    drs_max = ee.Number(min_max.get("DRS_max"))

    # Clamp values to the [DRSmin, DRSmax] range.
    adjusted_drs = drs.clamp(drs_min, drs_max)

    # Calculate NDRS using the min and max values.
    ndrs = adjusted_drs.expression(
        "(DRS - DRSmin) / (DRSmax - DRSmin)",
        {
            "DRS": adjusted_drs,
            "DRSmin": drs_min,
            "DRSmax": drs_max,
        },
    ).rename("NDRS")

    # Determine the band-name suffix.
    if len(forest_types) == 1:
        if forest_types[0] == 210:
            suffix = "_coni"
        elif forest_types[0] == 220:
            suffix = "_deci"
        else:
            suffix = "_mixed"
    else:
        suffix = "_mixed"

    # Append the suffix to the NDRS band name.
    ndrs = ndrs.rename(
        ndrs.bandNames().map(
            lambda band_name: ee.String(band_name).cat(suffix)
        )
    )

    # Collapse combined suffixes to '_mixed'.
    renamed_bands = ndrs.bandNames().map(
        lambda band_name: ee.String(band_name).replace(
            "NDRS_coni_deci_mixed", "NDRS_mixed"
        )
    )
    ndrs = ndrs.rename(renamed_bands)

    return image.addBands(ndrs)


def add_savi(image):
    """Add a Soil Adjusted Vegetation Index (SAVI) band.

    SAVI = ((NIR - Red) / (NIR + Red + 0.428)) * 1.428

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the SAVI band added.
    """
    savi = image.expression(
        "((NIR - R) / (NIR + R + 0.428)) * (1.428)",
        {
            "NIR": image.select("SR_B4"),
            "R": image.select("SR_B3"),
        },
    ).rename("SAVI")
    return image.addBands([savi])


def add_si(image):
    """Add a Shadow Index (SI) band to an image.

    SI = (1 - Blue) * (1 - Green) * (1 - Red)

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the SI band added.
    """
    si = image.expression(
        "(1 - blue) * (1 - green) * (1 - red)",
        {
            "blue": image.select("SR_B1"),
            "green": image.select("SR_B2"),
            "red": image.select("SR_B3"),
        },
    ).rename("SI")
    return image.addBands([si])


def mask_cloud_snow(image):
    """Mask clouds and snow from a Landsat image.

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with clouds and snow masked.
    """
    qa = image.select("QA_PIXEL")
    clouds_bit_mask = 1 << 3
    cloud_shadow_bit_mask = 1 << 4
    snow_bit_mask = 1 << 5
    mask = (
        qa.bitwiseAnd(clouds_bit_mask)
        .eq(0)
        .And(qa.bitwiseAnd(cloud_shadow_bit_mask).eq(0))
        .And(qa.bitwiseAnd(snow_bit_mask).eq(0))
    )
    return image.updateMask(mask)


def mask_cloud(image):
    """Mask clouds from a Landsat image.

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with clouds masked.
    """
    qa = image.select("QA_PIXEL")
    clouds_bit_mask = 1 << 3
    cloud_shadow_bit_mask = 1 << 4
    mask = (
        qa.bitwiseAnd(clouds_bit_mask)
        .eq(0)
        .And(qa.bitwiseAnd(cloud_shadow_bit_mask).eq(0))
    )
    return image.updateMask(mask)


def mask_negative_surface_reflectance(image):
    """Mask negative surface reflectance values.

    Masks pixels where any of SR_B1 to SR_B7 is negative.

    Parameters
    ----------
    image : ee.Image
        The input image to process.

    Returns
    -------
    ee.Image
        The image with negative values masked.
    """
    bands_to_mask = [
        "SR_B1",
        "SR_B2",
        "SR_B3",
        "SR_B4",
        "SR_B5",
        "SR_B7",
    ]
    mask = (
        image.select(bands_to_mask)
        .reduce(ee.Reducer.min())
        .gte(0)
    )
    return image.updateMask(mask).copyProperties(
        image, image.propertyNames()
    )


def add_snow(image):
    """Add a snow band based on NDSI values.

    The snow band is true where NDSI > 0.4.

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the snow band added.
    """
    snow = (
        image.normalizedDifference(["SR_B2", "SR_B5"])
        .gt(0.4)
        .rename("snow")
    )
    return image.addBands([snow])


def mask_qa9(image):
    """Mask a Landsat image using QA_RADSAT Bit 9.

    Bit 9: 0 = pixel present, 1 = detector has no value.

    Parameters
    ----------
    image : ee.Image
        The Landsat image to mask.

    Returns
    -------
    ee.Image
        The masked image.
    """
    qa_band = image.select("QA_RADSAT")
    mask = qa_band.bitwiseAnd(1 << 9).eq(0)
    return image.updateMask(mask)


def mask_fill(image):
    """Mask fill pixels using QA_PIXEL Bit 0.

    Bit 0: 0 = valid pixel, 1 = fill pixel.

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with fill pixels masked.
    """
    qa = image.select("QA_PIXEL")
    mask = qa.bitwiseAnd(1).eq(0)
    return image.updateMask(mask)


def ndrs_stressed(image):
    """Add a binary mask of stressed forest pixels.

    Masks non-forest pixels, then thresholds the NDRS band.
    Pixels above the threshold are considered stressed.

    Parameters
    ----------
    image : ee.Image
        The image to process (must contain an NDRS band).

    Returns
    -------
    ee.Image
        The image with an 'NDRS_stressed' band added.
    """
    # Mask non-forest pixels, replace with zero for a
    # continuous raster.
    masked_image = masks.mask_by_landcover(image).unmask(0)

    # Threshold the NDRS band to identify stressed pixels.
    band = masked_image.select(BAND_NAME)
    binary_mask = band.gt(THRESHOLD).rename("NDRS_stressed")

    return image.addBands(binary_mask)


def apply_scale_factors(image):
    """Apply scaling factors to Landsat bands.

    Optical bands are scaled by 0.0000275 and offset by
    -0.2. The thermal band (ST_B6) is scaled by 0.00341802
    and offset by 149.0.

    Parameters
    ----------
    image : ee.Image
        The input Landsat image to be scaled.

    Returns
    -------
    ee.Image
        The image with optical and thermal bands scaled.
    """
    optical_bands = (
        image.select("SR_B.").multiply(0.0000275).add(-0.2)
    )
    thermal_band = (
        image.select("ST_B6").multiply(0.00341802).add(149.0)
    )
    return image.addBands(
        optical_bands, None, True
    ).addBands(thermal_band, None, True)


# End of script ----
