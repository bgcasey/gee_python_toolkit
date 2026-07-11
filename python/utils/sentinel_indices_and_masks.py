# ---
# title:   Sentinel-2 Indices and Masks Functions
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Functions to calculate spectral indices and apply masks
#   to a time series of Sentinel-2 images. Indices include
#   vegetation, moisture, and stress-related measures. Masks
#   handle cloud, snow, and QA filtering.
# ---

import ee

from utils import annual_forest_land_cover as forest_lc
from utils import masks

# Band name and threshold used to classify stressed forest
# pixels from the NDRS index.
BAND_NAME = "NDRS"
THRESHOLD = 0.5


def add_cre(image):
    """Add a Red Edge Chlorophyll Index (CRE) band.

    CRE = (RedEdge3 / RedEdge1) - 1.
    Gitelson et al. (2003).

    Parameters
    ----------
    image : ee.Image
        The input image.

    Returns
    -------
    ee.Image
        The image with the CRE band added.
    """
    cre = image.expression(
        "(RedEdge3 / RedEdge1) - 1",
        {
            "RedEdge1": image.select("B5"),
            "RedEdge3": image.select("B7"),
        },
    ).rename("CRE")
    return image.addBands([cre])


def add_dswi(image):
    """Add a Disease Stress Water Index (DSWI) band.

    DSWI = (NIR + Green) / (Red + SWIR)

    Parameters
    ----------
    image : ee.Image
        The input image.

    Returns
    -------
    ee.Image
        The image with the DSWI band added.
    """
    dswi = image.expression(
        "(NIR + Green) / (Red + SWIR)",
        {
            "NIR": image.select("B8"),
            "Green": image.select("B3"),
            "Red": image.select("B4"),
            "SWIR": image.select("B11"),
        },
    ).rename("DSWI")
    return image.addBands([dswi])


def add_drs(image):
    """Add a Distance Red & SWIR (DRS) band to an image.

    DRS = sqrt((RED^2) + (SWIR^2))

    Parameters
    ----------
    image : ee.Image
        The input image.

    Returns
    -------
    ee.Image
        The image with the DRS band added.
    """
    drs = image.expression(
        "sqrt(((RED) * (RED)) + ((SWIR) * (SWIR)))",
        {
            "SWIR": image.select("B11"),
            "RED": image.select("B4"),
        },
    ).rename("DRS")
    return image.addBands([drs]).copyProperties(
        image, ["system:time_start"]
    )


def add_evi(image):
    """Add an Enhanced Vegetation Index (EVI) band.

    EVI = 2.5 * ((NIR - RED) /
          (NIR + 6 * RED - 7.5 * BLUE + 1))

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
            "NIR": image.select("B8"),
            "RED": image.select("B4"),
            "BLUE": image.select("B2"),
        },
    ).rename("EVI")
    return image.addBands([evi])


def add_gndvi(image):
    """Add a Green NDVI (GNDVI) band to an image.

    GNDVI = (B8 - B3) / (B8 + B3).
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
    gndvi = image.normalizedDifference(["B8", "B3"]).rename(
        "GNDVI"
    )
    return image.addBands([gndvi])


def add_lai(image):
    """Add a Leaf Area Index (LAI) band to an image.

    LAI = 3.618 * EVI - 0.118

    Parameters
    ----------
    image : ee.Image
        The image to process.

    Returns
    -------
    ee.Image
        The image with the LAI band added.
    """
    lai = image.expression(
        "3.618 * (EVI) - 0.118",
        {
            "EVI": image.expression(
                "2.5 * ((NIR - RED) / "
                "(NIR + 6 * RED - 7.5 * BLUE + 1))",
                {
                    "NIR": image.select("B8"),
                    "RED": image.select("B4"),
                    "BLUE": image.select("B2"),
                },
            )
        },
    ).rename("LAI")
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
            "NIR": image.select("B8"),
            "SWIR2": image.select("B12"),
        },
    ).rename("NBR")
    return image.addBands([nbr])


def add_ndre1(image):
    """Add a Normalized Difference Red-edge 1 (NDRE1) band.

    NDRE1 = (RedEdge2 - RedEdge1) / (RedEdge2 + RedEdge1)

    Parameters
    ----------
    image : ee.Image
        The input image.

    Returns
    -------
    ee.Image
        The image with the NDRE1 band added.
    """
    ndre1 = image.expression(
        "(RedEdge2 - RedEdge1) / (RedEdge2 + RedEdge1)",
        {
            "RedEdge2": image.select("B6"),
            "RedEdge1": image.select("B5"),
        },
    ).rename("NDRE1")
    return image.addBands([ndre1])


def add_ndre2(image):
    """Add a Normalized Difference Red-edge 2 (NDRE2) band.

    NDRE2 = (RedEdge3 - RedEdge1) / (RedEdge3 + RedEdge1)

    Parameters
    ----------
    image : ee.Image
        The input image.

    Returns
    -------
    ee.Image
        The image with the NDRE2 band added.
    """
    ndre2 = image.expression(
        "(RedEdge3 - RedEdge1) / (RedEdge3 + RedEdge1)",
        {
            "RedEdge3": image.select("B7"),
            "RedEdge1": image.select("B5"),
        },
    ).rename("NDRE2")
    return image.addBands([ndre2])


def add_ndre3(image):
    """Add a Normalized Difference Red-edge 3 (NDRE3) band.

    NDRE3 = (RedEdge4 - RedEdge3) / (RedEdge4 + RedEdge3).
    Checks for required bands (B8A, B7) before applying the
    calculation; returns the original image if they are
    missing.

    Parameters
    ----------
    image : ee.Image
        The input image.

    Returns
    -------
    ee.Image
        The image with the NDRE3 band added if bands exist.
    """
    bands = ["B8A", "B7"]
    has_bands = ee.List(bands).map(
        lambda b: image.bandNames().contains(b)
    )
    has_all = ee.List(has_bands).reduce(ee.Reducer.min())

    return ee.Image(
        ee.Algorithms.If(
            ee.Number(has_all).eq(1),
            image.addBands(
                image.expression(
                    "(RedEdge4 - RedEdge3) / "
                    "(RedEdge4 + RedEdge3)",
                    {
                        "RedEdge4": image.select("B8A"),
                        "RedEdge3": image.select("B7"),
                    },
                ).rename("NDRE3")
            ),
            image,
        )
    )


def add_ndvi(image):
    """Add a Normalized Difference Vegetation Index (NDVI).

    NDVI = (B8 - B4) / (B8 + B4)

    Parameters
    ----------
    image : ee.Image
        The input image.

    Returns
    -------
    ee.Image
        The image with the NDVI band added.
    """
    ndvi = image.normalizedDifference(["B8", "B4"]).rename(
        "NDVI"
    )
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
            "NIR": image.select("B8"),
            "Green": image.select("B3"),
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


def add_rdi(image):
    """Add a Ratio Drought Index (RDI) band to an image.

    RDI = SWIR2 / RedEdge4. Checks for required bands
    (B12, B8A) before applying the calculation; returns the
    original image if they are missing.

    Parameters
    ----------
    image : ee.Image
        The input image.

    Returns
    -------
    ee.Image
        The image with the RDI band added if bands exist.
    """
    bands = ["B12", "B8A"]
    has_bands = ee.List(bands).map(
        lambda b: image.bandNames().contains(b)
    )
    has_all = ee.List(has_bands).reduce(ee.Reducer.min())

    return ee.Image(
        ee.Algorithms.If(
            ee.Number(has_all).eq(1),
            image.addBands(
                image.expression(
                    "SWIR2 / RedEdge4",
                    {
                        "SWIR2": image.select("B12"),
                        "RedEdge4": image.select("B8A"),
                    },
                ).rename("RDI")
            ),
            image,
        )
    )


def create_binary_mask(image):
    """Add a binary mask of stressed forest pixels.

    Masks non-forest pixels, then thresholds the NDRS band.
    Pixels above the threshold are considered stressed.

    Parameters
    ----------
    image : ee.Image
        The input image (must contain an NDRS band).

    Returns
    -------
    ee.Image
        The image with an 'NDRS_stressed' band added.
    """
    masked_image = masks.mask_by_landcover(image).unmask(0)
    band = masked_image.select(BAND_NAME)
    binary_mask = band.gt(THRESHOLD).rename("NDRS_stressed")
    return image.addBands(binary_mask)


def mask_s2_clouds(image):
    """Mask clouds using the Sentinel-2 QA60 band.

    Bits 10 and 11 flag clouds and cirrus, respectively.
    Reflectance is scaled by dividing by 10000.

    Parameters
    ----------
    image : ee.Image
        The Sentinel-2 image to process.

    Returns
    -------
    ee.Image
        The cloud-masked, scaled Sentinel-2 image.
    """
    qa = image.select("QA60")
    cloud_bit_mask = 1 << 10
    cirrus_bit_mask = 1 << 11
    mask = (
        qa.bitwiseAnd(cloud_bit_mask)
        .eq(0)
        .And(qa.bitwiseAnd(cirrus_bit_mask).eq(0))
    )
    return image.updateMask(mask).divide(10000)


# End of script ----
