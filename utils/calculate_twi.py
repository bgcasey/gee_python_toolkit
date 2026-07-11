# ---
# title:   Calculate TWI from FABDEM
# author:  bgcasey
# created: 2026-07-10
# notes:
#   Calculate the Topographic Wetness Index (TWI) for a
#   given AOI: ln(a / tan(b)), where a is upslope drainage
#   area and b is slope.
#
#   FABDEM is a bare-earth DEM with no flow-accumulation
#   band, and Earth Engine has no native flow-accumulation
#   algorithm. This function therefore uses a hybrid:
#   slope from FABDEM (30 m, forests and buildings
#   removed) and upslope area from MERIT Hydro 'upa'
#   (~90 m, resampled on the fly).
#
#   Data citations:
#   Hawker, L., et al. (2022). A 30 m global map of
#   elevation with forests and buildings removed.
#   Environmental Research Letters, 17(2), 024016.
#   doi:10.1088/1748-9326/ac4d4f
#
#   Yamazaki, D., et al. (2019). MERIT Hydro: A
#   high-resolution global hydrography map based on latest
#   topography datasets. Water Resources Research, 55,
#   5053-5073. doi:10.1029/2019WR024873
# ---

import math

import ee


def calculate_twi_fabdem(aoi):
    """Calculate TWI using FABDEM slope and MERIT Hydro upa.

    TWI: ln(a / tan(b)), where a is upslope drainage area
    (m^2) and b is slope (radians).

    Args:
        aoi (ee.Geometry): Area of interest.

    Returns:
        ee.Image: Single-band 'twi' image.
    """
    # Load FABDEM, mosaic, and set a default projection so
    # terrain algorithms have a fixed 30 m scale
    elevation = (
        ee.ImageCollection("projects/sat-io/open-datasets/FABDEM")
        .mosaic()
        .setDefaultProjection("EPSG:3857", None, 30)
        .clip(aoi)
    )

    # Calculate slope from FABDEM elevation
    slope = ee.Terrain.slope(elevation)

    # Load upslope area from MERIT Hydro and convert
    # km^2 to m^2
    upslope_area = (
        ee.Image("MERIT/Hydro/v1_0_1")
        .select("upa")
        .clip(aoi)
        .multiply(1e6)
        .rename("upslope_area")
    )

    # Convert slope from degrees to radians
    slope_rad = slope.multiply(math.pi / 180).rename("slope_rad")

    # Floor tan(b) at a small value so flat areas (slope 0)
    # are not masked by division by zero
    tan_b = slope_rad.tan().max(0.001)

    # Calculate TWI: ln(a / tan(b))
    twi = upslope_area.divide(tan_b).log().rename("twi")

    return twi
