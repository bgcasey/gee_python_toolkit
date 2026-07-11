# ---
# title:   Get Canopy Height
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Load and combine global canopy height and its standard
#   deviation into a single image for a given AOI.
#
#   Data citation:
#   Lang, N., Jetz, W., Schindler, K., Wegner, J.D. (2022).
#   A high-resolution canopy height model of the Earth.
#   arXiv preprint arXiv:2204.08322.
# ---

import ee


def get_canopy_data(aoi):
    """Combine canopy height and its standard deviation.

    Args:
        aoi (ee.Geometry): Area of interest.

    Returns:
        ee.Image: Two-band image with 'canopy_height' and
        'canopy_standard_deviation', clipped to the AOI.
    """
    # Load canopy height, rename, and clip to AOI
    canopy_height = (
        ee.Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1")
        .rename("canopy_height")
        .clip(aoi)
    )

    # Load canopy height standard deviation, rename, and clip
    canopy_sd = (
        ee.Image("users/nlang/ETH_GlobalCanopyHeightSD_2020_10m_v1")
        .rename("canopy_standard_deviation")
        .clip(aoi)
    )

    # Combine canopy height and standard deviation bands
    canopy = canopy_height.addBands([canopy_sd])

    return canopy
