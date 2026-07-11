# ---
# title:   Annual Forest Land Cover
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Get annual land cover data from the High-resolution
#   Annual Forest Land Cover Maps for Canada's Forested
#   Ecosystems (1984-2019).
#
#   Data citation:
#   Hermosilla, T., Wulder, M.A., White, J.C., Coops, N.C.,
#   2022. Land cover classification in an era of big and
#   open data: Optimizing localized implementation and
#   training data selection to improve mapping outcomes.
#   Remote Sensing of Environment. No. 112780.
#   doi:10.1016/j.rse.2022.112780
# ---

import ee


def lc_fn(start_date, end_date, aoi=None):
    """Get annual forest land cover images for a date range.

    Filters the CA_FOREST_LC_VLCE2 collection to the given
    date range, optionally clips to an area of interest,
    and renames the land cover band.

    Args:
        start_date (str): Start date string (YYYY-MM-DD).
        end_date (str): End date string (YYYY-MM-DD).
        aoi (ee.Geometry): Optional area of interest. When
            provided, the collection is filtered and each
            image is clipped to this geometry.

    Returns:
        ee.ImageCollection: Land cover images with a single
        'forest_lc_class' band, clipped to the AOI when
        provided.
    """
    # Get the land cover collection for the date range
    lc_collection = ee.ImageCollection(
        "projects/sat-io/open-datasets/CA_FOREST_LC_VLCE2"
    ).filterDate(start_date, end_date)

    # Apply area of interest (AOI) filter if provided
    if aoi is not None:
        lc_collection = lc_collection.filterBounds(aoi)

    def _process(image):
        # Select and rename the band, clip to AOI if provided
        img = image.select("b1").rename("forest_lc_class")
        if aoi is not None:
            img = img.clip(aoi)
        return img.set(
            {
                "start_date": ee.Date(start_date).format(
                    "YYYY-MM-dd"
                ),
                "end_date": ee.Date(end_date).format(
                    "YYYY-MM-dd"
                ),
                "year": ee.Date(
                    image.get("system:time_start")
                ).get("year"),
            }
        )

    lc_collection = lc_collection.map(_process)

    return lc_collection
