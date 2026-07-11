# ---
# title:   Get Annual Forest Tree Species
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Get annual species data from the High-resolution Annual
#   Forest Species Maps for Canada's Forested Ecosystems.
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


def species_fn(start_date, end_date, aoi=None):
    """Get annual forest species images for a date range.

    Filters the CA_FOREST_SPECIES collection to the given
    date range, optionally clips to an area of interest,
    and renames the species band.

    Args:
        start_date (str): Start date string (YYYY-MM-DD).
        end_date (str): End date string (YYYY-MM-DD).
        aoi (ee.Geometry): Optional area of interest. When
            provided, the collection is filtered and each
            image is clipped to this geometry.

    Returns:
        ee.ImageCollection: Species images with a single
        'forest_species_class' band, clipped to the AOI
        when provided.
    """
    # Get the species collection for the date range
    species_collection = ee.ImageCollection(
        "projects/sat-io/open-datasets/CA_FOREST_SPECIES"
    ).filterDate(start_date, end_date)

    # Apply area of interest (AOI) filter if provided
    if aoi is not None:
        species_collection = species_collection.filterBounds(aoi)

    def _process(image):
        # Select and rename the band, clip to AOI if provided
        img = image.select("species").rename(
            "forest_species_class"
        )
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

    species_collection = species_collection.map(_process)

    return species_collection
