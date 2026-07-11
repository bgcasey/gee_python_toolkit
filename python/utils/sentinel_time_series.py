# ---
# title:   Get a Time Series of Sentinel-2 Images
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Processes Sentinel-2 imagery, calculates selected
#   vegetation indices, and merges the results into a single
#   image collection for a period and area of interest
#   (AOI).
#
#   Steps:
#   1. Retrieve the Sentinel-2 collection for the date range
#      and AOI.
#   2. Apply cloud masking to the images.
#   3. Calculate the selected indices for each image.
#   4. Merge results into a median composite per date range.
# ---

import ee

from utils import sentinel_indices_and_masks as indices

# Mapping of index code to the function that adds the band.
# NDMI maps to add_ndwi, matching the original source.
INDEX_FUNCTIONS = {
    "CRE": indices.add_cre,
    "DRS": indices.add_drs,
    "DSWI": indices.add_dswi,
    "EVI": indices.add_evi,
    "GNDVI": indices.add_gndvi,
    "LAI": indices.add_lai,
    "NBR": indices.add_nbr,
    "NDMI": indices.add_ndwi,
    "NDRE1": indices.add_ndre1,
    "NDRE2": indices.add_ndre2,
    "NDRE3": indices.add_ndre3,
    "NDVI": indices.add_ndvi,
    "NDRS": indices.add_ndrs,
    "NDWI": indices.add_ndwi,
    "RDI": indices.add_rdi,
}


def s2_fn(
    dates, interval, interval_type, aoi, selected_indices
):
    """Process Sentinel-2 images and merge into a series.

    Parameters
    ----------
    dates : list of str
        Date strings for image collection time ranges.
    interval : int
        Interval units to advance from each date.
    interval_type : str
        Type of interval ('days', 'weeks', 'months',
        'years').
    aoi : ee.Geometry
        Area of interest.
    selected_indices : list of str
        Indices to calculate (e.g., ['NDVI']).

    Returns
    -------
    ee.ImageCollection
        Processed images clipped to the AOI.
    """

    def s2_ts(d1):
        """Process images for a single date.

        Parameters
        ----------
        d1 : str
            Start date string for the image collection.

        Returns
        -------
        ee.Image
            Median image with selected indices.
        """
        start = ee.Date(d1)
        min_start_date = ee.Date("2019-06-01")

        # Clamp the start date to on or after 2019-06-01.
        start = ee.Date(
            ee.Algorithms.If(
                start.millis().lt(min_start_date.millis()),
                min_start_date,
                start,
            )
        )

        end = start.advance(interval, interval_type)

        # Get the Sentinel-2 collection for the date range.
        s2_collection = (
            ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
            .filterBounds(aoi)
            .filterDate(start, end)
            # Pre-filter to less cloudy granules.
            .filter(
                ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20)
            )
            .map(indices.mask_s2_clouds)
        )

        # Apply selected indices to the collection.
        for index in selected_indices:
            fn = INDEX_FUNCTIONS.get(index)
            if fn is not None:
                s2_collection = s2_collection.map(fn)

        # Median composite of the raw bands.
        raw_bands = s2_collection.select(
            [
                "B1", "B2", "B3", "B4", "B5", "B6", "B7",
                "B8", "B8A", "B9", "B11", "B12",
            ]
        ).median()

        # Median composite of the calculated indices.
        indices_composite = s2_collection.select(
            selected_indices
        ).median()

        # Combine raw bands and calculated indices.
        return indices_composite.addBands(raw_bands).set(
            {
                "start_date": start.format("YYYY-MM-dd"),
                "end_date": end.format("YYYY-MM-dd"),
                "month": start.get("month"),
                "year": start.get("year"),
            }
        )

    # Map over dates, clip to AOI, and return a collection.
    s2 = ee.ImageCollection.fromImages(
        [s2_ts(d).clip(aoi) for d in dates]
    )

    return s2


# End of script ----
