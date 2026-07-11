# ---
# title:   Get a Time Series of Landsat Images
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Processes Landsat imagery (Landsat 5, 7, 8, and 9),
#   harmonizes spectral reflectance across sensors,
#   calculates selected vegetation indices, and merges the
#   results into a single image collection. Landsat 5, 8,
#   and 9 are prioritized over Landsat 7 due to the latter's
#   Scan Line Corrector (SLC) failure.
#
#   Steps:
#   1. Retrieve and harmonize Landsat Surface Reflectance
#      (SR) collections for a period and area of interest
#      (AOI).
#   2. Combine harmonized Landsat 5, 7, 8, and 9 collections,
#      prioritizing 5, 8, and 9 over 7.
#   3. Composite selected vegetation indices and merge them
#      into a single image collection.
# ---

import ee

from utils import landsat_indices_and_masks as li

# Mapping of index code to the function that adds the band.
INDEX_FUNCTIONS = {
    "BSI": li.add_bsi,
    "DRS": li.add_drs,
    "DSWI": li.add_dswi,
    "EVI": li.add_evi,
    "GNDVI": li.add_gndvi,
    "LAI": li.add_lai,
    "NBR": li.add_nbr,
    "NDMI": li.add_ndmi,
    "NDSI": li.add_ndsi,
    "NDVI": li.add_ndvi,
    "NDWI": li.add_ndwi,
    "SAVI": li.add_savi,
    "SI": li.add_si,
}


def harmonize_oli_to_etm(image):
    """Harmonize Landsat 8/9 (OLI) to Landsat 7 (ETM+).

    Uses reduced major axis (RMA) regression coefficients.

    Citation: Roy, D.P., et al. (2016). Characterization of
    Landsat-7 to Landsat-8 reflective wavelength and
    normalized difference vegetation index continuity.
    Remote Sensing of Environment, 185, 57-70.
    doi:10.1016/j.rse.2015.12.024; Table 2.

    Parameters
    ----------
    image : ee.Image
        The input Landsat 8 or 9 image.

    Returns
    -------
    ee.Image
        The harmonized image.
    """
    # Slopes and intercepts for each band.
    slopes = ee.Image.constant(
        [0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949]
    )
    itcp = ee.Image.constant(
        [-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029]
    )

    # Apply the harmonization transformation.
    harmonized = (
        image.select(
            ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6",
             "SR_B7"],
            ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5",
             "SR_B7"],
        )
        .resample("bicubic")
        .subtract(itcp.multiply(10000))
        .divide(slopes)
        .set(
            "system:time_start",
            image.get("system:time_start"),
        )
    )

    # Preserve the QA_PIXEL band.
    qa_pixel = image.select("QA_PIXEL")
    return harmonized.addBands(qa_pixel, None, True)


def get_harmonized_ls_collection(
    start_date, end_date, sensor, aoi
):
    """Retrieve and harmonize a Landsat SR collection.

    Parameters
    ----------
    start_date : str
        The start date for the collection.
    end_date : str
        The end date for the collection.
    sensor : str
        The Landsat sensor code (e.g., 'LC08').
    aoi : ee.Geometry
        The area of interest.

    Returns
    -------
    ee.ImageCollection
        The harmonized image collection.
    """
    ls_collection = (
        ee.ImageCollection(
            "LANDSAT/" + sensor + "/C02/T1_L2"
        )
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .map(li.mask_cloud_snow)
    )

    # Apply harmonization to Landsat 8 or 9 images.
    if sensor in ("LC08", "LC09"):
        ls_collection = ls_collection.map(
            harmonize_oli_to_etm
        )

    # Apply scaling to all images.
    def scale_reflectance(image):
        return image.multiply(0.0000275).add(-0.2)

    ls_collection = ls_collection.map(scale_reflectance)

    # Apply the negative-value mask.
    ls_collection = ls_collection.map(
        li.mask_negative_surface_reflectance
    )

    # Select relevant bands and re-add the QA_PIXEL band.
    def select_bands(img):
        qa_pixel = img.select("QA_PIXEL")
        return img.select(
            ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5",
             "SR_B7"]
        ).addBands(qa_pixel)

    return ls_collection.map(select_bands)


def get_combined_harmonized_collection(
    start_date, end_date, aoi
):
    """Combine harmonized Landsat 5, 7, 8, and 9 collections.

    Prioritizes Landsat 5, 8, and 9 over 7 due to banding
    from the failure of Landsat 7's Scan Line Corrector
    (SLC).

    Parameters
    ----------
    start_date : str
        The start date for the collection.
    end_date : str
        The end date for the collection.
    aoi : ee.Geometry
        The area of interest.

    Returns
    -------
    ee.ImageCollection
        The combined harmonized collection.
    """
    # Retrieve harmonized collections for each sensor.
    lt5 = get_harmonized_ls_collection(
        start_date, end_date, "LT05", aoi
    )
    le7 = get_harmonized_ls_collection(
        start_date, end_date, "LE07", aoi
    )
    lc8 = get_harmonized_ls_collection(
        start_date, end_date, "LC08", aoi
    )
    lc9 = get_harmonized_ls_collection(
        start_date, end_date, "LC09", aoi
    )

    # Collection sizes for priority handling.
    lt5_size = lt5.size()
    lc8_size = lc8.size()
    lc9_size = lc9.size()

    # Combine collections based on availability and priority.
    combined_collection = ee.ImageCollection(
        ee.Algorithms.If(
            lt5_size.gt(0).And(
                lc8_size.gt(0).Or(lc9_size.gt(0))
            ),
            lt5.merge(lc8).merge(lc9),
            ee.Algorithms.If(
                lt5_size.gt(0),
                lt5,
                ee.Algorithms.If(
                    lc8_size.gt(0).Or(lc9_size.gt(0)),
                    lc8.merge(lc9),
                    le7,
                ),
            ),
        )
    )

    return combined_collection


def ls_fn(
    dates,
    interval,
    interval_type,
    aoi,
    selected_indices,
    statistic,
):
    """Process Landsat images and merge them into a series.

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
    statistic : str
        Statistic to apply ('mean', 'median', 'max', etc.).

    Returns
    -------
    ee.ImageCollection
        Processed images clipped to the AOI.
    """

    def ls_ts(d1):
        """Process images for a single date.

        Parameters
        ----------
        d1 : str
            Start date string for the image collection.

        Returns
        -------
        ee.Image
            Image reduced by the specified statistic.
        """
        start = ee.Date(d1)
        end = start.advance(interval, interval_type)

        # Get the combined Landsat collection for the range.
        combined_collection = (
            get_combined_harmonized_collection(
                start, end, aoi
            )
        )

        # Apply selected indices to the collection.
        for index in selected_indices:
            fn = INDEX_FUNCTIONS.get(index)
            if fn is not None:
                combined_collection = combined_collection.map(
                    fn
                )

        # Dynamically apply the specified statistic.
        reducer = getattr(ee.Reducer, statistic)()
        reduced_image = combined_collection.reduce(reducer)

        # Rename bands to remove the "_statistic" suffix.
        renamed_bands = reduced_image.bandNames().map(
            lambda band_name: ee.String(band_name).replace(
                "_" + statistic + "$", ""
            )
        )
        reduced_image = reduced_image.rename(renamed_bands)

        # Mask for missing pixels (where mask is 0).
        missing_pixels_mask = reduced_image.mask().Not()

        # Gaussian-spline interpolation of missing pixels.
        interpolated_image = reduced_image.focal_mean(
            kernel=ee.Kernel.gaussian(
                radius=3, sigma=1, units="pixels"
            ),
            iterations=1,
        ).updateMask(missing_pixels_mask)

        # Combine original and interpolated images.
        combined_image = reduced_image.unmask(
            interpolated_image
        )

        # Set metadata for the reduced image.
        return combined_image.set(
            {
                "start_date": start.format("YYYY-MM-dd"),
                "end_date": end.format("YYYY-MM-dd"),
                "month": start.get("month"),
                "year": start.get("year"),
            }
        )

    # Map over dates, clip to AOI, and return a collection.
    ls = ee.ImageCollection.fromImages(
        [ls_ts(d).clip(aoi) for d in dates]
    )

    return ls


# End of script ----
