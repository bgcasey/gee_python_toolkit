# ---
# title:   Land Cover Proportion Functions
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Calculate proportions of remapped land cover classes
#   annually across a given AOI, per-pixel or within a
#   neighborhood kernel.
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

# Remapped land cover class names (in target class order)
CLASS_NAMES = [
    "unclassified",
    "water",
    "snow_ice",
    "rock_rubble",
    "exposed_barren_land",
    "bryoids",
    "shrubs",
    "wetland",
    "wetland_treed",
    "herbs",
    "coniferous",
    "broadleaf",
    "mixedwood",
]

# Original class codes and their remapped target codes
FROM = [0, 20, 31, 32, 33, 40, 50, 80, 81, 100, 210, 220, 230]
TO = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]


def landcover_proportion(dates, interval, aoi):
    """Calculate per-pixel land cover class proportions.

    Args:
        dates (list): List of start date strings for each
            interval in the time series.
        interval (int): Interval length in months.
        aoi (ee.Geometry): Area of interest.

    Returns:
        ee.ImageCollection: Annual images of per-pixel
        proportions of remapped land cover classes.
    """
    ca_lc = ee.ImageCollection(
        "projects/sat-io/open-datasets/CA_FOREST_LC_VLCE2"
    )

    def _landcover_ts(d1):
        start = ee.Date(d1)
        end = ee.Date(d1).advance(interval, "month")
        date = ee.Date(d1)

        # Filter collection for the specific date range
        lc_image = (
            ee.Image(ca_lc.filterDate(start, end).first())
            .remap(FROM, TO)
            .clip(aoi)
        )

        # Calculate per-pixel proportions (0/1 masks)
        def _calculate_proportions(image):
            proportions = [
                image.eq(value).rename("Proportion_" + str(value))
                for value in TO
            ]
            return ee.Image(proportions)

        lc_proportions = (
            _calculate_proportions(lc_image).unmask(0).clip(aoi)
        )

        # Rename bands with readable class names
        renamed = ee.Image(
            [
                lc_proportions.select(index).rename(
                    CLASS_NAMES[index]
                )
                for index in range(len(TO))
            ]
        )

        return renamed.set("year", date.get("year"))

    # Generate annual land cover proportions
    return ee.ImageCollection(
        [_landcover_ts(d) for d in dates]
    ).map(lambda img: img.clip(aoi))


def landcover_proportion_focal(dates, interval, kernel_size, aoi):
    """Calculate focal land cover class proportions.

    Args:
        dates (list): List of start date strings for each
            interval in the time series.
        interval (int): Interval length in months.
        kernel_size (float): Kernel radius in meters.
        aoi (ee.Geometry): Area of interest.

    Returns:
        ee.ImageCollection: Annual images of focal
        proportions of remapped land cover classes, with the
        kernel size appended to each band name.
    """
    ca_lc = ee.ImageCollection(
        "projects/sat-io/open-datasets/CA_FOREST_LC_VLCE2"
    )

    def _landcover_ts(d1):
        start = ee.Date(d1)
        end = ee.Date(d1).advance(interval, "month")
        date = ee.Date(d1)

        # Filter collection for the specific date range
        lc_image = (
            ee.Image(ca_lc.filterDate(start, end).first())
            .remap(FROM, TO)
            .clip(aoi)
        )

        # Define kernel radius in meters and pixels
        radius_in_meters = kernel_size
        projection = lc_image.projection()
        radius_in_pixels = (
            ee.Number(radius_in_meters)
            .divide(projection.nominalScale())
            .round()
        )
        kernel = ee.Kernel.circle(radius_in_pixels, "pixels")

        # Calculate proportions within the kernel
        def _calculate_proportions(image):
            proportions = []
            for value in TO:
                class_count = image.updateMask(
                    image.eq(value)
                ).reduce(ee.Reducer.count())
                total_count = image.reduce(ee.Reducer.count())
                proportions.append(
                    class_count.divide(total_count).rename(
                        "Proportion_" + str(value)
                    )
                )
            return ee.Image(proportions)

        lc_proportions = (
            _calculate_proportions(
                lc_image.neighborhoodToBands(kernel)
            )
            .unmask(0)
            .clip(aoi)
        )

        # Rename bands with readable class names
        renamed = ee.Image(
            [
                lc_proportions.select(index).rename(
                    CLASS_NAMES[index]
                )
                for index in range(len(TO))
            ]
        )

        # Add kernel size suffix to band names
        band_names = renamed.bandNames()

        def _append_kernel_size(band_name):
            return (
                ee.String(band_name)
                .cat("_")
                .cat(str(radius_in_meters))
            )

        renamed = renamed.rename(
            band_names.map(_append_kernel_size)
        )

        return renamed.set("year", date.get("year"))

    # Generate focal land cover proportions
    return ee.ImageCollection(
        [_landcover_ts(d) for d in dates]
    ).map(lambda img: img.clip(aoi))
