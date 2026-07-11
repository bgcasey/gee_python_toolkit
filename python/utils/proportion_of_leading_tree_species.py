# ---
# title:   Leading Tree Species Proportion Functions
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Calculate proportions of leading tree species annually
#   across a given AOI, per-pixel or within a neighborhood
#   kernel.
#
#   Data citation:
#   Hermosilla, T., Wulder, M.A., White, J.C., Coops, N.C.,
#   Bater, C.W., Hobart, G.W., 2024. Characterizing
#   long-term tree species dynamics in Canada's forested
#   ecosystems using annual time series remote sensing
#   data. Forest Ecology and Management 572, 122313.
#   doi:10.1016/j.foreco.2024.122313
# ---

import ee

# Tree species class IDs and their snake_case names
SPECIES_VALUES = list(range(38))
SPECIES_NAMES = [
    "non_tree",
    "amabilis_fir",
    "balsam_fir",
    "subalpine_fir",
    "bigleaf_maple",
    "red_maple",
    "sugar_maple",
    "gray_alder",
    "red_alder",
    "yellow_birch",
    "white_birch",
    "yellow_cedar",
    "black_ash",
    "tamarack",
    "western_larch",
    "norway_spruce",
    "engelmann_spruce",
    "white_spruce",
    "black_spruce",
    "red_spruce",
    "sitka_spruce",
    "whitebark_pine",
    "jack_pine",
    "lodgepole_pine",
    "ponderosa_pine",
    "red_pine",
    "eastern_white_pine",
    "balsam_poplar",
    "largetooth_aspen",
    "trembling_aspen",
    "douglas_fir",
    "red_oak",
    "eastern_white_cedar",
    "western_redcedar",
    "eastern_hemlock",
    "western_hemlock",
    "mountain_hemlock",
    "white_elm",
]


def tree_species_proportion_focal(dates, interval, kernel_size, aoi):
    """Calculate focal proportions of leading tree species.

    Args:
        dates (list): List of start date strings for each
            interval in the time series.
        interval (int): Interval length in months.
        kernel_size (float): Kernel radius in meters.
        aoi (ee.Geometry): Area of interest.

    Returns:
        ee.ImageCollection: Annual images of focal
        proportions of leading tree species, with the kernel
        size appended to each band name.
    """

    def _species_ts(d1):
        start = ee.Date(d1)
        end = ee.Date(d1).advance(interval, "month")
        date = ee.Date(d1)

        # Load tree species image for the date range
        species = (
            ee.ImageCollection(
                "projects/sat-io/open-datasets/"
                "CA_FOREST/SPECIES-1984-2022"
            )
            .filterDate(start, end)
            .first()
            .clip(aoi)
        )

        # Define the kernel radius in meters and pixels
        radius_in_meters = kernel_size
        projection = species.projection()
        radius_in_pixels = (
            ee.Number(radius_in_meters)
            .divide(projection.nominalScale())
            .round()
        )
        kernel = ee.Kernel.circle(radius_in_pixels, "pixels")

        # Calculate proportions within the kernel
        def _calculate_species_proportions(image):
            proportions = []
            for value in SPECIES_VALUES:
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

        species_proportions = (
            _calculate_species_proportions(
                species.neighborhoodToBands(kernel)
            )
            .unmask(0)
            .clip(aoi)
        )

        # Rename bands with readable snake_case names
        renamed = ee.Image(
            [
                species_proportions.select(index).rename(
                    SPECIES_NAMES[index]
                )
                for index in range(len(SPECIES_VALUES))
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

    # Generate annual focal tree species proportions
    return ee.ImageCollection(
        [_species_ts(d) for d in dates]
    ).map(lambda img: img.clip(aoi))


def tree_species_proportion(dates, interval, aoi):
    """Calculate per-pixel proportions of tree species.

    Args:
        dates (list): List of start date strings for each
            interval in the time series.
        interval (int): Interval length in months.
        aoi (ee.Geometry): Area of interest.

    Returns:
        ee.ImageCollection: Annual images of per-pixel
        proportions of leading tree species.
    """

    def _species_ts(d1):
        start = ee.Date(d1)
        end = ee.Date(d1).advance(interval, "month")
        date = ee.Date(d1)

        # Load tree species image for the date range
        species = (
            ee.ImageCollection(
                "projects/sat-io/open-datasets/"
                "CA_FOREST/SPECIES-1984-2022"
            )
            .filterDate(start, end)
            .first()
            .clip(aoi)
        )

        # Calculate per-pixel proportions (0/1 masks)
        def _calculate_species_proportions(image):
            proportions = [
                image.eq(value).rename(
                    "Proportion_" + str(value)
                )
                for value in SPECIES_VALUES
            ]
            return ee.Image(proportions)

        species_proportions = (
            _calculate_species_proportions(species)
            .unmask(0)
            .clip(aoi)
        )

        # Rename bands with readable snake_case names
        renamed = ee.Image(
            [
                species_proportions.select(index).rename(
                    SPECIES_NAMES[index]
                )
                for index in range(len(SPECIES_VALUES))
            ]
        )

        return renamed.set("year", date.get("year"))

    # Generate annual tree species proportions
    return ee.ImageCollection(
        [_species_ts(d) for d in dates]
    ).map(lambda img: img.clip(aoi))
