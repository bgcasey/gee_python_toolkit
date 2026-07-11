# ---
# title:   Extract Geomorpho90m
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Load multiple geomorphometric variables from the
#   Geomorpho90m dataset, mosaic them, clip them to a given
#   area of interest, and combine them into a single
#   multiband image (with MERIT DEM elevation appended).
#
#   Data citation:
#   Amatulli, G., McInerney, D., Sethi, T., Strobl, P.,
#   Domisch, S. (2020). Geomorpho90m, empirical evaluation
#   and accuracy assessment of global high-resolution
#   geomorphometric layers. Scientific Data 7(1), 1-18.
# ---

import ee

# Base path for the Geomorpho90m collections
BASE_PATH = "projects/sat-io/open-datasets/Geomorpho90m/"

# Geomorphometric collection names to load and combine
COLLECTION_NAMES = [
    "aspect",  # Aspect
    "aspect-cosine",  # Aspect-Cosine
    "aspect-sine",  # Aspect-Sine
    "eastness",  # Eastness
    "northness",  # Northness
    "cti",  # Compound Topographic Index (CTI)
    "elev-stdev",  # Elevation Standard Deviation
    "vrm",  # Vector Ruggedness Measure (VRM)
    "roughness",  # Roughness
    "tri",  # Terrain Ruggedness Index (TRI)
    "tpi",  # Topographic Position Index (TPI)
    "dev-magnitude",  # Deviation Magnitude
    "dev-scale",  # Deviation Scale
    "rough-magnitude",  # Multiscale Roughness Magnitude
    "rough-scale",  # Multiscale Roughness Scale
]


def load_and_process(collection_name, aoi):
    """Load, mosaic, clip, and rename an image collection.

    Args:
        collection_name (str): Geomorpho90m collection name.
        aoi (ee.Geometry): Area of interest.

    Returns:
        ee.Image: Single-band image named after the
        collection, clipped to the AOI.
    """
    image = (
        ee.ImageCollection(BASE_PATH + collection_name)
        .mosaic()
        .clip(aoi)
        .rename(collection_name)
    )
    return image


def get_geomorpho90m(aoi):
    """Combine Geomorpho90m layers into one multiband image.

    Args:
        aoi (ee.Geometry): Area of interest.

    Returns:
        ee.Image: Multiband image of all Geomorpho90m layers
        with a MERIT DEM 'elev' band appended.
    """
    # MERIT DEM elevation clipped and renamed to 'elev'
    elevation = (
        ee.Image("MERIT/DEM/v1_0_3")
        .clip(aoi)
        .select(["dem"], ["elev"])
    )

    # Start from the first collection, then add the rest
    geomorpho90m = load_and_process(COLLECTION_NAMES[0], aoi)
    for name in COLLECTION_NAMES[1:]:
        geomorpho90m = geomorpho90m.addBands(
            load_and_process(name, aoi)
        )

    return geomorpho90m.addBands(elevation)
