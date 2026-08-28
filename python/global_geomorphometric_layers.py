# ---
# title:   Global Geomorphometric Layers (Geomorpho90m)
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - Geomorpho90m ImageCollections
#     (projects/sat-io/open-datasets/Geomorpho90m)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - Multiband Geomorpho90m GeoTIFF for Alberta, either
#     aligned to the ABMI 1 km reference grid or at
#     BASE_SCALE_M resolution, selected by EXPORT_TARGET
#     (exported to Google Drive). The 1 km product drops the
#     raw 'aspect' band (see notes).
# notes:
#   Loads the Geomorpho90m geomorphometric variables, mosaics
#   and clips each to the AOI, and stacks them into a single
#   multiband image. Geomorpho90m is a stored (pyramided)
#   dataset, so ~90 m -> 1 km is well under Earth Engine's
#   per-tile reprojection limit.
#
#   The raw 'aspect' band is dropped from the 1 km product:
#   averaging a circular angle (0-360 deg) is meaningless.
#   'aspect-cosine' and 'aspect-sine' carry aspect correctly
#   and mean cleanly.
#
#   Citation:
#   Amatulli, G., McInerney, D., Sethi, T., Strobl, P.,
#   Domisch, S. (2020). Geomorpho90m, empirical evaluation
#   and accuracy assessment of global high-resolution
#   geomorphometric layers. Scientific Data 7(1), 1-18.
#
#   Setup (once):
#     pip install earthengine-api
#     earthengine authenticate
#   Then set EE_PROJECT in _gee_config.py to your
#   registered Earth Engine cloud project and run the
#   script.
# ---

import math
import os
import sys

import ee

# Make utils importable regardless of the working
# directory VS Code runs the script from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gee_config import COARSE_SCALE, DRIVE_FOLDER, GRID_CRS
from utils.compute_report import ComputeReport
from utils.gee_utils import (
    define_study_area,
    export_image_to_drive,
    export_to_reference_grid,
    initialize_ee,
)

# 1. Setup ----

# 1.1 User parameters ----
TASK_PREFIX = "Geomorpho90m_AB"
FILE_PREFIX = "global_geomorphometric_layers"

BASE_PATH = "projects/sat-io/open-datasets/Geomorpho90m/"
# Collections to combine, in the order they are stacked into
# the multiband export image.
COLLECTION_NAMES = [
    "aspect",           # Aspect
    "aspect-cosine",    # Aspect-Cosine
    "aspect-sine",      # Aspect-Sine
    "convergence",      # Convergence Index
    "cti",              # Compound Topographic Index (CTI)
    "dev-magnitude",    # Deviation Magnitude
    "dev-scale",        # Deviation Scale
    "eastness",         # Eastness
    "elev-stdev",       # Elevation Standard Deviation
    "northness",        # Northness
    "rough-magnitude",  # Multiscale Roughness Magnitude
    "rough-scale",      # Multiscale Roughness Scale
    "roughness",        # Roughness
    "slope",            # Slope
    "spi",              # Stream Power Index
    "tpi",              # Topographic Position Index (TPI)
    "tri",              # Terrain Ruggedness Index (TRI)
    "vrm",              # Vector Ruggedness Measure (VRM)
]
# Dropped from the 1 km product; see notes.
GRID_EXCLUDE_BANDS = ["aspect"]

# Geomorpho90m's nominal resolution. native_scale reports
# 92.766 for the underlying grid if you prefer it exact.
BASE_SCALE_M = 90

# Each layer is read as a stored value, so no neighbourhood.
FOCAL_REACH_M = 0

# "native" skips the aggregation and writes BASE_SCALE_M
# pixels in the grid CRS, ungridded - useful for inspecting
# the input to the aggregation.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"

USE_TEST_AOI = True  # True: small test AOI; False: Alberta
PRINT_STATS = True  # value preview (slow for large AOIs)
COMPUTE_REPORT = True  # write EECU usage report (txt)
# Costs the full export runtime (hours province-wide), so
# turn it on only when profiling the test AOI.
WAIT_FOR_EXPORTS = False

# 1.2 Validate parameters ----
if EXPORT_TARGET not in ("native", "reference_grid"):
    raise ValueError(
        "Unknown EXPORT_TARGET: "
        f"{EXPORT_TARGET!r} (use 'native' or 'reference_grid')"
    )

# 1.3 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.4 Derived scales ----
# AGG_BUFFER_M is 2x the output scale, so a 1 km cell touching
# the aoi at a corner still reads a full diagonal (1414 m)
# beyond it; AGG_MAX_PIXELS is reduceResolution's per-cell
# input budget, +2 covering a cell that straddles a base pixel
# on each side.
AGG_BUFFER_M = 2 * (
    COARSE_SCALE
    if EXPORT_TARGET == "reference_grid"
    else BASE_SCALE_M
)
COMPUTE_BUFFER_M = max(FOCAL_REACH_M, AGG_BUFFER_M)
AGG_MAX_PIXELS = (
    math.ceil(COARSE_SCALE / BASE_SCALE_M) + 2
) ** 2

# 1.5 Set up run bookkeeping ----
# report profiles compute usage; export_tasks collects the
# export tasks it logs.
report = ComputeReport(FILE_PREFIX, enabled=COMPUTE_REPORT)
export_tasks = []

# 2. Define study area ----
# aoi is the export / crop boundary; aoi_compute adds the ring
# the source is read over.
aoi, aoi_compute = define_study_area(
    use_test_aoi=USE_TEST_AOI,
    buffer_m=COMPUTE_BUFFER_M,
)

# 3. Build the layer ----
# Reads over aoi_compute, never aoi; clipped back in section 4.


def load_band(collection_name, compute_aoi):
    """Mosaic one Geomorpho90m collection into a named band."""
    return (
        ee.ImageCollection(BASE_PATH + collection_name)
        .mosaic()
        .clip(compute_aoi)
        .rename(collection_name)
    )


layer = load_band(COLLECTION_NAMES[0], aoi_compute)
for name in COLLECTION_NAMES[1:]:
    layer = layer.addBands(load_band(name, aoi_compute))

# 3.1 Check layer values (optional) ----
# Earth Engine is lazy, so the profiler needs an evaluated
# computation to measure EECU usage; this also runs when
# COMPUTE_REPORT is on.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("Geomorpho90m min/max (reduceRegion)"):
        stats = layer.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=aoi,
            scale=BASE_SCALE_M,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print("Geomorpho90m min/max:", stats)

# 4. Aggregate to the grid and export ----
# setDefaultProjection pins the native base so
# reduceResolution knows the input resolution. Monitor
# progress at https://code.earthengine.google.com/tasks
target_suffix = (
    "abmi1km" if EXPORT_TARGET == "reference_grid" else "native"
)

if EXPORT_TARGET == "reference_grid":
    grid_bands = [
        n for n in COLLECTION_NAMES if n not in GRID_EXCLUDE_BANDS
    ]
    export_tasks.append(
        export_to_reference_grid(
            image=layer.select(grid_bands).setDefaultProjection(
                crs=GRID_CRS, scale=BASE_SCALE_M
            ),
            aoi=aoi,
            description=f"{TASK_PREFIX}_{target_suffix}",
            folder=DRIVE_FOLDER,
            file_name_prefix=f"{FILE_PREFIX}_{target_suffix}",
            agg_max_pixels=AGG_MAX_PIXELS,
            wait=False,
        )
    )
else:
    export_tasks.append(
        export_image_to_drive(
            image=layer.clip(aoi),
            description=f"{TASK_PREFIX}_{target_suffix}",
            region=aoi,
            folder=DRIVE_FOLDER,
            file_name_prefix=f"{FILE_PREFIX}_{target_suffix}",
            scale=BASE_SCALE_M,
            crs=GRID_CRS,
            max_pixels=1e13,
            wait=False,
        )
    )

# 5. Compute usage report ----
# Records each task's total EECU-seconds and writes the txt
# report to gee_compute_reports/.

if WAIT_FOR_EXPORTS:
    for task in export_tasks:
        report.log_task(task)
report.write()

# End of script ----
