# ---
# title:   Hydrologically Adjusted Elevations (HAND)
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - MERIT Hydro (MERIT/Hydro/v1_0_1)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - HAND GeoTIFF for Alberta, either aligned to the ABMI
#     1 km reference grid or at BASE_SCALE_M resolution,
#     selected by EXPORT_TARGET (exported to Google Drive)
# notes:
#   Height Above Nearest Drainage, the 'hnd' band of MERIT
#   Hydro. HAND is a stored (pyramided) band, so ~90 m -> 1 km
#   is well under Earth Engine's per-tile reprojection limit.
#
#   Data citation:
#   Yamazaki, D., et al. (2019). MERIT Hydro: A
#   high-resolution global hydrography map based on latest
#   topography datasets. Water Resources Research, 55,
#   5053-5073. doi:10.1029/2019WR024873
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
BAND_NAME = "HAND"
TASK_PREFIX = "HAND_AB"
FILE_PREFIX = "hydrologically_adjusted_elevations"

SOURCE_ASSET = "MERIT/Hydro/v1_0_1"
SOURCE_BAND = "hnd"

# MERIT Hydro's native grid spacing. native_scale(ee.Image(
# SOURCE_ASSET)) reports 92.766 if you prefer it exact.
BASE_SCALE_M = 92.77

# HAND is read as a stored value, so no neighbourhood.
FOCAL_REACH_M = 0

# Reach of the nearest-neighbour gap fill applied after
# aggregation, in 1 km cells; 0 disables it.
FILL_GAPS_PX = 20

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
if FILL_GAPS_PX < 0:
    raise ValueError(
        f"FILL_GAPS_PX must be >= 0, got {FILL_GAPS_PX!r}"
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
# setDefaultProjection pins the native base so reduceResolution
# knows the input resolution.
layer = (
    ee.Image(SOURCE_ASSET)
    .select(SOURCE_BAND)
    .clip(aoi_compute)
    .setDefaultProjection(crs=GRID_CRS, scale=BASE_SCALE_M)
    .rename(BAND_NAME)
)

# 3.1 Check layer values (optional) ----
# Earth Engine is lazy, so the profiler needs an evaluated
# computation to measure EECU usage; this also runs when
# COMPUTE_REPORT is on.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section(f"{BAND_NAME} min/max (reduceRegion)"):
        stats = layer.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=aoi,
            scale=BASE_SCALE_M,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print(f"{BAND_NAME} min/max:", stats)

# 4. Aggregate to the grid and export ----
# Monitor progress at
# https://code.earthengine.google.com/tasks
target_suffix = (
    "abmi1km" if EXPORT_TARGET == "reference_grid" else "native"
)

if EXPORT_TARGET == "reference_grid":
    export_tasks.append(
        export_to_reference_grid(
            image=layer,
            aoi=aoi,
            description=f"{TASK_PREFIX}_{target_suffix}",
            folder=DRIVE_FOLDER,
            file_name_prefix=f"{FILE_PREFIX}_{target_suffix}",
            agg_max_pixels=AGG_MAX_PIXELS,
            fill_gaps_px=FILL_GAPS_PX,
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
