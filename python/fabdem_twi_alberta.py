# ---
# title:   FABDEM Topographic Wetness Index for Alberta
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - FABDEM ImageCollection
#     (projects/sat-io/open-datasets/FABDEM)
#   - MERIT Hydro upslope area (MERIT/Hydro/v1_0_1)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - TWI GeoTIFF for Alberta, either aligned to the ABMI 1 km
#     reference grid or at BASE_SCALE_M resolution,
#     selected by EXPORT_TARGET (exported to Google Drive)
# notes:
#   Topographic Wetness Index, ln(a / tan(b)), where a is
#   upslope drainage area (m^2) and b is slope (radians).
#
#   FABDEM is a bare-earth DEM with no flow-accumulation
#   band, and Earth Engine has no native flow-accumulation
#   algorithm. This script therefore uses a hybrid: slope
#   from FABDEM (computed at BASE_SCALE_M) and upslope area
#   from MERIT Hydro 'upa' (~90 m).
#
#   Data citations:
#   Hawker, L., et al. (2022). A 30 m global map of
#   elevation with forests and buildings removed.
#   Environmental Research Letters, 17(2), 024016.
#   doi:10.1088/1748-9326/ac4d4f
#
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
    fabdem_elevation,
    initialize_ee,
)

# 1. Setup ----

# 1.1 User parameters ----
BAND_NAME = "twi"
TASK_PREFIX = "FABDEM_TWI_Alberta"
FILE_PREFIX = "fabdem_twi_alberta"

# Resolution slope is computed at before aggregating to 1 km.
BASE_SCALE_M = 30

# ee.Terrain.slope reads a 3x3 window, so one base pixel.
FOCAL_REACH_M = BASE_SCALE_M

# Reach of the nearest-neighbour gap fill applied after
# aggregation, in 1 km cells; 0 disables it.
FILL_GAPS_PX = 0

# "native" skips the aggregation and writes BASE_SCALE_M
# pixels in the grid CRS, ungridded - useful for inspecting
# the input to the aggregation.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"

USE_TEST_AOI = False  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt)
# Costs the full export runtime (hours province-wide), so
# turn it on only when profiling the test AOI.
WAIT_FOR_EXPORTS = False

# Derived: keeps test-AOI tasks and files
# distinguishable from full-extent runs.
TEST_SUFFIX = "_test" if USE_TEST_AOI else ""
    
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
elevation = fabdem_elevation(aoi_compute, base_m=BASE_SCALE_M)
slope_rad = ee.Terrain.slope(elevation).multiply(math.pi / 180)

# MERIT 'upa' is a stored (pyramided) dataset, so it needs no
# coarse base. km^2 -> m^2.
upslope_area = (
    ee.Image("MERIT/Hydro/v1_0_1")
    .select("upa")
    .clip(aoi_compute)
    .multiply(1e6)
)

# Floor tan(b) so flat areas (slope 0) are not masked by
# division by zero.
tan_b = slope_rad.tan().max(0.001)
layer = upslope_area.divide(tan_b).log().rename(BAND_NAME)

# 4. Aggregate to the grid and export ----
# Monitor progress at
# https://code.earthengine.google.com/tasks
target_suffix = (
    "abmi1km" if EXPORT_TARGET == "reference_grid" else "native"
) + TEST_SUFFIX

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
