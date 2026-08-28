# ---
# title:   Annual Forest Land Cover Time Series
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - High-resolution Annual Forest Land Cover Maps for
#     Canada (projects/sat-io/open-datasets/
#     CA_FOREST_LC_VLCE2)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - One annual land cover GeoTIFF (forest_lc_class) per
#     year, either aggregated to the ABMI 1 km reference grid
#     by modal class or at NATIVE_SCALE_M, selected by
#     EXPORT_TARGET (exported to Google Drive)
# notes:
#   Land cover is categorical, so the 1 km product aggregates
#   by modal class rather than by mean, and no focal-mean
#   smoothing is applied (unlike the continuous Landsat and
#   MODIS time series).
#
#   BASE_SCALE_M is the base each year is pinned to before
#   aggregating, keeping the 30 m -> 1 km jump under Earth
#   Engine's per-tile reprojection limit. Because mode needs
#   fine pixels, the base subsamples classes at BASE_SCALE_M
#   spacing, so the 1 km class is the dominant class of that
#   subsample.
#
#   Data citation: Hermosilla, T., Wulder, M.A., White,
#   J.C., Coops, N.C., 2022. Land cover classification in
#   an era of big and open data. Remote Sensing of
#   Environment. No. 112780. doi:10.1016/j.rse.2022.112780
#
#   Setup (once):
#     pip install earthengine-api
#     earthengine authenticate
#   Then set EE_PROJECT in _gee_config.py and run.
# ---

import math
import os
import sys

import ee

# Make utils importable regardless of the working
# directory VS Code runs the script from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gee_config import COARSE_SCALE, DRIVE_FOLDER, GRID_CRS
from utils.annual_forest_land_cover import lc_fn
from utils.compute_report import ComputeReport
from utils.gee_helpers import export_image_collection
from utils.gee_utils import (
    define_study_area,
    export_collection_to_reference_grid,
    initialize_ee,
)

# 1. Setup ----

# 1.1 User parameters ----
FILE_PREFIX = "forest_lc_class"

LC_START_DATE = "2000-01-01"  # first land cover year
LC_END_DATE = "2019-12-31"  # last available year (2019)

# Base each year is pinned to before aggregating; see notes.
BASE_SCALE_M = 90
NATIVE_SCALE_M = 30  # resolution of the "native" export

# Each pixel is read as a stored class, so no neighbourhood.
FOCAL_REACH_M = 0

# "native" skips the aggregation and writes NATIVE_SCALE_M
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
    else NATIVE_SCALE_M
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

# 3. Build the collection ----
# Annual land cover images for the date range, each with a
# single 'forest_lc_class' band. Read over aoi_compute, never
# aoi; clipped back in section 4.
lc = lc_fn(LC_START_DATE, LC_END_DATE, aoi_compute)

# 3.1 Check layer values (optional) ----
# Earth Engine is lazy, so the profiler needs an evaluated
# computation to measure EECU usage; this also runs when
# COMPUTE_REPORT is on.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("Land cover class frequencies"):
        stats = (
            lc.first()
            .reduceRegion(
                reducer=ee.Reducer.frequencyHistogram(),
                geometry=aoi,
                scale=NATIVE_SCALE_M,
                maxPixels=1e13,
                bestEffort=True,
            )
            .getInfo()
        )
    print("Land cover class frequencies:", stats)

# 4. Export the time series ----
# One export task per year. Monitor progress at
# https://code.earthengine.google.com/tasks
target_suffix = (
    "abmi1km" if EXPORT_TARGET == "reference_grid" else "native"
)


def land_cover_file_name(img):
    """File name for one annual land cover export."""
    year = img.get("year").getInfo() or "unknown"
    return f"{FILE_PREFIX}_{year}_{target_suffix}"


if EXPORT_TARGET == "reference_grid":
    export_tasks += export_collection_to_reference_grid(
        lc,
        aoi,
        land_cover_file_name,
        folder=DRIVE_FOLDER,
        reducer=ee.Reducer.mode(),
        agg_base_m=BASE_SCALE_M,
        agg_max_pixels=AGG_MAX_PIXELS,
        round_values=True,
    )
else:
    export_tasks += export_image_collection(
        lc,
        aoi,
        DRIVE_FOLDER,
        NATIVE_SCALE_M,
        GRID_CRS,
        land_cover_file_name,
    )

# 5. Compute usage report ----
# Records each task's total EECU-seconds and writes the txt
# report to gee_compute_reports/.

if WAIT_FOR_EXPORTS:
    for task in export_tasks:
        report.log_task(task)
report.write()

# End of script ----
