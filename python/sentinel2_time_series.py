# ---
# title:   Sentinel-2 Time Series Analysis
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - Sentinel-2 Surface Reflectance collection
#     (COPERNICUS/S2_SR_HARMONIZED)
#   - CA_FOREST_LC_VLCE2 land cover (for NDRS masks)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - One annual multiband spectral-index GeoTIFF per date,
#     either aggregated to the ABMI 1 km reference grid or at
#     NATIVE_SCALE_M, selected by EXPORT_TARGET (exported to
#     Google Drive)
# notes:
#   Builds an annual date list, computes the selected spectral
#   indices via the shared s2_fn helper, adds NDRS bands for
#   coniferous (210), broadleaf (220), and mixedwood (all)
#   forest, and casts to Float32.
#
#   BASE_SCALE_M is the base each year is pinned to before
#   aggregating, keeping the 10 m -> 1 km jump under Earth
#   Engine's per-tile reprojection limit.
#
#   Deviation: the shared s2_fn helper expects a client-side
#   list of date strings, so the ee.List produced by
#   create_date_list is materialized with getInfo() before
#   being passed in.
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
from utils import sentinel_indices_and_masks as indices
from utils.compute_report import ComputeReport
from utils.gee_helpers import (
    create_date_list,
    export_image_collection,
)
from utils.gee_utils import (
    define_study_area,
    export_collection_to_reference_grid,
    initialize_ee,
)
from utils.sentinel_time_series import s2_fn

# 1. Setup ----

# 1.1 User parameters ----
FILE_PREFIX = "sentinel2_multiband"

S2_START_DATE = "2023-06-01"  # first time-series date
S2_END_DATE = "2024-06-01"  # last time-series date
S2_DATE_INTERVAL = 1  # step between series start dates
S2_DATE_INTERVAL_TYPE = "years"  # units for the step
S2_WINDOW = 121  # compositing window length
S2_WINDOW_TYPE = "days"  # units for the window
S2_INDICES = [
    "CRE", "DRS", "DSWI", "EVI", "GNDVI", "LAI", "NBR",
    "NDRE1", "NDRE2", "NDRE3", "NDVI", "NDWI", "RDI",
]

# Base each year is pinned to before aggregating; see notes.
BASE_SCALE_M = 50
NATIVE_SCALE_M = 10  # resolution of the "native" export

# Indices are per-pixel, so no neighbourhood is read.
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
# Reads over aoi_compute, never aoi; clipped back in section 4.
# s2_fn iterates a client-side list, so the dates are
# materialized as YYYY-MM-dd strings.
date_list = create_date_list(
    ee.Date(S2_START_DATE),
    ee.Date(S2_END_DATE),
    S2_DATE_INTERVAL,
    S2_DATE_INTERVAL_TYPE,
)
start_dates = date_list.map(
    lambda d: ee.Date(d).format("YYYY-MM-dd")
).getInfo()

s2 = s2_fn(
    start_dates,
    S2_WINDOW,
    S2_WINDOW_TYPE,
    aoi_compute,
    S2_INDICES,
)
s2 = (
    s2.map(lambda img: indices.add_ndrs(img, [210]))
    .map(lambda img: indices.add_ndrs(img, [220]))
    .map(lambda img: indices.add_ndrs(img))
    .map(lambda img: img.toFloat())
)

# 3.1 Check layer values (optional) ----
# Earth Engine is lazy, so the profiler needs an evaluated
# computation to measure EECU usage; this also runs when
# COMPUTE_REPORT is on.
if PRINT_STATS or COMPUTE_REPORT:
    reducer = (
        ee.Reducer.min()
        .combine(ee.Reducer.max(), "", True)
        .combine(ee.Reducer.stdDev(), "", True)
    )
    with report.section("Sentinel-2 first-image stats"):
        stats = (
            s2.first()
            .reduceRegion(
                reducer=reducer,
                geometry=aoi,
                scale=COARSE_SCALE,
                bestEffort=True,
                maxPixels=1e13,
            )
            .getInfo()
        )
    print("Sentinel-2 first-image stats:", stats)

# 4. Export the time series ----
# One export task per image. Monitor progress at
# https://code.earthengine.google.com/tasks
target_suffix = (
    "abmi1km" if EXPORT_TARGET == "reference_grid" else "native"
)


def sentinel_file_name(img):
    """File name for one multiband export.

    Parameters
    ----------
    img : ee.Image
        Annual composite carrying a 'year' property.

    Returns
    -------
    str
        File name, e.g. 'sentinel2_multiband_2023_abmi1km'.

    Examples
    --------
    >>> sentinel_file_name(s2.first())
    'sentinel2_multiband_2023_abmi1km'
    """
    year = img.get("year").getInfo() or "unknown"
    return f"{FILE_PREFIX}_{year}_{target_suffix}"


if EXPORT_TARGET == "reference_grid":
    export_tasks += export_collection_to_reference_grid(
        s2,
        aoi,
        sentinel_file_name,
        folder=DRIVE_FOLDER,
        reducer=ee.Reducer.mean(),
        agg_base_m=BASE_SCALE_M,
        agg_max_pixels=AGG_MAX_PIXELS,
    )
else:
    export_tasks += export_image_collection(
        s2,
        aoi,
        DRIVE_FOLDER,
        NATIVE_SCALE_M,
        GRID_CRS,
        sentinel_file_name,
    )

# 5. Compute usage report ----
# Records each task's total EECU-seconds and writes the txt
# report to gee_compute_reports/.

if WAIT_FOR_EXPORTS:
    for task in export_tasks:
        report.log_task(task)
report.write()

# End of script ----
