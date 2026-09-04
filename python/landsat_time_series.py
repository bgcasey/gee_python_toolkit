# ---
# title:   Landsat Time Series Analysis
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - Landsat 5/7/8/9 Surface Reflectance collections
#     (LANDSAT/*/C02/T1_L2)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - One annual multiband spectral-index GeoTIFF per date,
#     either aggregated to the ABMI 1 km reference grid or at
#     NATIVE_SCALE_M, selected by EXPORT_TARGET
#   - Focal (neighbourhood) GeoTIFFs at 0/150/250 m when
#     RUN_FOCAL is True, written to the same target as the
#     time series (EXPORT_TARGET)
#   - Per-band min/max summary CSV (image_stats)
# notes:
#   Builds an annual date list, computes the selected spectral
#   indices via the shared ls_fn helper, drops the QA_PIXEL
#   band, and casts to Float32.
#
#   BASE_SCALE_M is the base each year is pinned to before
#   aggregating, keeping the 30 m -> 1 km jump under Earth
#   Engine's per-tile reprojection limit.
#
#   NDRS, enabled by RUN_NDRS, is normalized against the
#   forest pixels of each annual composite (CA_FOREST_LC_VLCE2
#   land cover, capped at 2019), so it costs one aoi-wide
#   min/max per band per year. LS_NDRS_FOREST_TYPES gives one
#   band per forest class group (NDRS_coni, NDRS_deci,
#   NDRS_mixed).
#
#   Deviation: the shared ls_fn helper expects a client-side
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
from utils.compute_report import ComputeReport
from utils.gee_helpers import (
    calculate_image_collection_stats,
    create_date_list,
    export_image_collection,
    export_stats_to_csv,
    focal_stats,
)
from utils.gee_utils import (
    define_study_area,
    export_collection_to_reference_grid,
    initialize_ee,
)
from utils.landsat_time_series import AVAILABLE_INDICES, ls_fn

# 1. Setup ----

# 1.1 User parameters ----
FILE_PREFIX = "landsat_multiband"

LS_START_DATE = "2024-06-01"  # first time-series date
LS_END_DATE = "2025-06-01"  # last time-series date
LS_DATE_INTERVAL = 1  # step between series start dates
LS_DATE_INTERVAL_TYPE = "years"  # units for the step
LS_WINDOW = 121  # compositing window length
LS_WINDOW_TYPE = "days"  # units for the window
LS_STATISTIC = "mean"  # 'mean', 'median', 'max', etc.
LS_INDICES = [
    "BSI", "DRS", "DSWI", "EVI", "GNDVI",
    "LAI", "NBR", "NDMI", "NDSI", "NDVI",
    "NDWI", "SAVI", "SI",
]

# NDRS bands (section 3), each costing an aoi-wide DRS
# min/max per year on top of the indices above.
RUN_NDRS = False  # False: skip the NDRS bands
# One NDRS band per entry, normalized within those forest
# classes: 210 coniferous (NDRS_coni), 220 broadleaf
# (NDRS_deci), 230 mixedwood; any other combination gives
# NDRS_mixed.
LS_NDRS_FOREST_TYPES = [[210], [220], [210, 220, 230]]

# Base each year is pinned to before aggregating; see notes.
BASE_SCALE_M = 30
NATIVE_SCALE_M = 30  # resolution of the "native" export

# Separate focal product (section 5), written to the same
# target as the time series.
RUN_FOCAL = False  # False: skip the focal exports
FOCAL_KERNELS = [150, 250]  # focal radii (m), circle

FOCAL_REACH_M = max(FOCAL_KERNELS) if RUN_FOCAL else 0

# Reach of the nearest-neighbour gap fill applied after
# aggregation, in 1 km cells; 0 disables it.
FILL_GAPS_PX = 20

# "native" skips the aggregation and writes NATIVE_SCALE_M
# pixels in the grid CRS, ungridded - useful for inspecting
# the input to the aggregation.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"

USE_TEST_AOI = True  # True: small test AOI; False: Alberta
PRINT_STATS = False  # value preview (slow for large AOIs)
COMPUTE_REPORT = False  # write EECU usage report (txt)
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
unknown_indices = [
    index
    for index in LS_INDICES
    if index not in AVAILABLE_INDICES
]
if unknown_indices:
    raise ValueError(
        f"Unknown LS_INDICES: {unknown_indices}; "
        f"available: {list(AVAILABLE_INDICES)}"
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
# ls_fn iterates a client-side list, so the dates are
# materialized as YYYY-MM-dd strings.
date_list = create_date_list(
    ee.Date(LS_START_DATE),
    ee.Date(LS_END_DATE),
    LS_DATE_INTERVAL,
    LS_DATE_INTERVAL_TYPE,
)
start_dates = date_list.map(
    lambda d: ee.Date(d).format("YYYY-MM-dd")
).getInfo()

ls = ls_fn(
    start_dates,
    LS_WINDOW,
    LS_WINDOW_TYPE,
    aoi_compute,
    LS_INDICES + (["NDRS"] if RUN_NDRS else []),
    LS_STATISTIC,
    ndrs_forest_types=LS_NDRS_FOREST_TYPES,
)


def drop_qa_and_cast(image):
    """Drop the QA_PIXEL band and cast bands to Float32."""
    keep = image.bandNames().filter(
        ee.Filter.neq("item", "QA_PIXEL")
    )
    return image.select(keep).toFloat()


ls = ls.map(drop_qa_and_cast)

# 3.1 Check layer values (optional) ----
# Earth Engine is lazy, so the profiler needs an evaluated
# computation to measure EECU usage; this also runs when
# COMPUTE_REPORT is on. The full per-band table also goes out
# as a CSV.
if PRINT_STATS or COMPUTE_REPORT:
    reducer = ee.Reducer.min().combine(
        ee.Reducer.max(), "", True
    )
    collection_stats = calculate_image_collection_stats(
        ls, aoi, COARSE_SCALE, 1e13, reducer
    )
    with report.section("Landsat band min/max stats"):
        stats = collection_stats.first().toDictionary().getInfo()
    print("Landsat first-image stats:", stats)
    export_stats_to_csv(
    collection_stats, f"image_stats{TEST_SUFFIX}"
)

# 4. Export the time series ----
# One export task per image. Monitor progress at
# https://code.earthengine.google.com/tasks
target_suffix = (
    "abmi1km" if EXPORT_TARGET == "reference_grid" else "native"
) + TEST_SUFFIX


def landsat_file_name(img):
    """File name for one multiband export."""
    year = img.get("year").getInfo() or "unknown"
    return f"{FILE_PREFIX}_{year}_{target_suffix}"


def export_to_target(collection, file_name_fn):
    """Export a collection to the EXPORT_TARGET grid."""
    if EXPORT_TARGET == "reference_grid":
        return export_collection_to_reference_grid(
            collection,
            aoi,
            file_name_fn,
            folder=DRIVE_FOLDER,
            reducer=ee.Reducer.mean(),
            agg_base_m=BASE_SCALE_M,
            agg_max_pixels=AGG_MAX_PIXELS,
            fill_gaps_px=FILL_GAPS_PX,
        )
    return export_image_collection(
        collection,
        aoi,
        DRIVE_FOLDER,
        NATIVE_SCALE_M,
        GRID_CRS,
        file_name_fn,
    )


export_tasks += export_to_target(ls, landsat_file_name)

# 5. Focal analysis ----
# An optional separate product, controlled by RUN_FOCAL: focal
# (neighbourhood) statistics at 0/150/250 m, exported to the
# same target as section 4. The 0 m case appends a "_0" band
# suffix but applies no smoothing.

if RUN_FOCAL:

    def make_focal_file_name(kernel_size):
        """Build the file-name function for one focal radius."""

        def focal_file_name(img):
            year = img.get("year").getInfo() or "unknown"
            return (
                f"{FILE_PREFIX}_{kernel_size}_{year}"
                f"_{target_suffix}"
            )

        return focal_file_name

    def rename_zero_focal(img):
        """Append a "_0" suffix to every band name."""
        new_names = img.bandNames().map(
            lambda name: ee.String(name).cat("_0")
        )
        return img.rename(new_names)

    export_tasks += export_to_target(
        ls.map(rename_zero_focal), make_focal_file_name(0)
    )

    for kernel_size in FOCAL_KERNELS:
        ls_focal = ls.map(
            lambda img, k=kernel_size: focal_stats(
                img, k, "circle", ["year"]
            )
        )
        export_tasks += export_to_target(
            ls_focal, make_focal_file_name(kernel_size)
        )

# 6. Compute usage report ----
# Records each task's total EECU-seconds and writes the txt
# report to gee_compute_reports/.

if WAIT_FOR_EXPORTS:
    for task in export_tasks:
        report.log_task(task)
report.write()

# End of script ----
