# ---
# title:   Annual Forest Land Cover Time Series
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - High-resolution Annual Forest Land Cover Maps for
#     Canada (projects/sat-io/open-datasets/
#     CA_FOREST_LC_VLCE2)
#   - AB2020 provincial boundary (EE asset)
# outputs:
#   - Annual land cover GeoTIFFs (forest_lc_class) exported to
#     Google Drive, at native (30 m) or aggregated to the ABMI
#     1 km reference grid by modal class, per EXPORT_TARGET.
# notes:
#   Runnable land cover time-series script for the Earth
#   Engine Python API. The original
#   land_cover_time_series.js source was empty, so this
#   script mirrors the structure of the other time-series
#   scripts and uses the shared lc_fn helper to retrieve
#   annual forest land cover for a date range, then exports
#   one GeoTIFF per year.
#
#   Deviation: land cover is categorical, so no focal-mean
#   smoothing is applied (unlike the continuous Landsat and
#   MODIS time series). Only native-resolution annual
#   images are exported.
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

import os
import sys

import ee

# Make utils importable regardless of the working
# directory VS Code runs the script from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gee_config import DRIVE_FOLDER, PROVINCIAL_BOUNDARY_ASSET
from utils.annual_forest_land_cover import lc_fn
from utils.compute_report import ComputeReport
from utils.gee_helpers import export_image_collection
from utils.gee_utils import (
    export_collection_to_reference_grid,
    initialize_ee,
)

# 1. Setup ----

# 1.1 User parameters ----
LC_START_DATE = "2000-01-01"  # first land cover year
LC_END_DATE = "2019-12-31"  # last available year (2019)
EXPORT_SCALE = 30  # native land cover resolution (m)
EXPORT_CRS = "EPSG:3400"  # native export CRS (AB 10-TM)
# Export target for the annual surfaces. "native" writes each
# year at 30 m; "reference_grid" aggregates each year onto the
# ABMI 1 km grid by MODAL class (land cover is categorical, so
# mean is meaningless). AGG_BASE_M is the coarse base each year
# is pinned to before aggregating so the 30 m -> 1 km jump stays
# under the reprojection limit; because mode needs fine pixels,
# the base subsamples classes at AGG_BASE_M spacing, so the
# 1 km class is the dominant class of that subsample.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"
AGG_BASE_M = 90  # aggregation base (m) for the grid path
PRINT_STATS = True  # summary check (slow for large AOIs)
USE_TEST_AOI = True  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt)
# Block until every export task finishes so its batch
# EECU-seconds land in the compute report. Costs the full
# export runtime (hours for a province-wide run), so keep it
# False for production runs and turn it on when profiling a
# test AOI.
WAIT_FOR_EXPORTS = False

# Export tasks started below, for the optional per-task EECU
# logging in the compute-report section at the end.
export_tasks = []

# 1.2 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.3 Set up compute usage report ----
# Profiles EECU usage per section. Best used with
# USE_TEST_AOI = True to find choke points cheaply.
report = ComputeReport(
    "land_cover_time_series",
    enabled=COMPUTE_REPORT,
)

# 2. Define study area ----
# Uses a small test polygon when USE_TEST_AOI is True;
# otherwise uses the AB2020 provincial boundary asset.

if USE_TEST_AOI:
    # Small aoi for testing purposes
    aoi = ee.Geometry.Polygon([
        [-113.5, 55.5],  # Top-left corner
        [-113.5, 55.0],  # Bottom-left corner
        [-112.8, 55.0],  # Bottom-right corner
        [-112.8, 55.5],  # Top-right corner
    ])
else:
    aoi = ee.FeatureCollection(
        PROVINCIAL_BOUNDARY_ASSET
    ).geometry()

# 3. Land cover time-series processing ----
# Retrieves annual forest land cover images for the date
# range, each with a single 'forest_lc_class' band clipped
# to the AOI.
lc = lc_fn(LC_START_DATE, LC_END_DATE, aoi)

# 3.1 Check calculated bands (optional) ----
# Earth Engine is lazy, so the profiler needs an evaluated
# computation to measure per-algorithm EECU usage.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("Land cover class frequencies"):
        histogram = (
            lc.first()
            .reduceRegion(
                reducer=ee.Reducer.frequencyHistogram(),
                geometry=aoi,
                scale=EXPORT_SCALE,
                maxPixels=1e13,
                bestEffort=True,
            )
            .getInfo()
        )
    print("Land cover class frequencies:", histogram)

# 4. Export time series to Google Drive ----
# Exports each annual land cover image as a GeoTIFF, one
# export task per image.


def land_cover_file_name(img):
    """File name for the native-resolution export."""
    year = img.get("year").getInfo() or "unknown"
    return "forest_lc_class_" + str(year)


if EXPORT_TARGET == "reference_grid":
    export_tasks += export_collection_to_reference_grid(
        lc,
        aoi,
        lambda img: land_cover_file_name(img) + "_abmi1km",
        folder=DRIVE_FOLDER,
        reducer=ee.Reducer.mode(),
        agg_base_m=AGG_BASE_M,
        round_values=True,
    )
elif EXPORT_TARGET == "native":
    export_tasks += export_image_collection(
        lc,
        aoi,
        DRIVE_FOLDER,
        EXPORT_SCALE,
        EXPORT_CRS,
        lambda img: land_cover_file_name(img) + "_native",
    )
else:
    raise ValueError(
        "Unknown EXPORT_TARGET: "
        f"{EXPORT_TARGET!r} (use 'native' or 'reference_grid')"
    )

# 5. Compute usage report ----
# Writes the profiled sections to gee_compute_reports/.
# Collection exports start many batch tasks, so per-task
# EECU totals are not logged here; monitor progress at
# https://code.earthengine.google.com/tasks
if WAIT_FOR_EXPORTS:
    for task in export_tasks:
        report.log_task(task)
report.write()

# End of script ----
