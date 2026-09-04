# ---
# title:   Summarize Landsat Time Series to Polygons
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - Landsat 5/7/8/9 Surface Reflectance collections
#     (LANDSAT/*/C02/T1_L2)
#   - Summary polygons (test FeatureCollection)
#   - AB2020 provincial boundary (EE asset)
# outputs:
#   - Per-polygon-per-date spectral-index summary CSV
#     (ls_poly_summary) exported to Google Drive
# notes:
#   Python port of landsat_time_series_to_poly.js for the
#   Earth Engine Python API. Builds an annual date list,
#   computes user-selected spectral indices via the shared
#   ls_fn helper, drops the QA_PIXEL band, casts to
#   Float32, and reduces each image over polygons with
#   image_collection_to_features, exporting the result as a
#   CSV.
#
#   The original Map.addLayer calls and debug print()
#   blocks are omitted.
#
#   Deviation: the shared ls_fn helper expects a client-
#   side list of date strings, so the ee.List produced by
#   create_date_list is materialized with getInfo() before
#   being passed in.
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

from utils.compute_report import ComputeReport
from utils.gee_helpers import create_date_list
from utils.gee_utils import define_study_area, initialize_ee
from utils.image_collection_to_features import (
    image_collection_to_features,
)
from utils.landsat_time_series import ls_fn

# 1. Setup ----

# 1.1 User parameters ----
LS_START_DATE = "2009-06-01"  # first time-series date
LS_END_DATE = "2024-06-01"  # last time-series date
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
SUMMARY_SCALE = 30  # reduction scale (m)
SUMMARY_CRS = "EPSG:4326"  # reduction CRS
SUMMARY_TILE_SCALE = 4  # tileScale for parallel reduction
SUMMARY_FILE_NAME = "ls_poly_summary"  # output CSV prefix
FILE_PREFIX = "landsat_time_series_to_poly"
PRINT_STATS = True  # value preview (slow for large AOIs)
USE_TEST_AOI = True  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt)

# Derived: keeps test-AOI tasks and files
# distinguishable from full-extent runs.
TEST_SUFFIX = "_test" if USE_TEST_AOI else ""

# 1.2 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.3 Set up run bookkeeping ----
# report profiles compute usage.
report = ComputeReport(FILE_PREFIX, enabled=COMPUTE_REPORT)

# 2. Define study area ----
# The summary polygons sit inside aoi, and the reduction is
# per polygon, so no compute ring is needed.
aoi, _ = define_study_area(use_test_aoi=USE_TEST_AOI)

# 2.1 Define summary polygons ----
# Five ~2 ha test polygons within the AOI, each tagged with
# an integer id.
poly1 = ee.Geometry.Polygon([
    [-113.48, 55.48],
    [-113.48, 55.47873],
    [-113.47779, 55.47873],
    [-113.47779, 55.48],
])
poly2 = ee.Geometry.Polygon([
    [-113.46, 55.47],
    [-113.46, 55.46873],
    [-113.45779, 55.46873],
    [-113.45779, 55.47],
])
poly3 = ee.Geometry.Polygon([
    [-113.44, 55.46],
    [-113.44, 55.45873],
    [-113.43779, 55.45873],
    [-113.43779, 55.46],
])
poly4 = ee.Geometry.Polygon([
    [-113.42, 55.45],
    [-113.42, 55.44873],
    [-113.41779, 55.44873],
    [-113.41779, 55.45],
])
poly5 = ee.Geometry.Polygon([
    [-113.40, 55.44],
    [-113.40, 55.43873],
    [-113.39779, 55.43873],
    [-113.39779, 55.44],
])

polys_fc = ee.FeatureCollection([
    ee.Feature(poly1, {"id": 1}),
    ee.Feature(poly2, {"id": 2}),
    ee.Feature(poly3, {"id": 3}),
    ee.Feature(poly4, {"id": 4}),
    ee.Feature(poly5, {"id": 5}),
])

# 3. Build the collection ----
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

# Selected spectral indices per interval, QA_PIXEL dropped,
# every band cast to Float32.
ls = ls_fn(
    start_dates,
    LS_WINDOW,
    LS_WINDOW_TYPE,
    aoi,
    LS_INDICES,
    LS_STATISTIC,
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
# COMPUTE_REPORT is on.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("Landsat band names"):
        band_names = ls.first().bandNames().getInfo()
    print("Landsat band names:", band_names)

# 4. Summarize the time series to polygons ----
# Applies the mean reducer to every band of each image over
# each polygon, producing a per-polygon-per-date summary
# table that is exported to Google Drive as a CSV.
ls_poly_summary = image_collection_to_features(
    ee.Reducer.mean(),
    polys_fc,
    aoi,
    ls,
    SUMMARY_CRS,
    SUMMARY_SCALE,
    SUMMARY_TILE_SCALE,
    f"{SUMMARY_FILE_NAME}{TEST_SUFFIX}",
)

# 5. Compute usage report ----
# Writes the profiled sections to gee_compute_reports/. The
# summary export runs as a batch table task; monitor progress
# at https://code.earthengine.google.com/tasks
report.write()

# End of script ----

# End of script ----
