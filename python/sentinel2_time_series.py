# ---
# title:   Sentinel-2 Time Series Analysis
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - Sentinel-2 Surface Reflectance collection
#     (COPERNICUS/S2_SR_HARMONIZED)
#   - CA_FOREST_LC_VLCE2 land cover (for NDRS masks)
#   - FAO GAUL province boundaries (Alberta)
# outputs:
#   - Annual multiband spectral-index GeoTIFFs exported to
#     Google Drive at native (10 m) resolution
# notes:
#   Python port of sentinel2_time_series.js for the Earth
#   Engine Python API. Builds an annual date list, computes
#   user-selected spectral indices via the shared s2_fn
#   helper, adds NDRS bands for coniferous (210), broadleaf
#   (220), and mixedwood (all) forest, casts to Float32,
#   and exports multiband images.
#
#   The original Map.addLayer/Map.centerObject calls, vis
#   parameters, and debug print() blocks are omitted.
#
#   Deviation: the shared s2_fn helper expects a client-
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

from _gee_config import DRIVE_FOLDER
from utils import sentinel_indices_and_masks as indices
from utils.compute_report import ComputeReport
from utils.gee_helpers import (
    create_date_list,
    export_image_collection,
)
from utils.gee_utils import initialize_ee
from utils.sentinel_time_series import s2_fn

# 1. Setup ----

# 1.1 User parameters ----
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
EXPORT_SCALE = 10  # native Sentinel-2 resolution (m)
EXPORT_CRS = "EPSG:4326"  # native export CRS
PRINT_STATS = True  # min/max check (slow for large AOIs)
USE_TEST_AOI = True  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt)

# 1.2 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.3 Set up compute usage report ----
# Profiles EECU usage per section. Best used with
# USE_TEST_AOI = True to find choke points cheaply.
report = ComputeReport(
    "sentinel2_time_series",
    enabled=COMPUTE_REPORT,
)

# 2. Define study area ----
# Uses a small test polygon when USE_TEST_AOI is True;
# otherwise filters the FAO GAUL provinces for Alberta.

if USE_TEST_AOI:
    # Small aoi for testing purposes
    aoi = ee.Geometry.Polygon([
        [-113.5, 55.5],  # Top-left corner
        [-113.5, 55.0],  # Bottom-left corner
        [-112.8, 55.0],  # Bottom-right corner
        [-112.8, 55.5],  # Top-right corner
    ])
else:
    aoi = (
        ee.FeatureCollection(
            "FAO/GAUL_SIMPLIFIED_500m/2015/level1"
        )
        .filter(ee.Filter.eq("ADM0_NAME", "Canada"))
        .filter(ee.Filter.eq("ADM1_NAME", "Alberta"))
        .geometry()
    )

# 3. Build the time-series date list ----
# create_date_list returns an ee.List; s2_fn iterates a
# client-side list, so the dates are materialized as
# YYYY-MM-dd strings.
date_list = create_date_list(
    ee.Date(S2_START_DATE),
    ee.Date(S2_END_DATE),
    S2_DATE_INTERVAL,
    S2_DATE_INTERVAL_TYPE,
)
start_dates = (
    date_list.map(
        lambda d: ee.Date(d).format("YYYY-MM-dd")
    ).getInfo()
)

# 4. Sentinel-2 time-series processing ----
# Computes the selected spectral indices for each interval,
# adds NDRS bands for coniferous (210), broadleaf (220),
# and mixedwood (all forest) pixels, and casts to Float32.
s2 = s2_fn(
    start_dates,
    S2_WINDOW,
    S2_WINDOW_TYPE,
    aoi,
    S2_INDICES,
)
s2 = (
    s2.map(lambda img: indices.add_ndrs(img, [210]))
    .map(lambda img: indices.add_ndrs(img, [220]))
    .map(lambda img: indices.add_ndrs(img))
    .map(lambda img: img.toFloat())
)

# 5. Check calculated bands (optional) ----
# Earth Engine is lazy, so the profiler needs an evaluated
# computation to measure per-algorithm EECU usage.
if PRINT_STATS or COMPUTE_REPORT:
    reducer = (
        ee.Reducer.min()
        .combine(ee.Reducer.max(), "", True)
        .combine(ee.Reducer.stdDev(), "", True)
    )
    with report.section("Sentinel-2 first-image stats"):
        stats_first = (
            s2.first()
            .reduceRegion(
                reducer=reducer,
                geometry=aoi,
                scale=1000,
                bestEffort=True,
                maxPixels=1e13,
            )
            .getInfo()
        )
    print("Sentinel-2 first-image stats:", stats_first)

# 6. Export time series to Google Drive ----
# Exports each image in the collection as a multiband
# GeoTIFF, one export task per image.


def sentinel_file_name(img):
    """File name for the native-resolution export."""
    year = img.get("year").getInfo() or "unknown"
    return "sentinel2_multiband_" + str(year)


export_image_collection(
    s2,
    aoi,
    DRIVE_FOLDER,
    EXPORT_SCALE,
    EXPORT_CRS,
    sentinel_file_name,
)

# 7. Compute usage report ----
# Writes the profiled sections to gee_compute_reports/.
# Collection exports start many batch tasks, so per-task
# EECU totals are not logged here; monitor progress at
# https://code.earthengine.google.com/tasks
report.write()

# End of script ----
