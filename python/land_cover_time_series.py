# ---
# title:   Annual Forest Land Cover Time Series
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - High-resolution Annual Forest Land Cover Maps for
#     Canada (projects/sat-io/open-datasets/
#     CA_FOREST_LC_VLCE2)
#   - FAO GAUL province boundaries (Alberta)
# outputs:
#   - Annual land cover GeoTIFFs (forest_lc_class) exported
#     to Google Drive at native (30 m) resolution
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

from _gee_config import DRIVE_FOLDER
from utils.annual_forest_land_cover import lc_fn
from utils.compute_report import ComputeReport
from utils.gee_helpers import export_image_collection
from utils.gee_utils import initialize_ee

# 1. Setup ----

# 1.1 User parameters ----
LC_START_DATE = "2000-01-01"  # first land cover year
LC_END_DATE = "2019-12-31"  # last available year (2019)
EXPORT_SCALE = 30  # native land cover resolution (m)
EXPORT_CRS = "EPSG:4326"  # native export CRS
PRINT_STATS = True  # summary check (slow for large AOIs)
USE_TEST_AOI = True  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt)

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


export_image_collection(
    lc,
    aoi,
    DRIVE_FOLDER,
    EXPORT_SCALE,
    EXPORT_CRS,
    land_cover_file_name,
)

# 5. Compute usage report ----
# Writes the profiled sections to gee_compute_reports/.
# Collection exports start many batch tasks, so per-task
# EECU totals are not logged here; monitor progress at
# https://code.earthengine.google.com/tasks
report.write()

# End of script ----
