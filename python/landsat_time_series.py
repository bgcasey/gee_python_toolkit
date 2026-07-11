# ---
# title:   Landsat Time Series Analysis
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - Landsat 5/7/8/9 Surface Reflectance collections
#     (LANDSAT/*/C02/T1_L2)
#   - FAO GAUL province boundaries (Alberta)
# outputs:
#   - Annual multiband spectral-index GeoTIFFs exported to
#     Google Drive at native (30 m) resolution and at
#     focal scales (0/150/250 m) in EPSG:3978
#   - Per-band min/max summary CSV (image_stats)
# notes:
#   Python port of landsat_time_series.js for the Earth
#   Engine Python API. Builds an annual date list, computes
#   user-selected spectral indices via the shared ls_fn
#   helper, drops the QA_PIXEL band, casts to Float32, and
#   exports multiband images plus focal derivatives.
#
#   The original Map.addLayer/Map.centerObject calls, vis
#   parameters, and debug print() blocks are omitted.
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

from _gee_config import DRIVE_FOLDER
from utils.compute_report import ComputeReport
from utils.gee_helpers import (
    calculate_image_collection_stats,
    create_date_list,
    export_image_collection,
    export_stats_to_csv,
    focal_stats,
)
from utils.gee_utils import initialize_ee
from utils.landsat_time_series import ls_fn

# 1. Setup ----

# 1.1 User parameters ----
LS_START_DATE = "2000-06-01"  # first time-series date
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
EXPORT_SCALE = 30  # native Landsat resolution (m)
EXPORT_CRS = "EPSG:4326"  # native export CRS
FOCAL_SCALE = 990  # focal export scale (m)
FOCAL_CRS = "EPSG:3978"  # focal export CRS
FOCAL_KERNELS = [150, 250]  # focal radii (m), circle
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
    "landsat_time_series",
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
# create_date_list returns an ee.List; ls_fn iterates a
# client-side list, so the dates are materialized as
# YYYY-MM-dd strings.
date_list = create_date_list(
    ee.Date(LS_START_DATE),
    ee.Date(LS_END_DATE),
    LS_DATE_INTERVAL,
    LS_DATE_INTERVAL_TYPE,
)
start_dates = (
    date_list.map(
        lambda d: ee.Date(d).format("YYYY-MM-dd")
    ).getInfo()
)

# 4. Landsat time-series processing ----
# Computes the selected spectral indices for each interval,
# drops the QA_PIXEL band, and casts every band to Float32.
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

# 5. Check calculated bands (optional) ----
# Computes per-band min/max for the collection, prints the
# first-image summary, and exports the full table as a CSV.
if PRINT_STATS or COMPUTE_REPORT:
    reducer = ee.Reducer.min().combine(
        ee.Reducer.max(), "", True
    )
    collection_stats = calculate_image_collection_stats(
        ls, aoi, 1000, 1e13, reducer
    )
    with report.section("Landsat band min/max stats"):
        summary = (
            collection_stats.first()
            .toDictionary()
            .getInfo()
        )
    print("Landsat first-image stats:", summary)
    export_stats_to_csv(collection_stats, "image_stats")

# 6. Export time series to Google Drive ----
# Exports each image in the collection as a multiband
# GeoTIFF, one export task per image.


def landsat_file_name(img):
    """File name for the native-resolution export."""
    year = img.get("year").getInfo() or "unknown"
    return "landsat_multiband_" + str(year)


export_image_collection(
    ls,
    aoi,
    DRIVE_FOLDER,
    EXPORT_SCALE,
    EXPORT_CRS,
    landsat_file_name,
)

# 7. Focal analysis ----
# Exports focal (neighbourhood) statistics at 0/150/250 m
# in EPSG:3978. The 0 m case renames bands with a "_0"
# suffix but applies no smoothing.

# 7.1 Zero-metre focal (no smoothing) ----


def rename_zero_focal(img):
    """Append a "_0" suffix to every band name."""
    new_names = img.bandNames().map(
        lambda name: ee.String(name).cat("_0")
    )
    return img.rename(new_names)


ls_0 = ls.map(rename_zero_focal)


def landsat_file_name_0(img):
    """File name for the 0 m focal export."""
    year = img.get("year").getInfo() or "unknown"
    return "landsat_multiband_0_" + str(year)


export_image_collection(
    ls_0,
    aoi,
    DRIVE_FOLDER,
    FOCAL_SCALE,
    FOCAL_CRS,
    landsat_file_name_0,
)

# 7.2 Circular focal means (150 m, 250 m) ----
for kernel_size in FOCAL_KERNELS:
    ls_focal = ls.map(
        lambda img, k=kernel_size: focal_stats(
            img, k, "circle", ["year"]
        )
    )

    def make_focal_file_name(k):
        def focal_file_name(img):
            year = img.get("year").getInfo() or "unknown"
            return "landsat_multiband_" + str(k) + "_" + str(
                year
            )

        return focal_file_name

    export_image_collection(
        ls_focal,
        aoi,
        DRIVE_FOLDER,
        FOCAL_SCALE,
        FOCAL_CRS,
        make_focal_file_name(kernel_size),
    )

# 8. Compute usage report ----
# Writes the profiled sections to gee_compute_reports/.
# Collection exports start many batch tasks, so per-task
# EECU totals are not logged here; monitor progress at
# https://code.earthengine.google.com/tasks
report.write()

# End of script ----
