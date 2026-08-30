# ---
# title:   FABDEM Elevation for US and Canada
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - FABDEM ImageCollection
#     (projects/sat-io/open-datasets/FABDEM)
#   - FAO GAUL country boundaries (FAO/GAUL/2015/level0)
# outputs:
#   - DEM GeoTIFF for US and Canada (exported to Google
#     Drive)
# notes:
#   This script downloads the FABDEM digital elevation
#   model for the United States and Canada. It mosaics the
#   collection, clips to the US + Canada boundaries, and
#   exports a GeoTIFF to Google Drive. Clipping directly to
#   the feature collection avoids geometry edge limits. Map
#   visualization layers (hillshade, ocean mask, elevation
#   palette) from the original GEE JavaScript are dropped.
#
#   Setup (once):
#     pip install earthengine-api
#     earthengine authenticate
#   Then set EE_PROJECT in _gee_config.py to your
#   registered Earth Engine cloud project and run the
#   script.
# ---

import os
import sys

import ee

# Make utils importable regardless of the working
# directory VS Code runs the script from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gee_config import DRIVE_FOLDER
from utils.compute_report import ComputeReport
from utils.gee_utils import export_image_to_drive, initialize_ee

# 1. Setup ----

# 1.1 User parameters ----
TASK_PREFIX = "FABDEM_US_Canada"
FILE_PREFIX = "fabdem_us_canada"

BASE_SCALE_M = 30  # FABDEM's nominal resolution (m)
EXPORT_CRS = "EPSG:4326"
PRINT_STATS = False  # min/max check (slow for large AOIs)
USE_TEST_AOI = True  # True: small test AOI; False: US+Canada
COMPUTE_REPORT = True  # write EECU usage report (txt)
# Block until every export task finishes so its batch
# EECU-seconds land in the compute report. Costs the full
# export runtime (hours for a province-wide run), so keep it
# False for production runs and turn it on when profiling a
# test AOI.
WAIT_FOR_EXPORTS = False

# 1.2 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.3 Set up run bookkeeping ----
# report profiles compute usage; export_tasks collects the
# export tasks it logs.
report = ComputeReport(FILE_PREFIX, enabled=COMPUTE_REPORT)
export_tasks = []

# 2. Define study area ----
# This script covers the US and Canada, not Alberta, so it
# does not use define_study_area or the ABMI reference grid.
# Clipping directly to the feature collection avoids geometry
# edge limits.

if USE_TEST_AOI:
    # Small aoi for testing purposes
    aoi = ee.Geometry.Polygon([
        [-113.5, 55.5],  # Top-left corner
        [-113.5, 55.0],  # Bottom-left corner
        [-112.8, 55.0],  # Bottom-right corner
        [-112.8, 55.5],  # Top-right corner
    ])
else:
    # US and Canada country features
    aoi = ee.FeatureCollection("FAO/GAUL/2015/level0").filter(
        ee.Filter.Or(
            ee.Filter.eq("ADM0_NAME", "Canada"),
            ee.Filter.eq(
                "ADM0_NAME", "United States of America"
            ),
        )
    )

# The export region and stats reducer both need an
# ee.Geometry; aoi is a FeatureCollection when USE_TEST_AOI
# is False.
aoi_geom = aoi if USE_TEST_AOI else aoi.geometry()

# 3. Build the layer ----
# Mosaicked FABDEM clipped to the study area.
layer = (
    ee.ImageCollection("projects/sat-io/open-datasets/FABDEM")
    .mosaic()
    .setDefaultProjection("EPSG:3402", None, BASE_SCALE_M)
    .clip(aoi)
)

# 3.1 Check layer values (optional) ----
# Earth Engine is lazy, so the profiler needs an evaluated
# computation to measure EECU usage; this also runs when
# COMPUTE_REPORT is on.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("Elevation min/max (reduceRegion)"):
        stats = layer.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=aoi_geom,
            scale=BASE_SCALE_M,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print("Elevation min/max:", stats)

# 4. Export ----
# Monitor progress at
# https://code.earthengine.google.com/tasks
export_tasks.append(
    export_image_to_drive(
        image=layer,
        description=TASK_PREFIX,
        region=aoi_geom,
        folder=DRIVE_FOLDER,
        file_name_prefix=FILE_PREFIX,
        scale=BASE_SCALE_M,
        crs=EXPORT_CRS,
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
