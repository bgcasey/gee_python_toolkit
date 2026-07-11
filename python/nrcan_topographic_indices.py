# ---
# title:   NRCan/CDEM Terrain Metrics
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - NRCan/CDEM DEM ImageCollection
#   - FAO GAUL province boundaries
# outputs:
#   - Terrain metrics GeoTIFF for Alberta (exported to
#     Google Drive): elevation, slope, aspect, northness,
#     eastness
# notes:
#   This script calculates terrain metrics (slope, aspect,
#   northness, eastness) from the NRCan/CDEM dataset. The
#   DEM collection is mosaicked and clipped to the AOI,
#   metrics are combined into a multiband image, and the
#   result is exported to Google Drive. Map visualization
#   layers from the original GEE JavaScript are dropped.
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

from _gee_config import DRIVE_FOLDER
from utils.compute_report import ComputeReport
from utils.gee_utils import export_image_to_drive, initialize_ee

# 1. Setup ----

# 1.1 User parameters ----
EXPORT_SCALE = 30  # meters
EXPORT_CRS = "EPSG:3348"  # NRCan/CDEM native projection
PRINT_STATS = True  # min/max check (slow for large AOIs)
USE_TEST_AOI = True  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt);
# blocks until the export task finishes

# 1.2 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.3 Set up compute usage report ----
# Profiles EECU usage per section and per export task.
# Best used with USE_TEST_AOI = True to find choke
# points cheaply before a full-province run.
report = ComputeReport(
    "nrcan_topographic_indices",
    out_dir=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "gee_compute_reports",
    ),
    enabled=COMPUTE_REPORT,
)

# 2. Define study area ----
# This section defines the export geometry. It uses a
# small test polygon when USE_TEST_AOI is True; otherwise
# it filters the FAO GAUL provinces for Alberta.

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

# 3. Terrain metrics calculation ----
# This section mosaics and clips the DEM, then derives
# slope, aspect, and northness. Northness converts aspect
# (degrees) to radians and takes the cosine. It produces a
# single multiband terrain image.

# 3.1 Load and mosaic DEM ----
dem = (
    ee.ImageCollection("NRCan/CDEM")
    .mosaic()
    .clip(aoi)
    .toFloat()
    .setDefaultProjection("EPSG:3348", None, 23.19)
)

# 3.2 Slope and aspect ----
slope = ee.Terrain.slope(dem).rename("slope")
aspect = ee.Terrain.aspect(dem).rename("aspect")

# 3.3 Northness ----
# Convert aspect degrees to radians before taking cosine
northness = (
    aspect.multiply(math.pi)
    .divide(180)
    .cos()
    .rename("northness")
    .toFloat()
)

# 3.4 Combine terrain metrics ----
# Band order matches the original script: elevation, slope,
# northness, aspect.
terrain = (
    dem.rename("elevation")
    .addBands(slope)
    .addBands(northness)
    .addBands(aspect)
)

# 3.5 Check min and max values (optional) ----
# Also runs when COMPUTE_REPORT is on: Earth Engine is
# lazy, so the profiler needs an evaluated computation
# (getInfo) to measure per-algorithm EECU usage.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("Terrain min/max (reduceRegion)"):
        stats = terrain.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=aoi,
            scale=EXPORT_SCALE,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print("Terrain min and max values:", stats)

# 4. Export data ----
# This section exports the terrain image to Google Drive
# as a GeoTIFF. Set wait=True to block until the task
# finishes; otherwise monitor progress at
# https://code.earthengine.google.com/tasks

task = export_image_to_drive(
    image=terrain,
    description="terrain_metrics_export",
    region=aoi,
    folder=DRIVE_FOLDER,
    file_name_prefix="terrain_metrics",
    scale=EXPORT_SCALE,
    crs=EXPORT_CRS,
    max_pixels=1e13,
    wait=False,
)

# 5. Compute usage report ----
# This section waits for the export to finish, records
# its total EECU-seconds, and writes the txt report to
# gee_compute_reports/. Note: a full-province export can
# take hours; for a quick profile use the test AOI.

report.log_task(task)
report.write()

# End of script ----
