# ---
# title:   Hydrologically Adjusted Elevations (HAND)
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - MERIT Hydro (MERIT/Hydro/v1_0_1)
#   - FAO GAUL province boundaries
# outputs:
#   - HAND GeoTIFF for Alberta (exported to Google Drive)
# notes:
#   This script extracts the hydrologically adjusted
#   elevations (Height Above Nearest Drainage - HAND) from
#   the MERIT Hydro dataset, clips to the AOI, and exports
#   the HAND band to Google Drive. Map visualization layers
#   from the original GEE JavaScript are dropped.
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
EXPORT_SCALE = 92.77  # meters (MERIT Hydro native ~90 m)
EXPORT_CRS = "EPSG:4326"
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
    "hydrologically_adjusted_elevation",
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

# 3. Extract HAND band ----
# This section clips the MERIT Hydro image to the AOI and
# selects the 'hnd' (Height Above Nearest Drainage) band,
# renaming it to 'HAND'. It produces a single-band image.

merit_hydro = ee.Image("MERIT/Hydro/v1_0_1")
hand = merit_hydro.clip(aoi).select("hnd").rename("HAND")

# 3.1 Check min and max values (optional) ----
# Also runs when COMPUTE_REPORT is on: Earth Engine is
# lazy, so the profiler needs an evaluated computation
# (getInfo) to measure per-algorithm EECU usage.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("HAND min/max (reduceRegion)"):
        stats = hand.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=aoi,
            scale=EXPORT_SCALE,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print("HAND min and max values:", stats)

# 4. Export data ----
# This section exports the HAND image to Google Drive as a
# GeoTIFF. Set wait=True to block until the task finishes;
# otherwise monitor progress at
# https://code.earthengine.google.com/tasks

task = export_image_to_drive(
    image=hand,
    description="HAND_Export",
    region=aoi,
    folder=DRIVE_FOLDER,
    file_name_prefix="hydrologically_adjusted_elevations",
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
