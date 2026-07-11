# ---
# title:   Download FABDEM for US and Canada
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
EXPORT_SCALE = 100  # meters
EXPORT_CRS = "EPSG:4326"
PRINT_STATS = True  # min/max check (slow for large AOIs)
USE_TEST_AOI = False  # True: small test AOI; False: US+Canada
COMPUTE_REPORT = True  # write EECU usage report (txt);
# blocks until the export task finishes

# 1.2 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.3 Set up compute usage report ----
# Profiles EECU usage per section and per export task.
# Best used with USE_TEST_AOI = True to find choke
# points cheaply before a full US + Canada run.
report = ComputeReport(
    "fabdem",
    out_dir=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "gee_compute_reports",
    ),
    enabled=COMPUTE_REPORT,
)

# 2. Define study area ----
# This section defines the export geometry. It uses a
# small test polygon when USE_TEST_AOI is True; otherwise
# it filters the FAO GAUL countries for the US and Canada.
# Clipping directly to the feature collection avoids
# geometry edge limits.

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
# ee.Geometry; derive one from the AOI (which may be a
# FeatureCollection when USE_TEST_AOI is False).
if USE_TEST_AOI:
    aoi_geom = aoi
else:
    aoi_geom = aoi.geometry()

# 3. Create and clip elevation mosaic ----
# This section mosaics the FABDEM collection and clips it
# to the study area. It produces a clipped elevation image.

fabdem = ee.ImageCollection(
    "projects/sat-io/open-datasets/FABDEM"
)
elev = fabdem.mosaic().setDefaultProjection(
    "EPSG:3857", None, 30
)
elev_clipped = elev.clip(aoi)

# 3.1 Check min and max values (optional) ----
# Also runs when COMPUTE_REPORT is on: Earth Engine is
# lazy, so the profiler needs an evaluated computation
# (getInfo) to measure per-algorithm EECU usage.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("Elevation min/max (reduceRegion)"):
        stats = elev_clipped.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=aoi_geom,
            scale=EXPORT_SCALE,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print("Elevation min and max values:", stats)

# 4. Export data ----
# This section exports the clipped elevation image to
# Google Drive as a GeoTIFF. Set wait=True to block until
# the task finishes; otherwise monitor progress at
# https://code.earthengine.google.com/tasks

task = export_image_to_drive(
    image=elev_clipped,
    description="FABDEM_US_Canada",
    region=aoi_geom,
    folder=DRIVE_FOLDER,
    file_name_prefix="fabdem_us_canada",
    scale=EXPORT_SCALE,
    crs=EXPORT_CRS,
    max_pixels=1e13,
    wait=False,
)

# 5. Compute usage report ----
# This section waits for the export to finish, records
# its total EECU-seconds, and writes the txt report to
# gee_compute_reports/. Note: a full US + Canada export
# can take hours; for a quick profile use the test AOI.

report.log_task(task)
report.write()

# End of script ----
