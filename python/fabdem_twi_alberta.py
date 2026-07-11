# ---
# title:   FABDEM Topographic Wetness Index for Alberta
# author:  bgcasey
# created: 2026-07-10
# inputs:
#   - FABDEM ImageCollection
#     (projects/sat-io/open-datasets/FABDEM)
#   - MERIT Hydro upslope area (MERIT/Hydro/v1_0_1)
#   - FAO GAUL province boundaries
# outputs:
#   - TWI GeoTIFF for Alberta (exported to Google Drive)
# notes:
#   This script calculates the Topographic Wetness Index
#   (TWI) as ln(a / tan(b)), where a is upslope drainage
#   area (m^2) and b is slope (radians).
#
#   FABDEM is a bare-earth DEM with no flow-accumulation
#   band, and Earth Engine has no native flow-accumulation
#   algorithm. This script therefore uses a hybrid: slope
#   from FABDEM (30 m, forests and buildings removed) and
#   upslope area from MERIT Hydro 'upa' (~90 m, resampled on
#   the fly). It exports a GeoTIFF to Google Drive.
#
#   Data citations:
#   Hawker, L., et al. (2022). A 30 m global map of
#   elevation with forests and buildings removed.
#   Environmental Research Letters, 17(2), 024016.
#   doi:10.1088/1748-9326/ac4d4f
#
#   Yamazaki, D., et al. (2019). MERIT Hydro: A
#   high-resolution global hydrography map based on latest
#   topography datasets. Water Resources Research, 55,
#   5053-5073. doi:10.1029/2019WR024873
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
    "fabdem_twi_alberta",
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

# 3. TWI calculation ----
# This section calculates TWI from FABDEM slope and MERIT
# Hydro upslope area. It produces a single-band TWI image.

# Load FABDEM, mosaic, and set a default projection so
# terrain algorithms have a fixed 30 m metric scale
# (EPSG:3402, Alberta 10-TM)
elevation = (
    ee.ImageCollection("projects/sat-io/open-datasets/FABDEM")
    .mosaic()
    .setDefaultProjection("EPSG:3402", None, 30)
    .clip(aoi)
)

# Calculate slope from FABDEM elevation
slope = ee.Terrain.slope(elevation)

# Load upslope area from MERIT Hydro and convert km^2 to m^2
upslope_area = (
    ee.Image("MERIT/Hydro/v1_0_1")
    .select("upa")
    .clip(aoi)
    .multiply(1e6)
    .rename("upslope_area")
)

# Convert slope from degrees to radians
slope_rad = slope.multiply(math.pi / 180).rename("slope_rad")

# Floor tan(b) at a small value so flat areas (slope 0) are
# not masked by division by zero
tan_b = slope_rad.tan().max(0.001)

# Calculate TWI: ln(a / tan(b))
twi = upslope_area.divide(tan_b).log().rename("twi")

# 3.1 Check min and max values (optional) ----
# Also runs when COMPUTE_REPORT is on: Earth Engine is
# lazy, so the profiler needs an evaluated computation
# (getInfo) to measure per-algorithm EECU usage.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("TWI min/max (reduceRegion)"):
        stats = twi.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=aoi,
            scale=EXPORT_SCALE,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print("TWI min and max values:", stats)

# 4. Export data ----
# This section exports the TWI image to Google Drive as a
# GeoTIFF. Set wait=True to block until the task finishes;
# otherwise monitor progress at
# https://code.earthengine.google.com/tasks

task = export_image_to_drive(
    image=twi,
    description="FABDEM_TWI_Alberta",
    region=aoi,
    folder=DRIVE_FOLDER,
    file_name_prefix="fabdem_twi_alberta",
    scale=EXPORT_SCALE,
    crs="EPSG:4326",
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
