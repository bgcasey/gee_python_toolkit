# ---
# title:   FABDEM Topographic Wetness Index for Alberta
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - FABDEM ImageCollection
#     (projects/sat-io/open-datasets/FABDEM)
#   - MERIT Hydro upslope area (MERIT/Hydro/v1_0_1)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - TWI GeoTIFF for Alberta, aligned to the ABMI 1 km
#     reference grid (exported to Google Drive)
# notes:
#   This script calculates the Topographic Wetness Index
#   (TWI) as ln(a / tan(b)), where a is upslope drainage
#   area (m^2) and b is slope (radians).
#
#   FABDEM is a bare-earth DEM with no flow-accumulation
#   band, and Earth Engine has no native flow-accumulation
#   algorithm. This script therefore uses a hybrid: slope
#   from FABDEM (computed at FOCAL_BASE_M) and upslope area
#   from MERIT Hydro 'upa' (~90 m). The result is aggregated
#   to the ABMI 1 km reference grid and exported so it stacks
#   with the other grid layers.
#
#   Slope is computed at FOCAL_BASE_M (50 m) rather than the
#   native 30 m so the 1 km aggregation stays under Earth
#   Engine's per-tile reprojection limit (see
#   utils.gee_utils.to_reference_grid). The grid / boundary /
#   aggregation / export plumbing lives in utils/gee_utils.py.
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
from utils.gee_utils import (
    define_study_area,
    export_to_reference_grid,
    fabdem_elevation,
    initialize_ee,
)

# 1. Setup ----

# 1.1 User parameters ----
# FOCAL_BASE_M is the resolution slope is computed at before
# aggregating to 1 km (>= ~50 m for a full-province run). Slope
# uses a 3x3 neighborhood, so the compute buffer only needs
# about one base pixel.
FOCAL_BASE_M = 50
USE_TEST_AOI = True  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt);
# blocks until the export task finishes

# 1.2 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.3 Set up compute usage report ----
# Records total EECU-seconds for each export task.
# Best used with USE_TEST_AOI = True to gauge compute
# cost cheaply before a full-province run.
report = ComputeReport(
    "fabdem_twi_alberta",
    enabled=COMPUTE_REPORT,
)

# 2. Define study area ----
# aoi is the export / crop boundary; aoi_compute adds a one-
# pixel ring so the 3x3 slope kernel is unbiased at the edge.
aoi, aoi_compute = define_study_area(
    use_test_aoi=USE_TEST_AOI,
    buffer_m=FOCAL_BASE_M,
)

# 3. TWI calculation ----
# TWI from FABDEM slope and MERIT Hydro upslope area. It
# produces a single-band TWI image.

# FABDEM at the FOCAL_BASE_M base projection, then slope.
elevation = fabdem_elevation(aoi_compute, base_m=FOCAL_BASE_M)
slope = ee.Terrain.slope(elevation)

# Upslope area from MERIT Hydro, converted km^2 -> m^2. MERIT
# 'upa' is a stored (pyramided) dataset, so it needs no coarse
# base; clip to the buffered AOI to match the elevation extent.
upslope_area = (
    ee.Image("MERIT/Hydro/v1_0_1")
    .select("upa")
    .clip(aoi_compute)
    .multiply(1e6)
    .rename("upslope_area")
)

# Convert slope from degrees to radians
slope_rad = slope.multiply(math.pi / 180).rename("slope_rad")

# Floor tan(b) at a small value so flat areas (slope 0) are
# not masked by division by zero
tan_b = slope_rad.tan().max(0.001)

# Calculate TWI: ln(a / tan(b)). Continuous, so no rounding.
twi = upslope_area.divide(tan_b).log().rename("twi")

# 4. Aggregate to the grid and export ----
# export_to_reference_grid aggregates TWI to the 1 km ABMI grid
# by area mean and exports it on the grid's exact CRS and
# transform. Set wait=True to block; otherwise monitor progress
# at https://code.earthengine.google.com/tasks
task = export_to_reference_grid(
    image=twi,
    aoi=aoi,
    description="FABDEM_TWI_Alberta_1km",
    folder=DRIVE_FOLDER,
    file_name_prefix="fabdem_twi_alberta_1km",
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
