# ---
# title:   FABDEM Elevation for Alberta
# author:  Brendan Casey
# created: 2026-08-05
# inputs:
#   - FABDEM ImageCollection
#     (projects/sat-io/open-datasets/FABDEM)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - Mean-elevation GeoTIFF for Alberta, either aligned to the ABMI
#     1 km reference grid or at FOCAL_BASE_M resolution,
#     selected by EXPORT_TARGET (exported to Google Drive)
# notes:
#   This script produces mean elevation per 1 km cell from the
#   FABDEM bare-earth DEM (30 m, forests and buildings
#   removed), as a covariate that stacks with the other ABMI
#   1 km terrain layers (slope, TPI, DEV, TRI, TWI).
#
#   Elevation is a raw value (no neighborhood operation), so
#   there is no edge bias and no compute buffer is needed. It
#   is served at FOCAL_BASE_M from FABDEM's pyramids (an area
#   mean of the native 30 m), then aggregated to the 1 km grid
#   by area mean -- i.e. mean elevation per cell. The grid /
#   boundary / aggregation / export plumbing lives in
#   utils/gee_utils.py.
#
#   For the full-resolution source DEM over a larger extent,
#   see fabdem.py (US + Canada, native 30 m, not grid-aligned).
#
#   Data citation:
#   Hawker, L., et al. (2022). A 30 m global map of
#   elevation with forests and buildings removed.
#   Environmental Research Letters, 17(2), 024016.
#   doi:10.1088/1748-9326/ac4d4f
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

# Make utils importable regardless of the working
# directory VS Code runs the script from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gee_config import DRIVE_FOLDER, GRID_CRS
from utils.compute_report import ComputeReport
from utils.gee_utils import (
    define_study_area,
    export_image_to_drive,
    export_to_reference_grid,
    fabdem_elevation,
    initialize_ee,
)

# 1. Setup ----

# 1.1 User parameters ----
# FOCAL_BASE_M is the resolution elevation is served at before
# aggregating to 1 km. Elevation has no neighborhood op, so any
# base under the reprojection limit works (>= ~50 m for a
# full-province run); 50 m is the area mean of the native 30 m.
FOCAL_BASE_M = 50
# Raster export target. "reference_grid" aggregates elevation (area
# mean) onto the ABMI 1 km grid so it stacks with the other 1 km
# covariates. "native" skips the aggregation and exports the
# image at FOCAL_BASE_M in the grid CRS, ungridded - useful for
# inspecting the input to the aggregation.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"
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
    "fabdem_elev_alberta",
    enabled=COMPUTE_REPORT,
)

# 2. Define study area ----
# No neighborhood op, so no compute ring is needed (buffer_m=0
# makes aoi_compute == aoi).
aoi, aoi_compute = define_study_area(
    use_test_aoi=USE_TEST_AOI,
    buffer_m=0,
)

# 3. Prepare the DEM ----
# FABDEM served at the FOCAL_BASE_M base projection (pyramid
# area mean of the native 30 m) and clipped to the AOI.
elevation = fabdem_elevation(aoi_compute, base_m=FOCAL_BASE_M)

# 4. Aggregate to the grid and export ----
# EXPORT_TARGET picks the output. "reference_grid" hands elevation to
# export_to_reference_grid, which aggregates it to the 1 km ABMI
# grid by area mean and exports it on the grid's exact CRS and
# transform. "native" skips the aggregation and writes elevation at
# FOCAL_BASE_M in the grid CRS instead. Set wait=True to block;
# otherwise monitor progress at
# https://code.earthengine.google.com/tasks
if EXPORT_TARGET == "reference_grid":
    task = export_to_reference_grid(
        image=elevation.rename("elevation"),
        aoi=aoi,
        description="FABDEM_Elevation_Alberta_1km",
        folder=DRIVE_FOLDER,
        file_name_prefix="fabdem_elevation_alberta_1km",
        wait=False,
    )
elif EXPORT_TARGET == "native":
    task = export_image_to_drive(
        image=elevation.rename("elevation"),
        description=f"FABDEM_Elevation_Alberta_{FOCAL_BASE_M}m",
        region=aoi,
        folder=DRIVE_FOLDER,
        file_name_prefix=f"fabdem_elevation_alberta_{FOCAL_BASE_M}m",
        scale=FOCAL_BASE_M,
        crs=GRID_CRS,
        max_pixels=1e13,
        wait=False,
    )
else:
    raise ValueError(
        "Unknown EXPORT_TARGET: "
        f"{EXPORT_TARGET!r} (use 'native' or 'reference_grid')"
    )

# 5. Compute usage report ----
# This section waits for the export to finish, records
# its total EECU-seconds, and writes the txt report to
# gee_compute_reports/. Note: a full-province export can
# take hours; for a quick profile use the test AOI.

report.log_task(task)
report.write()

# End of script ----
