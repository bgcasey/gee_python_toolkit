# ---
# title:   FABDEM Slope for Alberta
# author:  Brendan Casey
# created: 2026-07-11
# inputs:
#   - FABDEM ImageCollection
#     (projects/sat-io/open-datasets/FABDEM)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - Slope GeoTIFF for Alberta, either aligned to the ABMI 1 km
#     reference grid or at FOCAL_BASE_M resolution,
#     selected by EXPORT_TARGET (exported to Google Drive)
# notes:
#   This script calculates slope (in degrees) from the FABDEM
#   bare-earth DEM (30 m, forests and buildings removed). The
#   collection is mosaicked and pinned to a FOCAL_BASE_M metric
#   projection, slope is computed, then aggregated to the ABMI
#   1 km reference grid and exported so it stacks with the other
#   grid layers.
#
#   All of the grid / boundary / aggregation / export plumbing
#   lives in utils/gee_utils.py (define_study_area,
#   fabdem_elevation, export_to_reference_grid) and the grid
#   constants in _gee_config.py, so this script only holds the
#   slope-specific bits. See fabdem_tpi_alberta.py for the same
#   pattern with a focal calculation.
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

import ee

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
# FOCAL_BASE_M is the resolution slope is computed at before
# aggregating to 1 km (>= ~50 m for a full-province run; see
# to_reference_grid). Slope uses a 3x3 neighborhood, so the
# compute buffer only needs about one base pixel.
FOCAL_BASE_M = 50
# Raster export target. "reference_grid" aggregates slope (area
# mean) onto the ABMI 1 km grid so it stacks with the other 1 km
# covariates. "native" skips the aggregation and exports the
# image at FOCAL_BASE_M in the grid CRS, ungridded - useful for
# inspecting the input to the aggregation.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"

# Compute ring grown around the aoi before the source is
# clipped, sized at 2x the output scale. Every output pixel -
# a 1 km grid cell or a native pixel - is then built from a
# full neighbourhood rather than one truncated at the aoi
# edge; a 1 km cell can touch the aoi at a corner and still
# reach a full diagonal (1414 m) beyond it. The exported
# image is clipped back to the plain aoi, so the ring never
# widens the output.
COARSE_SCALE = 1000  # ABMI reference grid cell (m)
AGG_BUFFER_M = 2 * (
    COARSE_SCALE
    if EXPORT_TARGET == "reference_grid"
    else FOCAL_BASE_M
)
BUFFER_MAX_ERROR_M = 100
USE_TEST_AOI = True  # True: small test AOI; False: Alberta
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

# 1.3 Set up compute usage report ----
# Records total EECU-seconds for each export task.
# Best used with USE_TEST_AOI = True to gauge compute
# cost cheaply before a full-province run.
report = ComputeReport(
    "fabdem_slope_alberta",
    enabled=COMPUTE_REPORT,
)

# 2. Define study area ----
# aoi is the export / crop boundary; aoi_compute adds a one-
# pixel ring so the 3x3 slope kernel is unbiased at the edge.
aoi, aoi_compute = define_study_area(
    use_test_aoi=USE_TEST_AOI,
    buffer_m=max(FOCAL_BASE_M, AGG_BUFFER_M),
)

# 3. Slope calculation ----
# FABDEM at the FOCAL_BASE_M base projection, then slope in
# degrees. It produces a single-band slope image.
elevation = fabdem_elevation(aoi_compute, base_m=FOCAL_BASE_M)
slope = ee.Terrain.slope(elevation).rename("slope")

# 4. Aggregate to the grid and export ----
# EXPORT_TARGET picks the output. "reference_grid" hands slope to
# export_to_reference_grid, which aggregates it to the 1 km ABMI
# grid by area mean and exports it on the grid's exact CRS and
# transform. "native" skips the aggregation and writes slope at
# FOCAL_BASE_M in the grid CRS instead. Set wait=True to block;
# otherwise monitor progress at
# https://code.earthengine.google.com/tasks
if EXPORT_TARGET == "reference_grid":
    task = export_to_reference_grid(
        image=slope,
        aoi=aoi,
        description="FABDEM_Slope_Alberta_abmi1km",
        folder=DRIVE_FOLDER,
        file_name_prefix="fabdem_slope_alberta_abmi1km",
        wait=False,
    )
elif EXPORT_TARGET == "native":
    task = export_image_to_drive(
        image=slope.clip(aoi),
        description="FABDEM_Slope_Alberta_native",
        region=aoi,
        folder=DRIVE_FOLDER,
        file_name_prefix="fabdem_slope_alberta_native",
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

if WAIT_FOR_EXPORTS:
    report.log_task(task)
report.write()

# End of script ----
