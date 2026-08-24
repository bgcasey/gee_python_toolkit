# ---
# title:   FABDEM DEV (Deviation from Mean Elevation) for Alberta
# author:  Brendan Casey
# created: 2026-07-11
# inputs:
#   - FABDEM ImageCollection
#     (projects/sat-io/open-datasets/FABDEM)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - One DEV GeoTIFF per focal radius for Alberta, either
#     aligned to the ABMI 1 km reference grid or at
#     FOCAL_BASE_M resolution, selected by EXPORT_TARGET
#     (exported to Google Drive)
# notes:
#   This script calculates DEV (deviation from mean
#   elevation) from the FABDEM bare-earth DEM (30 m, forests
#   and buildings removed), following De Reu et al. (2013):
#
#     DEV = (z - mean_z) / SD_z
#
#   where mean_z and SD_z are the mean and standard deviation
#   of elevation within a focal window. The numerator is the
#   Topographic Position Index (TPI); dividing by
#   SD_z standardizes it by local relief, so DEV is expressed
#   in standard-deviation units rather than metres. A 5 m rise
#   on flat muskeg and a 300 m ridge in the Rockies can then
#   score alike, and values are comparable across radii by
#   construction. DEV is computed at each focal radius in
#   DEV_RADII, aggregated to the 1 km reference grid, and
#   exported to Google Drive per radius.
#
#   DEV is computed at FOCAL_BASE_M (50 m) rather than the
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
#   De Reu, J., et al. (2013). Application of the topographic
#   position index to heterogeneous landscapes. Geomorphology,
#   186, 39-49. doi:10.1016/j.geomorph.2012.12.015
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
# FOCAL_BASE_M is the resolution DEV is computed at before
# aggregating to 1 km (>= ~50 m for a full-province run; keep
# it <= ~min(DEV_RADII) / 10 so the focal window is well
# resolved).
FOCAL_BASE_M = 50
DEV_RADII = [250, 1000, 2000]  # one export per radius
# Raster export target, applied to every radius. "reference_grid"
# aggregates DEV (area mean) onto the ABMI 1 km grid so it stacks
# with the other 1 km covariates. "native" skips the aggregation
# and exports DEV at FOCAL_BASE_M in the grid CRS, ungridded.
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

if EXPORT_TARGET not in ("native", "reference_grid"):
    raise ValueError(
        "Unknown EXPORT_TARGET: "
        f"{EXPORT_TARGET!r} (use 'native' or 'reference_grid')"
    )
DEV_WINDOW_SHAPE = "circle"  # "circle" or "square"
DEV_UNITS = "meters"  # "meters" or "pixels"
SD_EPSILON = 0.001  # floor for SD_z to avoid divide-by-zero
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
    "fabdem_dev_alberta",
    enabled=COMPUTE_REPORT,
)

# 2. Define study area ----
# aoi is the export / crop boundary; aoi_compute adds a ring
# (the largest focal radius) so the focal mean and SD are
# unbiased at the true AOI edge.
aoi, aoi_compute = define_study_area(
    use_test_aoi=USE_TEST_AOI,
    buffer_m=max(max(DEV_RADII), AGG_BUFFER_M),
)

# 3. Prepare the DEM ----
# FABDEM at the FOCAL_BASE_M base projection so the focal
# radius maps to real ground distance. The same elevation
# image feeds every focal radius below.
elevation = fabdem_elevation(aoi_compute, base_m=FOCAL_BASE_M)

# 4. Compute, aggregate, and export DEV per focal radius ----
# For each radius in DEV_RADII, elevation mean and SD are
# reduced over one shared kernel, DEV = (z - mean) / SD is
# formed, aggregated to the 1 km grid, and exported to Google
# Drive as its own GeoTIFF. Larger radii use bigger kernels
# and cost proportionally more compute; the per-task batch
# EECU-seconds in the report show where. Set wait=True on the
# export to block; otherwise monitor progress at
# https://code.earthengine.google.com/tasks

tasks = []
for radius in DEV_RADII:
    # One kernel drives both the mean and SD reductions, so
    # numerator and denominator share the exact same window.
    if DEV_WINDOW_SHAPE == "circle":
        kernel = ee.Kernel.circle(radius=radius, units=DEV_UNITS)
    elif DEV_WINDOW_SHAPE == "square":
        kernel = ee.Kernel.square(radius=radius, units=DEV_UNITS)
    else:
        raise ValueError(
            f"Unsupported window shape: {DEV_WINDOW_SHAPE}"
        )

    # Mean and standard deviation of elevation in the window
    mean_z = elevation.reduceNeighborhood(
        reducer=ee.Reducer.mean(),
        kernel=kernel,
    )
    sd_z = elevation.reduceNeighborhood(
        reducer=ee.Reducer.stdDev(),
        kernel=kernel,
    )

    # DEV: TPI (z - mean) standardized by local relief (SD).
    # Floor SD at SD_EPSILON so near-flat windows (SD ~ 0) do
    # not blow up the ratio. DEV is unitless (SD units) and
    # continuous, so round_values stays False on export.
    dev = (
        elevation
        .subtract(mean_z)
        .divide(sd_z.max(SD_EPSILON))
        .rename(f"dev_{radius}")
    )

    # 4.1 Export this radius at the chosen EXPORT_TARGET ----
    # "reference_grid" aggregates DEV to the 1 km grid; "native"
    # writes it at FOCAL_BASE_M in the grid CRS instead.
    # EXPORT_TARGET is validated once at the top of the script.
    if EXPORT_TARGET == "reference_grid":
        task = export_to_reference_grid(
            image=dev,
            aoi=aoi,
            description=f"FABDEM_DEV_Alberta_r{radius}_abmi1km",
            folder=DRIVE_FOLDER,
            file_name_prefix=f"fabdem_dev_alberta_r{radius}_abmi1km",
            wait=False,
        )
    else:
        task = export_image_to_drive(
            image=dev.clip(aoi),
            description=(
                f"FABDEM_DEV_Alberta_r{radius}_native"
            ),
            region=aoi,
            folder=DRIVE_FOLDER,
            file_name_prefix=(
                f"fabdem_dev_alberta_r{radius}_native"
            ),
            scale=FOCAL_BASE_M,
            crs=GRID_CRS,
            max_pixels=1e13,
            wait=False,
        )
    tasks.append(task)

# 5. Compute usage report ----
# This section waits for each export to finish, records its
# total EECU-seconds, and writes the txt report to
# gee_compute_reports/. Note: a full-province export can
# take hours; for a quick profile use the test AOI.

if WAIT_FOR_EXPORTS:
    for task in tasks:
        report.log_task(task)
report.write()

# End of script ----
