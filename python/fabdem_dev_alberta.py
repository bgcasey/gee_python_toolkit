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
#     BASE_SCALE_M resolution, selected by EXPORT_TARGET
#     (exported to Google Drive)
# notes:
#   DEV (deviation from mean elevation) from the FABDEM
#   bare-earth DEM (30 m, forests and buildings removed),
#   following De Reu et al. (2013):
#
#     DEV = (z - mean_z) / SD_z
#
#   where mean_z and SD_z are the mean and standard deviation
#   of elevation within a focal window. The numerator is the
#   Topographic Position Index (TPI); dividing by SD_z
#   standardizes it by local relief, so DEV is expressed in
#   standard-deviation units rather than metres. A 5 m rise on
#   flat muskeg and a 300 m ridge in the Rockies can then score
#   alike, and values are comparable across radii by
#   construction. One export per radius in DEV_RADII.
#
#   DEV is computed at BASE_SCALE_M rather than the native
#   30 m so the 1 km aggregation stays under Earth Engine's
#   per-tile reprojection limit.
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

import math
import os
import sys

import ee

# Make utils importable regardless of the working
# directory VS Code runs the script from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gee_config import COARSE_SCALE, DRIVE_FOLDER, GRID_CRS
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
TASK_PREFIX = "FABDEM_DEV_Alberta"
FILE_PREFIX = "fabdem_dev_alberta"

# Resolution DEV is computed at before aggregating to 1 km.
# Keep it <= ~min(DEV_RADII) / 10 so the focal window is well
# resolved.
BASE_SCALE_M = 50

DEV_RADII = [250, 1000, 2000]  # one export per radius
DEV_WINDOW_SHAPE = "circle"  # "circle" or "square"
DEV_UNITS = "meters"  # "meters" or "pixels"
SD_EPSILON = 0.001  # floor for SD_z to avoid divide-by-zero

# The focal window reaches its radius; "pixels" radii convert
# to metres at the base scale.
FOCAL_REACH_M = max(DEV_RADII) * (
    BASE_SCALE_M if DEV_UNITS == "pixels" else 1
)

# Applied to every radius. "native" skips the aggregation and
# writes BASE_SCALE_M pixels in the grid CRS, ungridded -
# useful for inspecting the input to the aggregation.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"

USE_TEST_AOI = True  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt)
# Costs the full export runtime (hours province-wide), so
# turn it on only when profiling the test AOI.
WAIT_FOR_EXPORTS = False

# 1.2 Validate parameters ----
if EXPORT_TARGET not in ("native", "reference_grid"):
    raise ValueError(
        "Unknown EXPORT_TARGET: "
        f"{EXPORT_TARGET!r} (use 'native' or 'reference_grid')"
    )
if DEV_WINDOW_SHAPE not in ("circle", "square"):
    raise ValueError(
        f"Unsupported window shape: {DEV_WINDOW_SHAPE!r}"
    )

# 1.3 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.4 Derived scales ----
# AGG_BUFFER_M is 2x the output scale, so a 1 km cell touching
# the aoi at a corner still reads a full diagonal (1414 m)
# beyond it; AGG_MAX_PIXELS is reduceResolution's per-cell
# input budget, +2 covering a cell that straddles a base pixel
# on each side.
AGG_BUFFER_M = 2 * (
    COARSE_SCALE
    if EXPORT_TARGET == "reference_grid"
    else BASE_SCALE_M
)
COMPUTE_BUFFER_M = max(FOCAL_REACH_M, AGG_BUFFER_M)
AGG_MAX_PIXELS = (
    math.ceil(COARSE_SCALE / BASE_SCALE_M) + 2
) ** 2

# 1.5 Set up run bookkeeping ----
# report profiles compute usage; export_tasks collects the
# export tasks it logs.
report = ComputeReport(FILE_PREFIX, enabled=COMPUTE_REPORT)
export_tasks = []

# 2. Define study area ----
# aoi is the export / crop boundary; aoi_compute adds the ring
# the source is read over.
aoi, aoi_compute = define_study_area(
    use_test_aoi=USE_TEST_AOI,
    buffer_m=COMPUTE_BUFFER_M,
)

# 3. Prepare the DEM ----
# One elevation image feeds every focal radius below.
elevation = fabdem_elevation(aoi_compute, base_m=BASE_SCALE_M)

# 4. Compute, aggregate, and export DEV per focal radius ----
# Larger radii use bigger kernels and cost proportionally more
# compute; the per-task batch EECU-seconds in the report show
# where. DEV is unitless and continuous, so it is not rounded.
# Monitor progress at
# https://code.earthengine.google.com/tasks
target_suffix = (
    "abmi1km" if EXPORT_TARGET == "reference_grid" else "native"
)

for radius in DEV_RADII:
    # One kernel drives both reductions, so numerator and
    # denominator share the exact same window.
    kernel = (
        ee.Kernel.circle(radius=radius, units=DEV_UNITS)
        if DEV_WINDOW_SHAPE == "circle"
        else ee.Kernel.square(radius=radius, units=DEV_UNITS)
    )

    mean_z = elevation.reduceNeighborhood(
        reducer=ee.Reducer.mean(), kernel=kernel
    )
    sd_z = elevation.reduceNeighborhood(
        reducer=ee.Reducer.stdDev(), kernel=kernel
    )

    # Floor SD at SD_EPSILON so near-flat windows do not blow
    # up the ratio.
    dev = (
        elevation.subtract(mean_z)
        .divide(sd_z.max(SD_EPSILON))
        .rename(f"dev_{radius}")
    )

    description = f"{TASK_PREFIX}_r{radius}_{target_suffix}"
    file_name_prefix = f"{FILE_PREFIX}_r{radius}_{target_suffix}"

    if EXPORT_TARGET == "reference_grid":
        export_tasks.append(
            export_to_reference_grid(
                image=dev,
                aoi=aoi,
                description=description,
                folder=DRIVE_FOLDER,
                file_name_prefix=file_name_prefix,
                agg_max_pixels=AGG_MAX_PIXELS,
                wait=False,
            )
        )
    else:
        export_tasks.append(
            export_image_to_drive(
                image=dev.clip(aoi),
                description=description,
                region=aoi,
                folder=DRIVE_FOLDER,
                file_name_prefix=file_name_prefix,
                scale=BASE_SCALE_M,
                crs=GRID_CRS,
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
