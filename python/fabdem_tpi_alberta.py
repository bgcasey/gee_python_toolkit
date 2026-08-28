# ---
# title:   FABDEM TPI for Alberta
# author:  Brendan Casey
# created: 2026-07-11
# inputs:
#   - FABDEM ImageCollection
#     (projects/sat-io/open-datasets/FABDEM)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - One TPI GeoTIFF per focal radius for Alberta, either
#     aligned to the ABMI 1 km reference grid or at
#     BASE_SCALE_M resolution, selected by EXPORT_TARGET
#     (exported to Google Drive)
# notes:
#   Topographic Position Index from the FABDEM bare-earth DEM
#   (30 m, forests and buildings removed): elevation minus the
#   mean elevation of a surrounding neighbourhood, one export
#   per radius in TPI_RADII.
#
#   TPI is computed at BASE_SCALE_M rather than the native
#   30 m. Aggregating a computed 30 m layer straight to 1 km
#   over the whole province exceeds Earth Engine's per-tile
#   reprojection limit ("Reprojection output too large"):
#   filling one 1 km output tile forces the focal mean over a
#   ~256 km footprint, which at 30 m is ~8600 px per side, over
#   the ~8192 cap. A 50 m base drops it to ~5200 px. FABDEM is
#   served from its pyramids at that scale (an area mean of the
#   30 m elevation), so for radii much larger than the base the
#   coarser base costs almost nothing.
#
#   Reference grid:
#     \\ABMI-DATA2\science\spatial_data\temp\
#       GRID1SQKM_AB2020.gdb (layer Grid_1KM_revAB2020)
#     CRS EPSG:3400 (NAD83 / Alberta 10-TM Forest),
#     1000 m cells, registered at x = 616.1822 + k*1000 and
#     y = 532.4311 + m*1000 (metres).
#
#   Data citations:
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

import math
import os
import sys

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
TASK_PREFIX = "FABDEM_TPI_Alberta"
FILE_PREFIX = "fabdem_tpi_alberta"

# Resolution TPI is computed at before aggregating to 1 km.
BASE_SCALE_M = 40

TPI_RADII = [250]  # one export per radius
TPI_WINDOW_SHAPE = "circle"  # "circle" or "square"
TPI_UNITS = "meters"  # "meters" or "pixels"

# The focal mean reaches its radius; "pixels" radii convert to
# metres at the base scale.
FOCAL_REACH_M = max(TPI_RADII) * (
    BASE_SCALE_M if TPI_UNITS == "pixels" else 1
)

# Applied to every radius. "native" skips the aggregation and
# writes BASE_SCALE_M pixels in the grid CRS, ungridded -
# useful for inspecting the input to the aggregation.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"

USE_TEST_AOI = False  # True: small test AOI; False: Alberta
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

# 4. Compute, aggregate, and export TPI per focal radius ----
# Larger radii use bigger focal kernels and cost
# proportionally more compute; the per-task batch EECU-seconds
# in the report show where. round_values stores the 1 km
# product as integer metres, rounding after the area mean so
# the averaging keeps full precision. Monitor progress at
# https://code.earthengine.google.com/tasks
target_suffix = (
    "abmi1km" if EXPORT_TARGET == "reference_grid" else "native"
)

for radius in TPI_RADII:
    tpi = elevation.subtract(
        elevation.focalMean(radius, TPI_WINDOW_SHAPE, TPI_UNITS)
    ).rename(f"tpi_{radius}")

    description = f"{TASK_PREFIX}_r{radius}_{target_suffix}"
    file_name_prefix = f"{FILE_PREFIX}_r{radius}_{target_suffix}"

    if EXPORT_TARGET == "reference_grid":
        export_tasks.append(
            export_to_reference_grid(
                image=tpi,
                aoi=aoi,
                description=description,
                folder=DRIVE_FOLDER,
                file_name_prefix=file_name_prefix,
                agg_max_pixels=AGG_MAX_PIXELS,
                round_values=True,
                wait=False,
            )
        )
    else:
        # No rounding here: this is the raw input to the
        # aggregation, so it keeps sub-metre TPI values.
        export_tasks.append(
            export_image_to_drive(
                image=tpi.toFloat().clip(aoi),
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
