# ---
# title:   FABDEM TRI (Terrain Ruggedness Index) for Alberta
# author:  Brendan Casey
# created: 2026-07-11
# inputs:
#   - FABDEM ImageCollection
#     (projects/sat-io/open-datasets/FABDEM)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - TRI GeoTIFF for Alberta, either aligned to the ABMI
#     1 km reference grid or at BASE_SCALE_M resolution,
#     selected by EXPORT_TARGET (exported to Google Drive)
# notes:
#   Terrain Ruggedness Index from the FABDEM bare-earth DEM
#   (30 m, forests and buildings removed), following Riley et
#   al. (1999):
#
#     TRI = sqrt( sum( (z_i - z_0)^2 ) )
#
#   where z_0 is the centre cell and z_i are the cells in a
#   surrounding neighbourhood (classically the eight cells of
#   a 3x3 window). TRI is the root-summed-squared elevation
#   difference between a cell and its neighbours; it is high
#   in rugged terrain and near zero on smooth surfaces, in
#   the same units as elevation (metres).
#
#   At BASE_SCALE_M = 50 the 3x3 window is ~150 m across.
#   Ruggedness is scale-dependent, so this is a coarser
#   measure than a 30 m TRI.
#
#   Data citations:
#   Hawker, L., et al. (2022). A 30 m global map of
#   elevation with forests and buildings removed.
#   Environmental Research Letters, 17(2), 024016.
#   doi:10.1088/1748-9326/ac4d4f
#
#   Riley, S. J., DeGloria, S. D., & Elliot, R. (1999). A
#   terrain ruggedness index that quantifies topographic
#   heterogeneity. Intermountain Journal of Sciences,
#   5(1-4), 23-27.
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
    to_reference_grid,
)

# 1. Setup ----

# 1.1 User parameters ----
BAND_NAME = "tri"
TASK_PREFIX = "FABDEM_TRI_Alberta"
FILE_PREFIX = "fabdem_tri_alberta"

# Resolution TRI is computed at before aggregating to 1 km.
BASE_SCALE_M = 40
TRI_WINDOW_RADIUS = 1  # pixels; 1 = classic 3x3 Riley window

FOCAL_REACH_M = BASE_SCALE_M * TRI_WINDOW_RADIUS

# Reach of the nearest-neighbour gap fill applied after
# aggregation, in 1 km cells; 0 disables it.
FILL_GAPS_PX = 0

# "native" skips the aggregation and writes BASE_SCALE_M
# pixels in the grid CRS, ungridded - useful for inspecting
# the input to the aggregation.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"

USE_TEST_AOI = False  # True: small test AOI; False: Alberta
PRINT_STATS = True  # value preview (slow for large AOIs)
COMPUTE_REPORT = True  # write EECU usage report (txt)
# Costs the full export runtime (hours province-wide), so
# turn it on only when profiling the test AOI.
WAIT_FOR_EXPORTS = False

# Derived: keeps test-AOI tasks and files
# distinguishable from full-extent runs.
TEST_SUFFIX = "_test" if USE_TEST_AOI else ""

# 1.2 Validate parameters ----
if EXPORT_TARGET not in ("native", "reference_grid"):
    raise ValueError(
        "Unknown EXPORT_TARGET: "
        f"{EXPORT_TARGET!r} (use 'native' or 'reference_grid')"
    )
if FILL_GAPS_PX < 0:
    raise ValueError(
        f"FILL_GAPS_PX must be >= 0, got {FILL_GAPS_PX!r}"
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

# 3. Build the layer ----
# Reads over aoi_compute, never aoi; clipped back in section 4.
elevation = fabdem_elevation(aoi_compute, base_m=BASE_SCALE_M)

# normalize=False keeps each weight at 1, so the sum reducer
# returns true sums rather than means.
kernel = ee.Kernel.square(
    radius=TRI_WINDOW_RADIUS, units="pixels", normalize=False
)

sum_z = elevation.reduceNeighborhood(
    reducer=ee.Reducer.sum(), kernel=kernel
)
sum_z2 = elevation.multiply(elevation).reduceNeighborhood(
    reducer=ee.Reducer.sum(), kernel=kernel
)
count = elevation.reduceNeighborhood(
    reducer=ee.Reducer.count(), kernel=kernel
)

# sum((z_i - z_0)^2) = sum(z_i^2) - 2*z_0*sum(z_i) + N*z_0^2.
# max(0) guards against tiny negative values from floating
# point on flat terrain before the square root.
ssd = (
    sum_z2.subtract(elevation.multiply(sum_z).multiply(2))
    .add(elevation.multiply(elevation).multiply(count))
)
tri = ssd.max(0).sqrt().rename(BAND_NAME)

# The 1 km product is aggregated here rather than at export so
# the preview below reads the exported values; the export then
# passes aggregate=False. TRI is continuous metres, so no
# rounding.
if EXPORT_TARGET == "reference_grid":
    layer = to_reference_grid(
        tri, aoi, agg_max_pixels=AGG_MAX_PIXELS
    )
    stats_scale = COARSE_SCALE
else:
    layer = tri.clip(aoi)
    stats_scale = BASE_SCALE_M

# 3.1 Check layer values (optional) ----
# Earth Engine is lazy, so the profiler needs an evaluated
# computation to measure EECU usage; this also runs when
# COMPUTE_REPORT is on.
if PRINT_STATS or COMPUTE_REPORT:
    stats = None
    with report.section(
        f"{BAND_NAME} min/max (reduceRegion)",
        raise_on_error=False,
    ):
        stats = layer.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=aoi,
            scale=stats_scale,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print(f"{BAND_NAME} min/max:", stats)

# 4. Export ----
# layer is already built for the chosen EXPORT_TARGET, so the
# grid path passes aggregate=False. Monitor progress at
# https://code.earthengine.google.com/tasks
target_suffix = (
    "abmi1km" if EXPORT_TARGET == "reference_grid" else "native"
) + TEST_SUFFIX

if EXPORT_TARGET == "reference_grid":
    export_tasks.append(
        export_to_reference_grid(
            image=layer,
            aoi=aoi,
            description=f"{TASK_PREFIX}_{target_suffix}",
            folder=DRIVE_FOLDER,
            file_name_prefix=f"{FILE_PREFIX}_{target_suffix}",
            aggregate=False,
            fill_gaps_px=FILL_GAPS_PX,
            wait=False,
        )
    )
else:
    export_tasks.append(
        export_image_to_drive(
            image=layer,
            description=f"{TASK_PREFIX}_{target_suffix}",
            region=aoi,
            folder=DRIVE_FOLDER,
            file_name_prefix=f"{FILE_PREFIX}_{target_suffix}",
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
