# ---
# title:   Template: Alberta Gridded Raster Export
# author:  Brendan Casey
# created: 2026-08-27
# inputs:
#   - [EDIT] Source dataset (Earth Engine image or collection)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - [EDIT] <layer> GeoTIFF for Alberta, either aggregated to
#     the ABMI 1 km reference grid or written at BASE_SCALE_M,
#     selected by EXPORT_TARGET (exported to Google Drive)
# notes:
#   TEMPLATE - not a product script. Copy it to
#   python/<layer>_alberta.py, work through the [EDIT]
#   markers, delete this paragraph and the "Sizing" note
#   below, and add the script to the Contents table in
#   README.md.
#
#   The source is read over a compute ring (aoi buffered by
#   COMPUTE_BUFFER_M) and the result clipped back to aoi, so
#   boundary pixels are built on a full neighbourhood without
#   widening the exported footprint.
#
#   Key details:
#
#   - Earth Engine writes masked pixels in an integer GeoTIFF
#     as 0 with no nodata flag, so they read downstream as
#     valid zeros. float32 makes them NaN -> NA; class codes
#     survive it exactly.
#   - A 30 m base fails a full-province run with
#     "Reprojection output too large": filling one 1 km output
#     tile forces the base compute over ~256 km, about 8600 px
#     at 30 m, past Earth Engine's ~8192 px per-tile cap.
#     50 m gives ~5200 px.
#
#   Sizing BASE_SCALE_M and FOCAL_REACH_M:
#     - Pre-made source layer: FOCAL_REACH_M = 0,
#       BASE_SCALE_M = None and MIN_BASE_SCALE_M = 0 to take
#       the source's native scale.
#     - Computed focal layer (slope, TPI, TRI, DEV):
#       FOCAL_REACH_M = the largest reach in metres,
#       MIN_BASE_SCALE_M = 50, and BASE_SCALE_M no coarser
#       than ~1/10 of the smallest focal radius.
#
#   [EDIT] Data citation:
#   <Author(s). (Year). Title. Journal, vol(issue), pages.
#   doi:...>
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

import ee  # earthengine-api (version: 1.7.34)

# Make utils importable regardless of the working
# directory VS Code runs the script from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gee_config import COARSE_SCALE, DRIVE_FOLDER, GRID_CRS
from utils.compute_report import ComputeReport
from utils.gee_utils import (
    define_study_area,
    export_image_to_drive,
    export_to_reference_grid,
    initialize_ee,
    native_scale,
)

# 1. Setup ----

# 1.1 User parameters ----
# [EDIT] Output names. target_suffix ("abmi1km" or "native")
# is appended to TASK_PREFIX and FILE_PREFIX at export.
BAND_NAME = "template_layer"
TASK_PREFIX = "Template_Layer_Alberta"
FILE_PREFIX = "template_layer_alberta"

# [EDIT] Source dataset. The placeholder is FABDEM, the 30 m
# bare-earth DEM with forests and buildings removed.
SOURCE_ASSET = "projects/sat-io/open-datasets/FABDEM"

# [EDIT] Base resolution reduceResolution reads from; named
# FOCAL_BASE_M or NATIVE_SCALE / EXPORT_SCALE in the sibling
# scripts. None takes the source's native scale, floored at
# MIN_BASE_SCALE_M. See "Sizing" in the header.
BASE_SCALE_M = None
MIN_BASE_SCALE_M = 50

# [EDIT] Largest neighbourhood reach in metres; 0 when the
# layer needs no neighbours. A 3x3 kernel reaches one base
# pixel, a focal mean its radius.
FOCAL_REACH_M = 0

# [EDIT] Class codes aggregate by mode, not by mean.
IS_CATEGORICAL = False

# "native" skips the aggregation and writes BASE_SCALE_M
# pixels in the grid CRS, ungridded - useful for inspecting
# the input to the aggregation.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"

USE_TEST_AOI = True  # True: small test AOI; False: Alberta
PRINT_STATS = True  # value preview (slow for large AOIs)
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

# 1.4 Resolve the base resolution and derived scales ----
# native_scale reads the projection from an unmosaicked image,
# so source is passed to build_layer rather than a mosaic.
# AGG_BUFFER_M is 2x the output scale, so a 1 km cell touching
# the aoi at a corner still reads a full diagonal (1414 m)
# beyond it; AGG_MAX_PIXELS is reduceResolution's per-cell
# input budget, +2 covering a cell that straddles a base pixel
# on each side.
source = ee.ImageCollection(SOURCE_ASSET)
if BASE_SCALE_M is None:
    BASE_SCALE_M = native_scale(source, floor_m=MIN_BASE_SCALE_M)
print(f"Base resolution: {BASE_SCALE_M:.2f} m")

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


def build_layer(source, compute_aoi, base_m):
    """Build the layer to export, on the compute ring.

    [EDIT] Add this layer's own logic to the mosaic below. The
    returned image must stay clipped to compute_aoi rather than
    aoi, and pinned to the base projection so reduceResolution
    knows the input resolution. Drop the mosaic() call for a
    source that is a single ee.Image. Bands are renamed by the
    caller, so a multi-band layer can pass a list.
    """
    return (
        source.mosaic()
        .setDefaultProjection(GRID_CRS, None, base_m)
        .clip(compute_aoi)
    )


layer = build_layer(source, aoi_compute, BASE_SCALE_M).rename(
    BAND_NAME
)

# 3.1 Check layer values (optional) ----
# Earth Engine is lazy, so the profiler needs an evaluated
# computation to measure EECU usage; this also runs when
# COMPUTE_REPORT is on. raise_on_error=False lets a failed
# preview through without stopping the export.
if PRINT_STATS or COMPUTE_REPORT:
    stats = None
    stats_label = (
        "class frequencies" if IS_CATEGORICAL else "min/max"
    )
    with report.section(
        f"{BAND_NAME} {stats_label} (reduceRegion)",
        raise_on_error=False,
    ):
        stats = layer.reduceRegion(
            reducer=(
                ee.Reducer.frequencyHistogram()
                if IS_CATEGORICAL
                else ee.Reducer.minMax()
            ),
            geometry=aoi,
            scale=BASE_SCALE_M,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print(f"{BAND_NAME} {stats_label}:", stats)

# 4. Aggregate to the grid and export ----
# toFloat() on the native path is what keeps masked pixels
# readable as NA; export_to_reference_grid does it internally.
# Monitor progress at
# https://code.earthengine.google.com/tasks
target_suffix = (
    "abmi1km" if EXPORT_TARGET == "reference_grid" else "native"
)

if EXPORT_TARGET == "reference_grid":
    export_tasks.append(
        export_to_reference_grid(
            image=layer,
            aoi=aoi,
            description=f"{TASK_PREFIX}_{target_suffix}",
            folder=DRIVE_FOLDER,
            file_name_prefix=f"{FILE_PREFIX}_{target_suffix}",
            reducer=(
                ee.Reducer.mode() if IS_CATEGORICAL else None
            ),
            agg_max_pixels=AGG_MAX_PIXELS,
            round_values=IS_CATEGORICAL,
            wait=False,
        )
    )
else:
    export_tasks.append(
        export_image_to_drive(
            image=layer.toFloat().clip(aoi),
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
