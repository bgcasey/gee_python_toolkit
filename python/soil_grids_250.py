# ---
# title:   SoilGrids 250m v2.0 Layers Export
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - ISRIC SoilGrids 250m v2.0 mean Images
#     (projects/soilgrids-isric/*_mean)
#   - AB2020 provincial boundary (EE asset)
#   - XY points (may include locations outside Alberta)
# outputs:
#   - Multiband SoilGrids image clipped to Alberta, exported
#     either at native (~250 m, EPSG:3400) or aggregated to the
#     ABMI 1 km reference grid (EPSG:3400), selected by
#     EXPORT_TARGET.
#   - Per-batch CSVs of point-level extracted soil values for
#     ALL points including those outside the AOI (optional,
#     gated by EXTRACT_XY_POINTS).
# notes:
#   SoilGrids 250m v2.0 is a globally consistent,
#   data-driven system that predicts soil properties at
#   six standard depths (0-5, 5-15, 15-30, 30-60, 60-100,
#   100-200 cm). Each *_mean asset is a multiband Image
#   with one band per depth.
#
#   Mapped units are integer-scaled; a per-variable
#   conversion factor is applied to recover conventional
#   units:
#     bdod     (cg/cm3)     / 100 -> kg/dm3
#     cec      (mmol(c)/kg) / 10  -> cmol(c)/kg
#     cfvo     (cm3/dm3)    / 10  -> cm3/100cm3 (vol %)
#     clay     (g/kg)       / 10  -> g/100g (%)
#     nitrogen (cg/kg)      / 100 -> g/kg
#     phh2o    (pH*10)      / 10  -> pH
#     sand     (g/kg)       / 10  -> g/100g (%)
#     silt     (g/kg)       / 10  -> g/100g (%)
#     soc      (dg/kg)      / 10  -> g/kg
#     ocd      (hg/dm3)     / 10  -> kg/dm3
#     ocs      (t/ha)       / 10  -> kg/m2
#
#   The 'ocs' (organic carbon stock) asset covers only the
#   0-30 cm depth (single band). All other variables
#   retain the six-depth structure. Native band names take
#   the form '<var>_<depth>_mean' (e.g. 'clay_0-5cm_mean').
#
#   The base SoilGrids image is built unclipped (global).
#   The AOI clip is applied only for raster aggregation and
#   export so that XY point extraction can return values
#   for points located anywhere with SoilGrids coverage.
#
#   When EXPORT_TARGET is "reference_grid", aggregation to
#   1 km uses area-mean reduction and lands on the ABMI
#   reference grid via export_to_reference_grid
#   (utils.gee_utils). setDefaultProjection pins the native
#   base so reduceResolution knows the input resolution.
#
#   Citation:
#   Poggio, L., de Sousa, L. M., Batjes, N. H., Heuvelink,
#   G. B. M., Kempen, B., Ribeiro, E., and Rossiter, D.:
#   SoilGrids 2.0: producing soil information for the
#   globe with quantified spatial uncertainty, SOIL, 7,
#   217-240, https://doi.org/10.5194/soil-7-217-2021, 2021.
#
#   Setup (once):
#     pip install earthengine-api
#     earthengine authenticate
#   Then set EE_PROJECT in _gee_config.py to your
#   registered Earth Engine cloud project and run.
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
    initialize_ee,
)

# 1. Setup ----

# 1.1 User parameters ----
TASK_PREFIX = "SoilGrids_AB"
FILE_PREFIX = "soilgrids_ab"

# SoilGrids' native resolution; native_scale confirms 250.
BASE_SCALE_M = 250

# Each value is read as a stored pixel, so no neighbourhood.
FOCAL_REACH_M = 0

# Reach of the nearest-neighbour gap fill applied after
# aggregation, in 1 km cells; 0 disables it.
FILL_GAPS_PX = 20

# Raster export target. "native" exports the ~250 m SoilGrids
# image in EPSG:3400 (ungridded). "reference_grid" aggregates
# to the ABMI 1 km reference grid (EPSG:3400) by area mean via
# export_to_reference_grid, so it stacks with the FABDEM
# terrain layers. 250 m -> 1 km is only a 4x factor, well under
# Earth Engine's per-tile reprojection limit.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"


# Base path for SoilGrids 250m v2.0 assets.
BASE_PATH = "projects/soilgrids-isric/"

# Optional band filter. Set to None (or empty list) to
# keep all bands. Otherwise provide a list of band names
# in <variable>_<depth>_mean format. Any variable not
# represented is skipped at load time; remaining variables
# are loaded fully and filtered after.
SELECTED_BANDS = [
    # "sand_0-5cm_mean",
    # "clay_0-5cm_mean",
    # "soc_0-5cm_mean",
    # "phh2o_0-5cm_mean",
    # "cfvo_0-5cm_mean",
    # "cec_0-5cm_mean",
]

# SoilGrids variables and their conversion factors. Mapped
# integer values are divided by the factor to recover
# conventional units (see notes in header).
VARIABLES = [
    {"name": "bdod", "factor": 100},
    {"name": "cec", "factor": 10},
    {"name": "cfvo", "factor": 10},
    {"name": "clay", "factor": 10},
    {"name": "nitrogen", "factor": 100},
    {"name": "phh2o", "factor": 10},
    {"name": "sand", "factor": 10},
    {"name": "silt", "factor": 10},
    {"name": "soc", "factor": 10},
    {"name": "ocd", "factor": 10},
    {"name": "ocs", "factor": 10},
]

# Point extraction (section 5). Set EXTRACT_XY_POINTS = False to
# skip the batched XY point-value extraction and its per-batch
# CSV exports (e.g. when you only need the raster output).
# True: extract SoilGrids values to XY points
EXTRACT_XY_POINTS = False

# XY points asset. Must contain a 'batch' property with
# integer values matching the loop range below (N_BATCHES).
XY_POINTS_ASSET = (
    "projects/ee-bgcasey-abmi/assets/non_abmi_sites_xy_batch"
)

# Batched extraction parameters.
EXTRACT_SCALE = BASE_SCALE_M  # 250 m (COARSE_SCALE for 1 km)
TILE_SCALE = 16  # higher -> more tiles, lower per-tile mem
N_BATCHES = 100  # Match the number of batches assigned in R

PRINT_STATS = False  # min/max check (slow for large AOIs)
USE_TEST_AOI = False  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt)
# Block until every export task finishes so its batch
# EECU-seconds land in the compute report. Costs the full
# export runtime (hours for a province-wide run), so keep it
# False for production runs and turn it on when profiling a
# test AOI.
WAIT_FOR_EXPORTS = False

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

# 3. Build SoilGrids image ----
# Load each variable, apply its conversion factor, and
# combine into a single multiband image. The image is NOT
# clipped to AOI here; clipping is applied later only for
# raster aggregation and export. Point extraction operates
# on the unclipped image so out-of-AOI points still get
# values. Native SoilGrids band names already include the
# variable name, so no renaming is needed.


def load_variable(name, factor):
    """Load one SoilGrids *_mean asset and rescale it.

    Native band names are preserved. Image is global
    (unclipped) and reprojected to EPSG:4326 to avoid
    sample loss from the native Mollweide projection.

    Args:
        name (str): Variable short name (e.g. 'clay').
        factor (float): Conversion factor (mapped / factor
            = conventional units).

    Returns:
        ee.Image: Rescaled, reprojected multiband image
        (global extent, EPSG:4326).
    """
    return (
        ee.Image(BASE_PATH + name + "_mean")
        .divide(factor)
        .toFloat()
        .reproject(crs="EPSG:4326", scale=250)
    )


# 3.1 Determine which variables are needed based on
# SELECTED_BANDS. If the filter is None/empty, load
# everything; otherwise load only variables that
# contribute to the requested bands.
use_filter = bool(SELECTED_BANDS)
if use_filter:
    var_set = {bn.split("_")[0] for bn in SELECTED_BANDS}
    needed_vars = [
        v for v in VARIABLES if v["name"] in var_set
    ]
else:
    needed_vars = VARIABLES

# 3.2 Combine needed variables into a single multiband
# image.
soilgrids = load_variable(
    needed_vars[0]["name"], needed_vars[0]["factor"]
)
for v in needed_vars[1:]:
    soilgrids = soilgrids.addBands(
        load_variable(v["name"], v["factor"])
    )

# 3.3 Apply the band filter to trim to exactly the
# requested bands.
if use_filter:
    soilgrids = soilgrids.select(SELECTED_BANDS)

# 4. Check bands (optional) ----
# Print band names and min/max stats over the AOI. Earth
# Engine is lazy, so the profiler needs an evaluated
# computation (getInfo) to measure per-algorithm EECU use.

if PRINT_STATS or COMPUTE_REPORT:
    with report.section("SoilGrids band names"):
        print(
            "SoilGrids bands:",
            soilgrids.bandNames().getInfo(),
        )

    # 4.1 Print min/max for a subset of bands over the AOI.
    with report.section("Sample band min/max"):
        sample_bands = [
            "clay_0-5cm_mean",
            "sand_0-5cm_mean",
            "soc_0-5cm_mean",
            "phh2o_0-5cm_mean",
        ]
        # Only reduce bands that survived the filter.
        available = soilgrids.bandNames().getInfo()
        for band in sample_bands:
            if band not in available:
                continue
            stats = (
                soilgrids.select(band)
                .reduceRegion(
                    reducer=ee.Reducer.minMax(),
                    geometry=aoi,
                    scale=1000,
                    maxPixels=1e13,
                    bestEffort=True,
                    tileScale=4,
                )
                .getInfo()
            )
            print(band + " Min and Max:", stats)

# 5. Extract SoilGrids values to XY points (batched) ----
# Use sampleRegions to extract the pixel value at each XY
# location. With ~13M points a single extraction exceeds
# GEE's per-tile memory cap, so the points asset is
# pre-tagged with a 'batch' column (set in R before
# upload) and this loop launches one export task per
# batch. Each batch exports a CSV named
# 'soilgrids_xy_batchNN'. Merge the CSVs in R afterward.
# The whole section is skipped when EXTRACT_XY_POINTS is False.

if EXTRACT_XY_POINTS:
    # 5.1 Load XY points.
    xy_points = ee.FeatureCollection(XY_POINTS_ASSET)

    # 5.2 Diagnostic: inspect the batch column to confirm type
    # and value range. If distinct batch values print as
    # strings (e.g. '1', '2', ...) instead of numbers, the
    # column is stored as character and the Filter.eq calls
    # below need to pass strings, e.g.
    #   ee.Filter.eq('batch', ee.Number(b).format())
    if PRINT_STATS or COMPUTE_REPORT:
        with report.section("Batch diagnostics"):
            print("Total points:", xy_points.size().getInfo())
            print(
                "First feature properties:",
                xy_points.first().getInfo(),
            )
            print(
                "Distinct batch values:",
                xy_points.aggregate_array("batch")
                .distinct()
                .sort()
                .getInfo(),
            )

    # 5.3 Launch one export task per batch. Loop runs
    # 1..N_BATCHES (inclusive) to match the 1-indexed batch
    # values assigned in R.
    for b in range(1, N_BATCHES + 1):
        batch_pts = xy_points.filter(ee.Filter.eq("batch", b))
        extracted = soilgrids.sampleRegions(
            collection=batch_pts,
            scale=EXTRACT_SCALE,
            tileScale=TILE_SCALE,
            geometries=False,
        )
        # Zero-pad batch number to 2 digits for tidy filenames.
        batch_str = str(b).zfill(2)
        task = ee.batch.Export.table.toDrive(
            collection=extracted,
            description="soilgrids_xy_batch" + batch_str,
            folder=DRIVE_FOLDER,
            fileNamePrefix="soilgrids_xy_batch" + batch_str,
            fileFormat="CSV",
        )
        task.start()
        export_tasks.append(task)
        print(
            "Started export task:",
            task.config["description"],
        )

# 6. Export raster output (Alberta only) ----
# Clip to the AOI, then export at the target chosen by
# EXPORT_TARGET. "native" writes the ~250 m image in EPSG:3400;
# "reference_grid" aggregates to the ABMI 1 km grid (EPSG:3400,
# area mean) so it stacks with the FABDEM terrain layers.
# Native-resolution exports over Alberta are large; monitor the
# Tasks tab and expect substantial processing time.

target_suffix = (
    "abmi1km" if EXPORT_TARGET == "reference_grid" else "native"
)
soilgrids_ab = soilgrids.clip(aoi_compute)

if EXPORT_TARGET == "reference_grid":
    # setDefaultProjection pins the native base so the
    # reduceResolution inside export_to_reference_grid knows
    # the input resolution before aggregating to 1 km.
    export_tasks.append(
        export_to_reference_grid(
            image=soilgrids_ab.setDefaultProjection(
                crs=GRID_CRS, scale=BASE_SCALE_M
            ),
            aoi=aoi,
            description=f"{TASK_PREFIX}_{target_suffix}",
            folder=DRIVE_FOLDER,
            file_name_prefix=f"{FILE_PREFIX}_{target_suffix}",
            agg_max_pixels=AGG_MAX_PIXELS,
            fill_gaps_px=FILL_GAPS_PX,
            wait=False,
        )
    )
else:
    export_tasks.append(
        export_image_to_drive(
            image=soilgrids_ab.clip(aoi),
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

# 7. Compute usage report ----
# Multiple export tasks are launched above, so this does
# not block on any single task; it writes the collected
# section profiles to gee_compute_reports/.
if WAIT_FOR_EXPORTS:
    for task in export_tasks:
        report.log_task(task)
report.write()

# End of script ----
