# ---
# title:   SoilGrids 250m v2.0 Layers Export
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - ISRIC SoilGrids 250m v2.0 mean Images
#     (projects/soilgrids-isric/*_mean)
#   - Alberta boundary (FAO GAUL level1)
#   - XY points (may include locations outside Alberta)
# outputs:
#   - Multiband SoilGrids image clipped to Alberta,
#     exported at native (~250 m) and 1000 m resolution.
#   - Per-batch CSVs of point-level extracted soil values
#     for ALL points (including those outside the AOI).
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
#   1000 m aggregation uses mean reduction on continuous
#   layers. setDefaultProjection is required before
#   reduceResolution when aggregating by a factor > 64.
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

import os
import sys

import ee

# Make utils importable regardless of the working
# directory VS Code runs the script from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gee_config import DRIVE_FOLDER
from utils.compute_report import ComputeReport
from utils.gee_utils import initialize_ee

# 1. Setup ----

# 1.1 User parameters ----
NATIVE_SCALE = 250  # Native resolution (m)
COARSE_SCALE = 1000  # Aggregated resolution (m)
CRS = "EPSG:4326"

# Base path for SoilGrids 250m v2.0 assets.
BASE_PATH = "projects/soilgrids-isric/"

# Optional band filter. Set to None (or empty list) to
# keep all bands. Otherwise provide a list of band names
# in <variable>_<depth>_mean format. Any variable not
# represented is skipped at load time; remaining variables
# are loaded fully and filtered after.
SELECTED_BANDS = [
    "sand_0-5cm_mean",
    "clay_0-5cm_mean",
    "soc_0-5cm_mean",
    "phh2o_0-5cm_mean",
    "cfvo_0-5cm_mean",
    "cec_0-5cm_mean",
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

# XY points asset. Must contain a 'batch' property with
# integer values matching the loop range below (N_BATCHES).
XY_POINTS_ASSET = (
    "projects/ee-bgcasey-abmi/assets/non_abmi_sites_xy_batch"
)

# Batched extraction parameters.
EXTRACT_SCALE = NATIVE_SCALE  # 250 m (COARSE_SCALE for 1 km)
TILE_SCALE = 16  # higher -> more tiles, lower per-tile mem
N_BATCHES = 100  # Match the number of batches assigned in R

PRINT_STATS = True  # min/max check (slow for large AOIs)
USE_TEST_AOI = True  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt)

# 1.2 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.3 Set up compute usage report ----
# Profiles EECU usage per section. Best used with
# USE_TEST_AOI = True to find choke points cheaply.
report = ComputeReport(
    "soil_grids_250",
    out_dir=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "gee_compute_reports",
    ),
    enabled=COMPUTE_REPORT,
)

# 2. Define study area ----
# Uses a small test polygon when USE_TEST_AOI is True;
# otherwise filters the FAO GAUL provinces for Alberta.

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
    print(
        "Started export task:",
        task.config["description"],
    )

# 6. Aggregate to 1000 m (Alberta only) ----
# Clip to AOI first so aggregation only operates on Alberta
# pixels. setDefaultProjection is required before
# reduceResolution when aggregating by more than a factor
# of 64.

soilgrids_ab = soilgrids.clip(aoi)

soilgrids_1km = (
    soilgrids_ab.setDefaultProjection(
        crs=CRS, scale=NATIVE_SCALE
    )
    .reduceResolution(
        reducer=ee.Reducer.mean(), maxPixels=1024
    )
    .reproject(crs=CRS, scale=COARSE_SCALE)
    .toFloat()
)

# 7. Export raster outputs ----
# Export native-resolution and 1000 m images (Alberta
# only) to Google Drive. Native-resolution exports over
# Alberta are large; monitor the Tasks tab and expect
# substantial processing time.

# 7.1 Native resolution (~250 m), Alberta only.
ee.batch.Export.image.toDrive(
    image=soilgrids_ab,
    description="SoilGrids_AB_250m",
    folder=DRIVE_FOLDER,
    fileNamePrefix="soilgrids_ab_250m",
    region=aoi,
    scale=NATIVE_SCALE,
    crs=CRS,
    maxPixels=1e13,
).start()

# 7.2 1000 m aggregated, Alberta only.
ee.batch.Export.image.toDrive(
    image=soilgrids_1km,
    description="SoilGrids_AB_1000m",
    folder=DRIVE_FOLDER,
    fileNamePrefix="soilgrids_ab_1000m",
    region=aoi,
    scale=COARSE_SCALE,
    crs=CRS,
    maxPixels=1e13,
).start()

# 8. Compute usage report ----
# Multiple export tasks are launched above, so this does
# not block on any single task; it writes the collected
# section profiles to gee_compute_reports/.
report.write()

# End of script ----
