# ---
# title:   HiHydroSoil v2.0 Layers Export
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - HiHydroSoil v2.0 ImageCollections
#     (FutureWater / sat-io)
#   - Hydrologic_Soil_Group_250m Image
#     (FutureWater / sat-io)
#   - Alberta boundary (FAO GAUL level1)
#   - XY points asset (may include locations outside AB)
# outputs:
#   - Multiband HiHydroSoil images clipped to Alberta,
#     exported at native (~250 m) and 1000 m resolution.
#   - Per-batch CSVs of point-level extracted values.
# notes:
#   HiHydroSoil v2.0 provides global soil hydraulic
#   properties at 250 m, derived from SoilGrids250m v2.0
#   by FutureWater. Most continuous layers are stored as
#   int16 * 10000 and are rescaled to physical units by
#   multiplying by 0.0001. The Soil Texture Class (stc)
#   and Hydrologic Soil Group (HSG) layers are
#   categorical and are exported without rescaling.
#
#   Most assets are ImageCollections representing the six
#   standard soil depths. They are collapsed to a
#   multiband image using .toBands(), producing band
#   names of the form <index>_<asset>. The
#   Hydrologic_Soil_Group_250m asset is a single Image.
#
#   1000 m exports use ee.Reducer.mean() for continuous
#   layers and ee.Reducer.mode() for categorical layers
#   (STC, HSG) to avoid producing meaningless averages of
#   class codes.
#
#   Citation:
#   Simons, G.W.H., R. Koster, P. Droogers. 2020.
#   HiHydroSoil v2.0 - A high resolution soil map of
#   global hydraulic properties. FutureWater Report 213.
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
CRS = "EPSG:4326"  # Alternative: 'EPSG:3400' (AB 10-TM)

# Base path for HiHydroSoil v2.0 assets.
BASE_PATH = "projects/sat-io/open-datasets/HiHydroSoilv2_0/"

# Optional asset filter. Set to None (or empty list) to
# keep all assets. Otherwise provide a list of asset short
# names from CONTINUOUS_COLLECTIONS and/or
# CATEGORICAL_COLLECTIONS (e.g. ['ksat'], ['stc']). Assets
# not listed here are skipped at load time.
SELECTED_ASSETS = ["ksat"]

# Optional depth filter (per-image, ImageCollection assets
# only). HiHydroSoil collections contain one image per
# soil depth (plus aggregated topsoil/subsoil layers). The
# exact system:index per image is provider-specific. The
# stats section prints available system:index values for
# every selected asset so you can copy the correct strings
# here on a follow-up run. Set to None to keep all images.
DEPTH_FILTER = ["Ksat_0-5cm_M_250m"]

# Continuous (float) ImageCollection assets. Rescaled by
# multiplying with 0.0001.
CONTINUOUS_COLLECTIONS = [
    "alpha",       # Mualem-van Genuchten alpha (1/cm)
    "crit-wilt",   # Water content pF3 - pF4.2 (m3/m3)
    "field-crit",  # Water content pF2 - pF3 (m3/m3)
    "ksat",        # Saturated hydraulic conductivity (cm/d)
    "N",           # Mualem-van Genuchten N (-)
    "ormc",        # Organic matter content (%)
    "sat-field",   # Water content sat - pF2 (m3/m3)
    "wcavail",     # Available water content (m3/m3)
    "wcpf2",       # Water content at pF2 (m3/m3)
    "wcpf3",       # Water content at pF3 (m3/m3)
    "wcpf4-2",     # Water content at pF4.2 (m3/m3)
    "wcres",       # Residual water content (m3/m3)
    "wcsat",       # Saturated water content (m3/m3)
]

# Categorical ImageCollection assets. NOT rescaled; use
# mode() for aggregation.
CATEGORICAL_COLLECTIONS = [
    "stc",  # Soil Texture Class (1-6)
]

# XY points asset. Must contain a 'batch' property with
# integer values matching the loop range below (N_BATCHES).
XY_POINTS_ASSET = (
    "projects/ee-bgcasey-abmi/assets/non_abmi_sites_xy_batch"
)

# Batched extraction parameters.
EXTRACT_SCALE = NATIVE_SCALE  # 250 m (COARSE_SCALE for 1 km)
TILE_SCALE = 16  # higher -> more tiles, lower per-tile mem
N_BATCHES = 50  # Match the number of batches assigned in R

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
    "hihydrosoil_v2",
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

# 2.1 Apply the asset filter. Empty/None = no filter. The
# Hydrologic_Soil_Group asset is loaded separately (single
# Image). Use the name 'hydrologic_soil_group' to include
# it in the filter.
use_filter = bool(SELECTED_ASSETS)
include_hsg = True
if use_filter:
    continuous_collections = [
        name for name in CONTINUOUS_COLLECTIONS
        if name in SELECTED_ASSETS
    ]
    categorical_collections = [
        name for name in CATEGORICAL_COLLECTIONS
        if name in SELECTED_ASSETS
    ]
    include_hsg = "hydrologic_soil_group" in SELECTED_ASSETS
else:
    continuous_collections = list(CONTINUOUS_COLLECTIONS)
    categorical_collections = list(CATEGORICAL_COLLECTIONS)

# 3. Build HiHydroSoil image ----
# Process each collection into a multiband image, rescale
# continuous layers, then combine all layers into a single
# multiband image.


def collection_to_image(asset_name, index_filter):
    """Load a collection and collapse it to a multiband image.

    Optionally filters by system:index, collapses to a
    multiband image via toBands(), and assigns clean band
    names. The image is NOT clipped to AOI here; clipping
    is applied later only for raster exports so XY
    extraction still gets values outside Alberta.

    Band naming rules:
      - With a filter, each band is renamed to the
        system:index of its source image.
      - Without a filter, bands are renamed to
        '<asset_name>_<system:index>' to disambiguate
        across multiple assets.
    The trailing '_b1' that toBands() inserts for
    single-band images is stripped in both cases.

    Args:
        asset_name (str): Short asset name (e.g. 'ksat').
        index_filter (list or None): system:index strings
            to keep. If None, all images are kept.

    Returns:
        ee.Image: Multiband image (global extent).
    """
    ic = ee.ImageCollection(BASE_PATH + asset_name)
    has_filter = bool(index_filter)
    if has_filter:
        ic = ic.filter(
            ee.Filter.inList("system:index", index_filter)
        )
    # toBands() creates bands '<system:index>_<origBand>'.
    img = ic.toBands()

    # Strip trailing '_b1' from single-band source images,
    # then optionally prefix with the asset name.
    def rename_band(bn):
        stripped = ee.String(bn).replace("_b1$", "")
        if has_filter:
            return stripped  # system:index already clean
        return ee.String(asset_name).cat("_").cat(stripped)

    band_names = img.bandNames().map(rename_band)
    return img.rename(band_names)


# 3.1 Continuous collections (rescale by 0.0001).
continuous_images = [
    collection_to_image(name, DEPTH_FILTER)
    .multiply(0.0001)
    .toFloat()
    for name in continuous_collections
]

# 3.2 Categorical collections (no rescale, Int16).
categorical_images = [
    collection_to_image(name, DEPTH_FILTER).toInt16()
    for name in categorical_collections
]

# 3.3 Hydrologic Soil Group (single Image, categorical),
# only if it passed the asset filter. Not clipped here.
hsg = None
if include_hsg:
    hsg = (
        ee.Image(BASE_PATH + "Hydrologic_Soil_Group_250m")
        .rename("hydrologic_soil_group")
        .toInt16()
    )

# 3.4 Combine into continuous and categorical multiband
# images. Either group may be empty after filtering;
# downstream sections guard with the has_* flags.
has_continuous = len(continuous_images) > 0
has_categorical = (
    len(categorical_images) > 0 or hsg is not None
)

hihydro_continuous = None
if has_continuous:
    hihydro_continuous = ee.Image(continuous_images[0])
    for img in continuous_images[1:]:
        hihydro_continuous = hihydro_continuous.addBands(img)

hihydro_categorical = None
if has_categorical:
    if len(categorical_images) > 0:
        hihydro_categorical = ee.Image(categorical_images[0])
        for img in categorical_images[1:]:
            hihydro_categorical = (
                hihydro_categorical.addBands(img)
            )
        if hsg is not None:
            hihydro_categorical = (
                hihydro_categorical.addBands(hsg)
            )
    elif hsg is not None:
        hihydro_categorical = hsg

# 3.5 Combine continuous and categorical stacks into a
# single extraction image. At bufferSize = 0 there is no
# reducer distinction, so sampleRegions just reads the
# pixel value for both. Either group may be empty.
hihydro_combined = None
if has_continuous and has_categorical:
    hihydro_combined = hihydro_continuous.addBands(
        hihydro_categorical
    )
elif has_continuous:
    hihydro_combined = hihydro_continuous
elif has_categorical:
    hihydro_combined = hihydro_categorical

# 4. Check bands (optional) ----
# Print band names, available system:index values, and
# min/max stats. Earth Engine is lazy, so the profiler
# needs an evaluated computation (getInfo) to measure
# per-algorithm EECU usage.

if PRINT_STATS or COMPUTE_REPORT:
    with report.section("HiHydroSoil band names"):
        if has_continuous:
            print(
                "HiHydroSoil Continuous bands:",
                hihydro_continuous.bandNames().getInfo(),
            )
        if has_categorical:
            print(
                "HiHydroSoil Categorical bands:",
                hihydro_categorical.bandNames().getInfo(),
            )

    # 4.1 Inspect collection contents (system:index).
    # Use this output to populate DEPTH_FILTER if you want
    # to keep only specific depth(s).
    with report.section("Inspect system:index values"):
        assets_to_inspect = (
            continuous_collections + categorical_collections
        )
        for name in assets_to_inspect:
            ic = ee.ImageCollection(BASE_PATH + name)
            ids = ic.aggregate_array("system:index").getInfo()
            print(name + " system:index values:", ids)

    # 4.2 Print min/max for all continuous bands.
    if has_continuous:
        with report.section("Continuous band min/max"):
            bands = hihydro_continuous.bandNames().getInfo()
            for band in bands:
                stats = (
                    hihydro_continuous.select(band)
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

    # 4.3 Confirm the native projection of a sample asset.
    with report.section("Projection check"):
        hh = ee.Image(
            "projects/sat-io/open-datasets/HiHydroSoilv2_0/"
            "ksat/Ksat_0-5cm_M_250m"
        )
        print(
            "HiHydroSoil projection:",
            hh.projection().getInfo(),
        )

# 5. Extract HiHydroSoil values to XY points (batched) ----
# Use sampleRegions to extract the pixel value at each XY
# location. With large point sets a single extraction
# exceeds GEE's per-tile memory cap, so the points asset is
# pre-tagged with a 'batch' column (set in R before upload)
# and this loop launches one export task per batch. Each
# batch exports a separate CSV named
# 'hihydrosoil_xy_batchNN'. Merge the CSVs in R afterward.

# 5.1 Load XY points.
xy_points = ee.FeatureCollection(XY_POINTS_ASSET)

# 5.2 Diagnostic: inspect the batch column. If distinct
# batch values print as strings (e.g. '1', '2', ...)
# instead of numbers, the column is stored as character
# and the Filter.eq calls below need to pass strings, e.g.
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
        # First rows of batch 1 for a sanity check.
        if hihydro_combined is not None:
            batch1 = xy_points.filter(
                ee.Filter.eq("batch", 1)
            )
            sample = hihydro_combined.sampleRegions(
                collection=batch1.limit(5),
                scale=EXTRACT_SCALE,
                tileScale=TILE_SCALE,
                geometries=False,
            )
            print(
                "Batch 1 sample extraction (first 5):",
                sample.getInfo(),
            )

# 5.3 Launch one export task per batch. Loop runs
# 1..N_BATCHES (inclusive) to match the 1-indexed batch
# values assigned in R.
if hihydro_combined is not None:
    for b in range(1, N_BATCHES + 1):
        batch_pts = xy_points.filter(
            ee.Filter.eq("batch", b)
        )
        extracted = hihydro_combined.sampleRegions(
            collection=batch_pts,
            scale=EXTRACT_SCALE,
            tileScale=TILE_SCALE,
            geometries=False,
        )
        # Zero-pad batch number to 2 digits for filenames.
        batch_str = str(b).zfill(2)
        task = ee.batch.Export.table.toDrive(
            collection=extracted,
            description="hihydrosoil_xy_batch" + batch_str,
            folder=DRIVE_FOLDER,
            fileNamePrefix="hihydrosoil_xy_batch" + batch_str,
            fileFormat="CSV",
        )
        task.start()
        print("Started export task:", task.config[
            "description"
        ])

# 6. Aggregate to 1000 m (Alberta only) ----
# Clip to AOI first so aggregation and export only operate
# on Alberta pixels. setDefaultProjection is required
# before reduceResolution when aggregating by more than a
# factor of 64.

hihydro_continuous_ab = None
hihydro_categorical_ab = None
hihydro_continuous_1km = None
hihydro_categorical_1km = None

# 6.1 Continuous: clip to AB and aggregate to 1 km (mean).
if has_continuous:
    hihydro_continuous_ab = hihydro_continuous.clip(aoi)
    hihydro_continuous_1km = (
        hihydro_continuous_ab.setDefaultProjection(
            crs=CRS, scale=NATIVE_SCALE
        )
        .reduceResolution(
            reducer=ee.Reducer.mean(), maxPixels=1024
        )
        .reproject(crs=CRS, scale=COARSE_SCALE)
        .toFloat()
    )

# 6.2 Categorical: clip to AB and aggregate to 1 km (mode).
if has_categorical:
    hihydro_categorical_ab = hihydro_categorical.clip(aoi)
    hihydro_categorical_1km = (
        hihydro_categorical_ab.setDefaultProjection(
            crs=CRS, scale=NATIVE_SCALE
        )
        .reduceResolution(
            reducer=ee.Reducer.mode(), maxPixels=1024
        )
        .reproject(crs=CRS, scale=COARSE_SCALE)
        .toInt16()
    )

# 7. Export raster outputs ----
# Export native-resolution and 1000 m images to Google
# Drive. These are large exports (especially native
# continuous). Monitor the Tasks tab and expect
# substantial processing time. Either group is skipped if
# no assets passed the filter.

if has_continuous:
    # 7.1 Continuous layers - native resolution (~250 m).
    ee.batch.Export.image.toDrive(
        image=hihydro_continuous_ab,
        description="HiHydroSoil_Continuous_AB_250m",
        folder=DRIVE_FOLDER,
        fileNamePrefix="hihydrosoil_continuous_ab_250m",
        region=aoi,
        scale=NATIVE_SCALE,
        crs=CRS,
        maxPixels=1e13,
    ).start()

    # 7.2 Continuous layers - 1000 m.
    ee.batch.Export.image.toDrive(
        image=hihydro_continuous_1km,
        description="HiHydroSoil_Continuous_AB_1000m",
        folder=DRIVE_FOLDER,
        fileNamePrefix="hihydrosoil_continuous_ab_1000m",
        region=aoi,
        scale=COARSE_SCALE,
        crs=CRS,
        maxPixels=1e13,
    ).start()

if has_categorical:
    # 7.3 Categorical layers - native resolution (~250 m).
    ee.batch.Export.image.toDrive(
        image=hihydro_categorical_ab,
        description="HiHydroSoil_Categorical_AB_250m",
        folder=DRIVE_FOLDER,
        fileNamePrefix="hihydrosoil_categorical_ab_250m",
        region=aoi,
        scale=NATIVE_SCALE,
        crs=CRS,
        maxPixels=1e13,
    ).start()

    # 7.4 Categorical layers - 1000 m.
    ee.batch.Export.image.toDrive(
        image=hihydro_categorical_1km,
        description="HiHydroSoil_Categorical_AB_1000m",
        folder=DRIVE_FOLDER,
        fileNamePrefix="hihydrosoil_categorical_ab_1000m",
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
