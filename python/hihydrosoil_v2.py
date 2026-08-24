# ---
# title:   HiHydroSoil v2.0 Layers Export
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - HiHydroSoil v2.0 ImageCollections
#     (FutureWater / sat-io)
#   - Hydrologic_Soil_Group_250m Image
#     (FutureWater / sat-io)
#   - AB2020 provincial boundary (EE asset)
#   - XY points asset (may include locations outside AB)
# outputs:
#   - Multiband HiHydroSoil raster clipped to Alberta, at
#     native (~250 m) or on the ABMI 1 km reference grid, per
#     EXPORT_TARGET. COMBINE_OUTPUTS puts the continuous and
#     categorical bands in one file (default) or two.
#   - Per-batch CSVs of point-level extracted values
#     (optional, gated by EXTRACT_XY_POINTS).
# notes:
#   HiHydroSoil v2.0 provides global soil hydraulic
#   properties at 250 m, derived from SoilGrids250m v2.0
#   by FutureWater. Most continuous layers are stored as
#   int16 * 10000 and are rescaled to physical units by
#   multiplying by 0.0001. The Soil Texture Class (stc)
#   and Hydrologic Soil Group (HSG) layers are
#   categorical and are exported without rescaling.
#
#   Aggregation to 1 km goes through export_to_reference_grid
#   (utils.gee_utils): mean() for continuous layers, mode() for
#   categorical (STC, HSG) to avoid averaging class codes.
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

# Make utils importable
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gee_config import DRIVE_FOLDER, PROVINCIAL_BOUNDARY_ASSET
from utils.compute_report import ComputeReport
from utils.gee_utils import (
    export_to_reference_grid,
    initialize_ee,
    to_reference_grid,
)

# 1. Setup ----

# 1.1 User parameters ----
NATIVE_SCALE = 250  # Native resolution (m)
COARSE_SCALE = 1000  # Aggregated resolution (m)
CRS = "EPSG:3400"  # AB 10-TM (Forest)

# Raster export target, applied to both stacks. "native" writes
# the ~250 m images in EPSG:3400; "reference_grid" aggregates
# them onto the ABMI 1 km grid.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"

# Compute ring grown around the aoi before the source is
# clipped, sized at 2x the output scale. Every output pixel -
# a 1 km grid cell or a native pixel - is then built from a
# full neighbourhood rather than one truncated at the aoi
# edge; a 1 km cell can touch the aoi at a corner and still
# reach a full diagonal (1414 m) beyond it. The exported
# image is clipped back to the plain aoi, so the ring never
# widens the output.
AGG_BUFFER_M = 2 * (
    COARSE_SCALE
    if EXPORT_TARGET == "reference_grid"
    else NATIVE_SCALE
)
BUFFER_MAX_ERROR_M = 100

# True: continuous and categorical bands share one raster.
# False: one file per group. Each group still aggregates with
# its own reducer either way.
COMBINE_OUTPUTS = True

if EXPORT_TARGET not in ("native", "reference_grid"):
    raise ValueError(
        "Unknown EXPORT_TARGET: "
        f"{EXPORT_TARGET!r} (use 'native' or 'reference_grid')"
    )

# Base path for HiHydroSoil v2.0 assets.
BASE_PATH = "projects/sat-io/open-datasets/HiHydroSoilv2_0/"

# Optional asset filter. Set to None (or empty list) to
# keep all assets. Otherwise provide a list of asset short
# names from CONTINUOUS_COLLECTIONS and/or
# CATEGORICAL_COLLECTIONS (e.g. ['ksat'], ['stc']). Assets
# not listed here are skipped at load time.
SELECTED_ASSETS = None  # None = keep all assets

# Soil depths present in every ImageCollection. Images are
# named '<VAR>_<depth>_M_250m' (e.g. 'Ksat_0-5cm_M_250m').
DEPTHS = [
    "0-5cm",
    "5-15cm",
    "15-30cm",
    "30-60cm",
    "60-100cm",
    "100-200cm",
]

# Optional depth filter: list depth tokens from DEPTHS, e.g.
# ['0-5cm']. None/empty keeps all depths. Use depth tokens, not
# full system:index values - an index like 'Ksat_0-5cm_M_250m'
# exists only in the ksat collection and would empty the rest.
DEPTH_FILTER = None

unknown_depths = [d for d in (DEPTH_FILTER or []) if d not in DEPTHS]
if unknown_depths:
    raise ValueError(
        f"Unknown DEPTH_FILTER value(s): {unknown_depths}. "
        f"Valid depths are {DEPTHS}."
    )

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

# Point extraction (section 5). Set EXTRACT_XY_POINTS = False to
# skip the batched XY point-value extraction and its per-batch
# CSV exports (e.g. when you only need the raster outputs).
EXTRACT_XY_POINTS = False  # True: extract HiHydroSoil values to XY points

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
# Block until every export task finishes so its batch
# EECU-seconds land in the compute report. Costs the full
# export runtime (hours for a province-wide run), so keep it
# False for production runs and turn it on when profiling a
# test AOI.
WAIT_FOR_EXPORTS = True

# Export tasks started below, for the optional per-task EECU
# logging in the compute-report section at the end.
export_tasks = []

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
# otherwise uses the AB2020 provincial boundary asset.

if USE_TEST_AOI:
    # Small aoi for testing purposes
    aoi = ee.Geometry.Polygon([
        [-113.5, 55.5],  # Top-left corner
        [-113.5, 55.0],  # Bottom-left corner
        [-112.8, 55.0],  # Bottom-right corner
        [-112.8, 55.5],  # Top-right corner
    ])
else:
    aoi = ee.FeatureCollection(
        PROVINCIAL_BOUNDARY_ASSET
    ).geometry()

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


def collection_to_image(asset_name, depth_filter):
    """Load a collection and collapse it to a multiband image.

    Not clipped to AOI here; clipping happens later for raster
    exports only, so XY extraction still gets values outside
    Alberta.

    Depths are matched on the '_<depth>_' token in
    system:index, the only naming component shared across
    assets ('ALFA_0-5cm_M_250m', 'Ksat_0-5cm_M_250m', ...).

    Bands are renamed '<asset_name>_<depth>' in lower case
    (e.g. 'alpha_0-5cm', 'wcpf4-2_100-200cm'), dropping the
    index's variable token, its '_M_250m' suffix and the '_b1'
    that toBands() adds for single-band images. Hyphens are
    kept, so the one underscore separates variable from depth.

    Args:
        asset_name (str): Short asset name (e.g. 'ksat').
        depth_filter (list or None): Depth tokens to keep
            (e.g. ['0-5cm']). If None/empty, all depths are
            kept.

    Returns:
        ee.Image: Multiband image (global extent).
    """
    ic = ee.ImageCollection(BASE_PATH + asset_name)
    if depth_filter:
        depth_filters = [
            ee.Filter.stringContains(
                "system:index", "_" + d + "_"
            )
            for d in depth_filter
        ]
        ic = ic.filter(ee.Filter.Or(*depth_filters))
    # toBands() creates bands '<system:index>_<origBand>'.
    img = ic.toBands()

    # Only the case is normalised ('N' -> 'n'); hyphens are
    # meaningful ('wcpf4-2' = pF 4.2, '0-5cm' = 0 to 5 cm).
    prefix = asset_name.lower()

    def rename_band(bn):
        """'ALFA_0-5cm_M_250m_b1' -> 'alpha_0-5cm'."""
        depth = (
            ee.String(bn)
            .replace("_b1$", "")       # toBands single-band tag
            .replace("^[^_]+_", "")    # provider variable token
            .replace("_M_250m$", "")   # constant suffix
            .toLowerCase()
        )
        return ee.String(prefix + "_").cat(depth)

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
    # Confirms the depth tokens embedded in each index; copy
    # the depth part (e.g. '30-60cm') into DEPTH_FILTER to
    # keep only specific depth(s).
    with report.section("Inspect system:index values"):
        assets_to_inspect = (
            continuous_collections + categorical_collections
        )
        for name in assets_to_inspect:
            ic = ee.ImageCollection(BASE_PATH + name)
            ids = ic.aggregate_array("system:index").getInfo()
            print(name + " system:index values:", ids)

    # 4.2 Print min/max for all continuous bands. One
    # reduceRegion over the whole stack, not one per band: many
    # separate getInfo calls are slow and can trip the EECU
    # profiler. minMax returns '<band>_min' / '<band>_max' keys.
    if has_continuous:
        with report.section("Continuous band min/max"):
            bands = hihydro_continuous.bandNames().getInfo()
            stats = hihydro_continuous.reduceRegion(
                reducer=ee.Reducer.minMax(),
                geometry=aoi,
                scale=1000,
                maxPixels=1e13,
                bestEffort=True,
                tileScale=4,
            ).getInfo()
            for band in bands:
                print(
                    band + " Min and Max:",
                    {
                        "min": stats.get(band + "_min"),
                        "max": stats.get(band + "_max"),
                    },
                )

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
# location. With large point sets a single extraction exceeds
# GEE's per-tile memory cap, so the points asset is pre-tagged
# with a 'batch' column (set in R before upload) and this loop
# launches one export task per batch. Each batch exports a
# separate CSV named 'hihydrosoil_xy_batchNN'. Merge the CSVs
# in R afterward. The whole section is skipped when
# EXTRACT_XY_POINTS is False, or when no assets passed the
# filters in section 2.1 / 3.

if EXTRACT_XY_POINTS and hihydro_combined is not None:
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
        export_tasks.append(task)
        print(
            "Started export task:",
            task.config["description"],
        )

# 6. Clip to the AOI (Alberta only) ----
# Each group is clipped once; section 7 exports the result
# directly or feeds it to the 1 km aggregation.
#
# The 1 km path clips to a BUFFERED aoi. Clipping at 250 m masks
# the pixels outside Alberta, and reduceResolution then averages
# only the unmasked ones, so a 1 km cell straddling the boundary
# would summarise just its covered part. Aggregating from the
# buffered image gives every such cell a full 1 km of input.
# to_reference_grid / export_to_reference_grid clip the
# aggregated result back to the unbuffered aoi, so the buffer
# never reaches the output - but boundary cells do now include
# values from outside Alberta, which is the point.
# The native path aggregates nothing, so it clips to aoi directly.

clip_geom = (
    aoi.buffer(AGG_BUFFER_M, BUFFER_MAX_ERROR_M)
    if AGG_BUFFER_M
    else aoi
)

hihydro_continuous_ab = None
hihydro_categorical_ab = None

if has_continuous:
    hihydro_continuous_ab = hihydro_continuous.clip(clip_geom)
if has_categorical:
    hihydro_categorical_ab = hihydro_categorical.clip(clip_geom)

# 7. Export raster outputs ----
# EXPORT_TARGET picks the resolution, COMBINE_OUTPUTS picks one
# file or two. A group is skipped if no assets passed the filter.
#
# Continuous bands aggregate to 1 km by area mean, categorical
# by modal class. reduceResolution takes one reducer per call,
# so the combined path reduces each group separately and stacks
# the results. setDefaultProjection pins the native base so
# reduceResolution knows the input resolution. These are large
# exports; monitor the Tasks tab.

def export_native(image, description, file_name_prefix):
    """Export an image at NATIVE_SCALE in the grid CRS."""
    task = ee.batch.Export.image.toDrive(
        image=image.clip(aoi).toFloat(),
        description=description,
        folder=DRIVE_FOLDER,
        fileNamePrefix=file_name_prefix,
        region=aoi,
        scale=NATIVE_SCALE,
        crs=CRS,
        maxPixels=1e13,
    )
    task.start()
    export_tasks.append(task)
    print("Started export task:", description)
    return task


def on_native_base(image):
    """Pin an image to the native base for reduceResolution."""
    return image.setDefaultProjection(crs=CRS, scale=NATIVE_SCALE)


if COMBINE_OUTPUTS:
    # 7.1 One file with every band.
    if EXPORT_TARGET == "native":
        # No aggregation, so the stacks just concatenate.
        combined_ab = hihydro_continuous_ab
        if combined_ab is None:
            combined_ab = hihydro_categorical_ab
        elif hihydro_categorical_ab is not None:
            combined_ab = combined_ab.addBands(
                hihydro_categorical_ab
            )
        export_native(
            combined_ab,
            "HiHydroSoil_AB_native",
            "hihydrosoil_ab_native",
        )
    else:
        # Reduce each group with its own reducer, then stack;
        # both are on the grid already, hence aggregate=False.
        parts = []
        if has_continuous:
            parts.append(
                to_reference_grid(
                    on_native_base(hihydro_continuous_ab),
                    aoi,
                    ee.Reducer.mean(),
                )
            )
        if has_categorical:
            parts.append(
                to_reference_grid(
                    on_native_base(hihydro_categorical_ab),
                    aoi,
                    ee.Reducer.mode(),
                )
            )
        combined_1km = parts[0]
        for part in parts[1:]:
            combined_1km = combined_1km.addBands(part)
        export_tasks.append(export_to_reference_grid(
            image=combined_1km,
            aoi=aoi,
            description="HiHydroSoil_AB_abmi1km",
            folder=DRIVE_FOLDER,
            file_name_prefix="hihydrosoil_ab_abmi1km",
            aggregate=False,
            wait=False,
        ))

elif EXPORT_TARGET == "native":
    # 7.2 Separate files, native resolution (~250 m).
    if has_continuous:
        export_native(
            hihydro_continuous_ab,
            "HiHydroSoil_Continuous_AB_native",
            "hihydrosoil_continuous_ab_native",
        )
    if has_categorical:
        export_native(
            hihydro_categorical_ab,
            "HiHydroSoil_Categorical_AB_native",
            "hihydrosoil_categorical_ab_native",
        )

else:
    # 7.3 Separate files, 1 km on the ABMI reference grid.
    if has_continuous:
        export_tasks.append(export_to_reference_grid(
            image=on_native_base(hihydro_continuous_ab),
            aoi=aoi,
            description="HiHydroSoil_Continuous_AB_abmi1km",
            folder=DRIVE_FOLDER,
            file_name_prefix="hihydrosoil_continuous_ab_abmi1km",
            aggregate=True,
            reducer=ee.Reducer.mean(),
            wait=False,
        ))
    if has_categorical:
        export_tasks.append(export_to_reference_grid(
            image=on_native_base(hihydro_categorical_ab),
            aoi=aoi,
            description="HiHydroSoil_Categorical_AB_abmi1km",
            folder=DRIVE_FOLDER,
            file_name_prefix="hihydrosoil_categorical_ab_abmi1km",
            aggregate=True,
            reducer=ee.Reducer.mode(),
            wait=False,
        ))

# 8. Compute usage report ----
# Multiple export tasks are launched above, so this does
# not block on any single task; it writes the collected
# section profiles to gee_compute_reports/.
if WAIT_FOR_EXPORTS:
    for task in export_tasks:
        report.log_task(task)
report.write()

# End of script ----
