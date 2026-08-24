# ---
# title:   Hydrologically Adjusted Elevations (HAND)
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - MERIT Hydro (MERIT/Hydro/v1_0_1)
#   - AB2020 provincial boundary (EE asset)
# outputs:
#   - HAND GeoTIFF for Alberta, at native (~90 m) or on the
#     ABMI 1 km reference grid, per EXPORT_TARGET (exported to
#     Google Drive).
# notes:
#   This script extracts the hydrologically adjusted
#   elevations (Height Above Nearest Drainage - HAND) from
#   the MERIT Hydro dataset, clips to the AOI, and exports
#   the HAND band to Google Drive. Map visualization layers
#   from the original GEE JavaScript are dropped.
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

from _gee_config import DRIVE_FOLDER, PROVINCIAL_BOUNDARY_ASSET
from utils.compute_report import ComputeReport
from utils.gee_utils import (
    export_image_to_drive,
    export_to_reference_grid,
    initialize_ee,
)

# 1. Setup ----

# 1.1 User parameters ----
EXPORT_SCALE = 92.77  # meters (MERIT Hydro native ~90 m)
EXPORT_CRS = "EPSG:3400"  # AB 10-TM (Forest)
# Raster export target. "native" exports the ~90 m HAND image;
# "reference_grid" aggregates (area mean) onto the ABMI 1 km
# grid so it stacks with the other 1 km covariates. HAND is a
# stored MERIT Hydro band (pyramided), so 90 m -> 1 km is well
# under Earth Engine's per-tile reprojection limit.
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
    else EXPORT_SCALE
)
BUFFER_MAX_ERROR_M = 100
PRINT_STATS = True  # min/max check (slow for large AOIs)
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
# Profiles EECU usage per section and per export task.
# Best used with USE_TEST_AOI = True to find choke
# points cheaply before a full-province run.
report = ComputeReport(
    "hydrologically_adjusted_elevation",
    enabled=COMPUTE_REPORT,
)

# 2. Define study area ----
# This section defines the export geometry. It uses a
# small test polygon when USE_TEST_AOI is True; otherwise
# it uses the AB2020 provincial boundary asset.

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

# 3. Extract HAND band ----
# This section clips the MERIT Hydro image to the AOI and
# selects the 'hnd' (Height Above Nearest Drainage) band,
# renaming it to 'HAND'. It produces a single-band image.

merit_hydro = ee.Image("MERIT/Hydro/v1_0_1")
# Aggregation reads from the ring (AGG_BUFFER_M); the 1 km
# result is clipped back to the plain aoi downstream.
clip_geom = (
    aoi.buffer(AGG_BUFFER_M, BUFFER_MAX_ERROR_M)
    if AGG_BUFFER_M
    else aoi
)

hand = merit_hydro.clip(clip_geom).select("hnd").rename("HAND")

# 3.1 Check min and max values (optional) ----
# Also runs when COMPUTE_REPORT is on: Earth Engine is
# lazy, so the profiler needs an evaluated computation
# (getInfo) to measure per-algorithm EECU usage.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("HAND min/max (reduceRegion)"):
        stats = hand.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=aoi,
            scale=EXPORT_SCALE,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print("HAND min and max values:", stats)

# 4. Export data ----
# Export at the target chosen by EXPORT_TARGET. "native" writes
# the ~90 m HAND image; "reference_grid" aggregates it (area
# mean) onto the ABMI 1 km grid. setDefaultProjection pins the
# native base so reduceResolution knows the input resolution.
# Set wait=True to block; otherwise monitor progress at
# https://code.earthengine.google.com/tasks

if EXPORT_TARGET == "reference_grid":
    task = export_to_reference_grid(
        image=hand.setDefaultProjection(
            crs=EXPORT_CRS, scale=EXPORT_SCALE
        ),
        aoi=aoi,
        description="HAND_AB_abmi1km",
        folder=DRIVE_FOLDER,
        file_name_prefix="hydrologically_adjusted_elevations_abmi1km",
        aggregate=True,
        wait=False,
    )
elif EXPORT_TARGET == "native":
    task = export_image_to_drive(
        image=hand.clip(aoi),
        description="HAND_AB_native",
        region=aoi,
        folder=DRIVE_FOLDER,
        file_name_prefix="hydrologically_adjusted_elevations_native",
        scale=EXPORT_SCALE,
        crs=EXPORT_CRS,
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
