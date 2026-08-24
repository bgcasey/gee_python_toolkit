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
#     1 km reference grid or at FOCAL_BASE_M resolution,
#     selected by EXPORT_TARGET (exported to Google Drive)
# notes:
#   This script calculates the Terrain Ruggedness Index
#   (TRI) from the FABDEM bare-earth DEM (30 m, forests and
#   buildings removed), following Riley et al. (1999):
#
#     TRI = sqrt( sum( (z_i - z_0)^2 ) )
#
#   where z_0 is the centre cell and z_i are the cells in a
#   surrounding neighborhood (classically the eight cells of
#   a 3x3 window). TRI is the root-summed-squared elevation
#   difference between a cell and its neighbours; it is high
#   in rugged terrain and near zero on smooth surfaces, in
#   the same units as elevation (metres).
#
#   TRI is computed at FOCAL_BASE_M (50 m) rather than the
#   native 30 m so the 1 km aggregation stays under Earth
#   Engine's per-tile reprojection limit (see
#   utils.gee_utils.to_reference_grid). The 3x3 window is
#   therefore ~150 m across; ruggedness is scale-dependent, so
#   this is a coarser measure than a 30 m TRI. The grid /
#   boundary / aggregation / export plumbing lives in
#   utils/gee_utils.py.
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
    to_reference_grid,
)

# 1. Setup ----

# 1.1 User parameters ----
# FOCAL_BASE_M is the resolution TRI is computed at before
# aggregating to 1 km (>= ~50 m for a full-province run).
FOCAL_BASE_M = 50
TRI_WINDOW_RADIUS = 1  # pixels; 1 = classic 3x3 Riley window
# Raster export target. "reference_grid" aggregates TRI (area
# mean) onto the ABMI 1 km grid so it stacks with the other 1 km
# covariates. "native" skips the aggregation and exports TRI at
# FOCAL_BASE_M in the grid CRS, ungridded - useful for
# inspecting the input to the aggregation.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"
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
    "fabdem_tri_alberta",
    enabled=COMPUTE_REPORT,
)

# 2. Define study area ----
# aoi is the export / crop boundary; aoi_compute adds a ring
# (the TRI window reach) so the neighbourhood sums are unbiased
# at the true AOI edge.
aoi, aoi_compute = define_study_area(
    use_test_aoi=USE_TEST_AOI,
    buffer_m=FOCAL_BASE_M * TRI_WINDOW_RADIUS,
)

# 3. TRI calculation ----
# FABDEM at the FOCAL_BASE_M base projection, then TRI from the
# neighbourhood sums of z, z^2, and the valid-pixel count.
elevation = fabdem_elevation(aoi_compute, base_m=FOCAL_BASE_M)

# Unweighted window (normalize=False keeps each weight at 1
# so the sum reducer returns true sums, not means)
kernel = ee.Kernel.square(
    radius=TRI_WINDOW_RADIUS, units="pixels", normalize=False
)

# Neighbourhood sums of z, z^2, and the valid-pixel count
sum_z = elevation.reduceNeighborhood(
    reducer=ee.Reducer.sum(), kernel=kernel
)
sum_z2 = (
    elevation.multiply(elevation)
    .reduceNeighborhood(reducer=ee.Reducer.sum(), kernel=kernel)
)
count = elevation.reduceNeighborhood(
    reducer=ee.Reducer.count(), kernel=kernel
)

# Sum of squared differences from the centre cell, via
# sum((z_i - z_0)^2) = sum(z_i^2) - 2*z_0*sum(z_i) + N*z_0^2.
# max(0) guards against tiny negative values from floating
# point on flat terrain before the square root.
ssd = (
    sum_z2
    .subtract(elevation.multiply(sum_z).multiply(2))
    .add(elevation.multiply(elevation).multiply(count))
)
tri = ssd.max(0).sqrt().rename("tri")

# Build the image EXPORT_TARGET asks for. "reference_grid"
# aggregates TRI to the 1 km grid (area mean, float32, clipped to
# Alberta; TRI is continuous metres, so no rounding). "native"
# keeps TRI at FOCAL_BASE_M and only clips. stats_scale follows
# so the min/max below is read at the exported resolution.
if EXPORT_TARGET == "reference_grid":
    tri_export = to_reference_grid(tri, aoi)
    stats_scale = 1000
elif EXPORT_TARGET == "native":
    tri_export = tri.clip(aoi)
    stats_scale = FOCAL_BASE_M
else:
    raise ValueError(
        "Unknown EXPORT_TARGET: "
        f"{EXPORT_TARGET!r} (use 'native' or 'reference_grid')"
    )

# 3.1 Check min and max values (optional) ----
# Also runs when COMPUTE_REPORT is on: Earth Engine is
# lazy, so the profiler needs an evaluated computation
# (getInfo) to measure per-algorithm EECU usage.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("TRI min/max (reduceRegion)"):
        stats = tri_export.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=aoi,
            scale=stats_scale,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print("TRI min and max values:", stats)

# 4. Export data ----
# tri_export was already built for the chosen EXPORT_TARGET
# above, so this only picks the matching export call. The 1 km
# path passes aggregate=False because tri_export is on the grid
# already. Set wait=True to block; otherwise monitor progress
# at https://code.earthengine.google.com/tasks
if EXPORT_TARGET == "reference_grid":
    task = export_to_reference_grid(
        image=tri_export,
        aoi=aoi,
        description="FABDEM_TRI_Alberta_abmi1km",
        folder=DRIVE_FOLDER,
        file_name_prefix="fabdem_tri_alberta_abmi1km",
        aggregate=False,
        wait=False,
    )
else:
    task = export_image_to_drive(
        image=tri_export,
        description="FABDEM_TRI_Alberta_native",
        region=aoi,
        folder=DRIVE_FOLDER,
        file_name_prefix="fabdem_tri_alberta_native",
        scale=FOCAL_BASE_M,
        crs=GRID_CRS,
        max_pixels=1e13,
        wait=False,
    )

# 5. Compute usage report ----
# This section waits for the export to finish, records its
# total EECU-seconds, and writes the txt report to
# gee_compute_reports/. Note: a full-province export can
# take hours; for a quick profile use the test AOI.

if WAIT_FOR_EXPORTS:
    report.log_task(task)
report.write()

# End of script ----
