# ---
# title:   FABDEM TRI (Terrain Ruggedness Index) for Alberta
# author:  Brendan Casey
# created: 2026-07-11
# inputs:
#   - FABDEM ImageCollection
#     (projects/sat-io/open-datasets/FABDEM)
#   - FAO GAUL province boundaries
# outputs:
#   - TRI GeoTIFF for Alberta (exported to Google Drive)
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

from _gee_config import DRIVE_FOLDER
from utils.compute_report import ComputeReport
from utils.gee_utils import export_image_to_drive, initialize_ee

# 1. Setup ----

# 1.1 User parameters ----
EXPORT_SCALE = 30  # meters
EXPORT_CRS = "EPSG:4326"
TRI_WINDOW_RADIUS = 1  # pixels; 1 = classic 3x3 Riley window
PRINT_STATS = True  # min/max check (slow for large AOIs)
USE_TEST_AOI = True  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt);
# blocks until the export task finishes

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
# This section defines the export geometry. It uses a
# small test polygon when USE_TEST_AOI is True; otherwise
# it filters the FAO GAUL provinces for Alberta.

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

# 3. TRI calculation ----
# This section mosaics the FABDEM collection, pins it to a
# metric projection (EPSG:3402, Alberta 10-TM) so the
# neighbourhood sits on a consistent 30 m grid, clips to the
# AOI, and derives TRI. It produces a single-band TRI image.

elevation = (
    ee.ImageCollection("projects/sat-io/open-datasets/FABDEM")
    .mosaic()
    .setDefaultProjection("EPSG:3402", None, 30)
    .clip(aoi)
    .double()
)

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

# 3.1 Check min and max values (optional) ----
# Also runs when COMPUTE_REPORT is on: Earth Engine is
# lazy, so the profiler needs an evaluated computation
# (getInfo) to measure per-algorithm EECU usage.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("TRI min/max (reduceRegion)"):
        stats = tri.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=aoi,
            scale=EXPORT_SCALE,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print("TRI min and max values:", stats)

# 4. Export data ----
# This section exports the TRI image to Google Drive as a
# GeoTIFF. Set wait=True to block until the task finishes;
# otherwise monitor progress at
# https://code.earthengine.google.com/tasks

task = export_image_to_drive(
    image=tri,
    description="FABDEM_TRI_Alberta",
    region=aoi,
    folder=DRIVE_FOLDER,
    file_name_prefix="fabdem_tri_alberta",
    scale=EXPORT_SCALE,
    crs=EXPORT_CRS,
    max_pixels=1e13,
    wait=False,
)

# 5. Compute usage report ----
# This section waits for the export to finish, records its
# total EECU-seconds, and writes the txt report to
# gee_compute_reports/. Note: a full-province export can
# take hours; for a quick profile use the test AOI.

report.log_task(task)
report.write()

# End of script ----
