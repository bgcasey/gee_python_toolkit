# ---
# title:   FABDEM TPI for Alberta
# author:  Brendan Casey
# created: 2026-07-11
# inputs:
#   - FABDEM ImageCollection
#     (projects/sat-io/open-datasets/FABDEM)
#   - AB2020 provincial boundary (Earth Engine asset;
#     PROVINCIAL_BOUNDARY_ASSET) for the full-province crop
# outputs:
#   - One 30 m TPI GeoTIFF per focal radius for Alberta,
#     cropped to the province and exported to Google Drive
#     (large; Earth Engine shards it into multiple tiles).
#     Aggregate to the ABMI 1 km reference grid downstream
#     with r/aggregate_tpi_to_grid.R.
# notes:
#   This script calculates the Topographic Position Index
#   (TPI) from the FABDEM bare-earth DEM (30 m, forests and
#   buildings removed): elevation minus the mean elevation of
#   a surrounding neighborhood. The collection is mosaicked,
#   the AOI is defined, and TPI is computed at each focal
#   radius in TPI_RADII at the DEM's native 30 m resolution.
#
#   The TPI surface is exported at native 30 m (EPSG:3400,
#   anchored to the ABMI grid origin) and is NOT aggregated
#   here. Aggregating a computed 30 m layer straight to 1 km
#   over the whole province exceeds Earth Engine's per-tile
#   reprojection limit ("Reprojection output too large"):
#   filling one 1 km output tile forces the 30 m focalMean
#   over a ~256 km footprint (~8600x8600 px). So the 30 m TPI
#   is exported as-is and the 30 m -> 1 km area mean, grid
#   snap, and province crop are done downstream in R
#   (r/aggregate_tpi_to_grid.R), where there is no such limit.
#
#   Reference grid the R step targets:
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
# Export grid: the 30 m TPI is exported in the ABMI grid's CRS
# and anchored to its lattice origin, so the downstream 1 km
# aggregation lands cleanly on the reference grid.
# EXPORT_CRS_TRANSFORM is [xScale, xShear, originX, yShear,
# yScale, originY].
GRID_CRS = "EPSG:3400"  # NAD83 / Alberta 10-TM (Forest)
EXPORT_SCALE = 30  # metres (native FABDEM resolution)
# TPI is exported at 30 m anchored to the ABMI 1 km grid's
# top-left lattice node (west 170616.1822, north 6659532.4311)
# so pixels sit on a fixed, reproducible 30 m lattice. Exact
# nesting into 1 km cells is not required: the downstream R
# step aggregates by area-weighted mean, which handles the
# non-integer 1000/30 ratio. yScale is negative (rows run
# north -> south).
EXPORT_CRS_TRANSFORM = [
    EXPORT_SCALE, 0, 170616.1822,
    0, -EXPORT_SCALE, 6659532.4311,
]
# Full-province runs are cropped to this boundary (an Earth
# Engine table asset). Upload AB2020_provincial_boundary.shp
# to your EE assets and set the ID here. Only used when
# USE_TEST_AOI is False.
PROVINCIAL_BOUNDARY_ASSET = (
    "projects/ee-bgcasey-abmi/assets/AB2020_provincial_boundary"
)
TPI_RADII = [1000]  # one export per radius
TPI_WINDOW_SHAPE = "circle"  # "circle" or "square"
TPI_UNITS = "meters"  # "meters" or "pixels"
# Compute buffer (metres). The DEM is clipped to the AOI grown
# by this ring before the focal mean, so cells at the true AOI
# edge see a full neighborhood of real elevation instead of
# masked pixels (which would bias TPI inward). The ring is
# discarded when the output is clipped back to the AOI, and
# compute stays bounded to AOI + ring rather than the whole
# DEM. It must be >= the largest focal reach; native pixels
# are 30 m, so "pixels" radii are converted to metres.
COMPUTE_BUFFER_M = max(TPI_RADII) * (30 if TPI_UNITS == "pixels" else 1)
USE_TEST_AOI = False  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt);
# blocks until the export task finishes

# 1.2 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.3 Set up compute usage report ----
# Records total EECU-seconds for each export task.
# Best used with USE_TEST_AOI = True to gauge compute
# cost cheaply before a full-province run.
report = ComputeReport(
    "fabdem_tpi_alberta",
    enabled=COMPUTE_REPORT,
)

# 2. Define study area ----
# This section defines two geometries: aoi, the export/crop
# boundary (a small test polygon when USE_TEST_AOI is True,
# otherwise the AB2020 provincial boundary asset so the
# full-province export is cropped to Alberta's exact outline);
# and aoi_compute, aoi grown by COMPUTE_BUFFER_M so the DEM
# carries real elevation just outside the true edge and the
# focal mean there is unbiased.

if USE_TEST_AOI:
    # Small aoi for testing purposes
    aoi = ee.Geometry.Polygon([
        [-113.5, 55.5],  # Top-left corner
        [-113.5, 55.0],  # Bottom-left corner
        [-112.8, 55.0],  # Bottom-right corner
        [-112.8, 55.5],  # Top-right corner
    ])
else:
    aoi = ee.FeatureCollection(PROVINCIAL_BOUNDARY_ASSET).geometry()

# Buffered AOI for computation only. maxError (100 m) keeps
# buffering the detailed boundary cheap; the ring never
# reaches the output, which is clipped back to aoi.
aoi_compute = aoi.buffer(COMPUTE_BUFFER_M, 100)

# 3. Prepare the DEM ----
# This section mosaics the FABDEM collection, pins it to the
# reference grid's metric projection (EPSG:3400, Alberta
# 10-TM) at the native 30 m resolution so the focal radius
# maps to real ground distance and the later 30 m -> 1 km
# aggregation stays within one CRS. It clips to the buffered
# AOI (aoi_compute) so the focal mean has real elevation on
# all sides of the true AOI edge. The same elevation image
# feeds every focal radius below.

elevation = (
    ee.ImageCollection("projects/sat-io/open-datasets/FABDEM")
    .mosaic()
    .setDefaultProjection(GRID_CRS, None, 30)
    .clip(aoi_compute)
    .double()
)

# 4. Compute and export 30 m TPI per focal radius ----
# For each radius in TPI_RADII, TPI is elevation minus the
# neighborhood mean elevation at the native 30 m resolution,
# cropped to the province and exported to Google Drive as its
# own (tiled) GeoTIFF. It is NOT aggregated to 1 km here (see
# the header note); r/aggregate_tpi_to_grid.R does that.
# Larger radii use bigger focal kernels and cost
# proportionally more compute; the per-task batch
# EECU-seconds in the report show where.
# Set wait=True on the export to block; otherwise monitor
# progress at https://code.earthengine.google.com/tasks

tasks = []
for radius in TPI_RADII:
    # TPI at native 30 m: elevation minus neighborhood mean
    # elevation. Kept as float so the downstream area mean
    # averages full-precision values (rounding to integer
    # metres happens after aggregation, in R). clip(aoi) crops
    # to the province; the focal mean itself stays unbiased at
    # the edge because it is computed from the buffered
    # elevation (aoi_compute).
    tpi = (
        elevation
        .subtract(
            elevation.focalMean(radius, TPI_WINDOW_SHAPE, TPI_UNITS)
        )
        .clip(aoi)
        .rename(f"tpi_{radius}")
    )

    # 4.1 Export this radius as a 30 m GeoTIFF ----
    # crs_transform pins pixels to the 30 m lattice; scale is
    # intentionally not passed. Large exports are sharded by
    # Earth Engine into multiple tile files that the R step
    # stitches back together.
    task = export_image_to_drive(
        image=tpi,
        description=f"FABDEM_TPI_Alberta_30m_r{radius}",
        region=aoi,
        folder=DRIVE_FOLDER,
        file_name_prefix=f"fabdem_tpi_alberta_30m_r{radius}",
        crs=GRID_CRS,
        crs_transform=EXPORT_CRS_TRANSFORM,
        max_pixels=1e13,
        wait=False,
    )
    tasks.append(task)

# 5. Compute usage report ----
# This section waits for each export to finish, records its
# total EECU-seconds, and writes the txt report to
# gee_compute_reports/. Note: a full-province export can
# take hours; for a quick profile use the test AOI.

for task in tasks:
    report.log_task(task)
report.write()

# End of script ----
