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
#   - One TPI GeoTIFF per focal radius for Alberta, either
#     aligned to the ABMI 1 km reference grid or at
#     FOCAL_BASE_M resolution, selected by EXPORT_TARGET
#     (exported to Google Drive)
# notes:
#   This script calculates the Topographic Position Index
#   (TPI) from the FABDEM bare-earth DEM (30 m, forests and
#   buildings removed): elevation minus the mean elevation of
#   a surrounding neighborhood. The collection is mosaicked,
#   the AOI is defined, and TPI is computed at each focal
#   radius in TPI_RADII at FOCAL_BASE_M resolution.
#
#   Each TPI surface is then aggregated to the 1 km ABMI
#   reference grid by area mean (reduceResolution) and exported
#   on that grid's exact CRS and affine transform (crsTransform,
#   not scale), so every output shares identical spatial
#   properties -- CRS, cell size, and origin -- and stacks
#   without any further interpolation or alignment.
#
#   TPI is computed at FOCAL_BASE_M (50 m) rather than the
#   native 30 m. Aggregating a computed 30 m layer straight to
#   1 km over the whole province exceeds Earth Engine's
#   per-tile reprojection limit ("Reprojection output too
#   large"): filling one 1 km output tile forces the focal mean
#   over a ~256 km footprint, which at 30 m is ~8600 px per
#   side -- over the ~8192 cap. A 50 m base drops it to
#   ~5200 px. FABDEM is served from its pyramids at 50 m (an
#   area mean of the 30 m elevation), so for radii >> the base
#   the coarser base costs almost nothing.
#
#   Reference grid:
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
# Output grid: exports are pinned to the ABMI 1 km reference
# grid so every raster shares its exact spatial properties and
# they stack without further processing. GRID_CRS_TRANSFORM is
# [xScale, xShear, originX, yShear, yScale, originY]; the origin
# is the grid-aligned top-left corner covering Alberta (west
# 170616.1822, north 6659532.4311), yScale negative because
# rows run north -> south.
GRID_CRS = "EPSG:3400"  # NAD83 / Alberta 10-TM (Forest)
GRID_CRS_TRANSFORM = [1000, 0, 170616.1822, 0, -1000, 6659532.4311]
# Focal base resolution (metres). # Focal-base pixels averaged into 
# each 1 km cell TPI is computed at this resolution, not the native 
# 30 m, so the jump to 1 km stays under Earth Engine's per-tile 
# reprojection limit. 
FOCAL_BASE_M = 40
AGG_MAX_PIXELS = 1024
# Full-province runs are cropped to this boundary (an Earth
# Engine table asset). Upload AB2020_provincial_boundary.shp
# to your EE assets and set the ID here. Only used when
# USE_TEST_AOI is False.
PROVINCIAL_BOUNDARY_ASSET = (
    "projects/ee-bgcasey-abmi/assets/AB2020_provincial_boundary"
)
TPI_RADII = [250]  # one export per radius
# Raster export target, applied to every radius.
# "reference_grid" aggregates TPI (area mean) onto the ABMI 1 km
# grid so it stacks with the other 1 km covariates. "native"
# skips the aggregation and exports TPI at FOCAL_BASE_M in the
# grid CRS, ungridded - useful for inspecting the input to the
# aggregation.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"

if EXPORT_TARGET not in ("native", "reference_grid"):
    raise ValueError(
        "Unknown EXPORT_TARGET: "
        f"{EXPORT_TARGET!r} (use 'native' or 'reference_grid')"
    )
TPI_WINDOW_SHAPE = "circle"  # "circle" or "square"
TPI_UNITS = "meters"  # "meters" or "pixels"
# Compute buffer (metres). The DEM is clipped to the AOI grown
# by this ring before the focal mean, so cells at the true AOI
# edge see a full neighborhood of real elevation instead of
# masked pixels (which would bias TPI inward). The ring is
# discarded when the output is clipped back to the AOI, and
# compute stays bounded to AOI + ring rather than the whole
# DEM. It must be >= the largest focal reach; focal-base pixels
# are FOCAL_BASE_M, so "pixels" radii are converted to metres.
COMPUTE_BUFFER_M = max(TPI_RADII) * (
    FOCAL_BASE_M if TPI_UNITS == "pixels" else 1
)
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
# This section mosaics the FABDEM collection and pins it to the
# reference grid's metric projection (EPSG:3400, Alberta 10-TM)
# at FOCAL_BASE_M so the focal radius maps to real ground
# distance and the aggregation stays within one CRS. At this
# scale FABDEM is served from its pyramids (an area mean of the
# 30 m elevation), which is what keeps the per-tile footprint
# under Earth Engine's reprojection limit. It clips to the
# buffered AOI (aoi_compute) so the focal mean has real
# elevation on all sides of the true AOI edge. The same
# elevation image feeds every focal radius below.

elevation = (
    ee.ImageCollection("projects/sat-io/open-datasets/FABDEM")
    .mosaic()
    .setDefaultProjection(GRID_CRS, None, FOCAL_BASE_M)
    .clip(aoi_compute)
    .double()
)

# 4. Compute, aggregate, and export TPI per focal radius ----
# For each radius in TPI_RADII, TPI is elevation minus the
# neighborhood mean elevation at FOCAL_BASE_M, then aggregated
# to the 1 km reference grid by area mean and exported to
# Google Drive as its own GeoTIFF on the grid's exact CRS and
# transform.
# Larger radii use bigger focal kernels and cost
# proportionally more compute; the per-task batch
# EECU-seconds in the report show where.
# Set wait=True on the export to block; otherwise monitor
# progress at https://code.earthengine.google.com/tasks

tasks = []
for radius in TPI_RADII:
    # TPI at FOCAL_BASE_M: elevation minus neighborhood mean
    # elevation. Kept as float here so the area mean below
    # averages full-precision values.
    tpi = elevation.subtract(
        elevation.focalMean(radius, TPI_WINDOW_SHAPE, TPI_UNITS)
    )

    # The "reference_grid" path (EXPORT_TARGET, validated at the
    # top of the script); the "native" path below skips all of
    # this and exports tpi at FOCAL_BASE_M instead.
    # Aggregate FOCAL_BASE_M -> 1 km by area mean, pinned to the
    # ABMI grid. reduceResolution averages every focal-base
    # pixel that falls in each 1 km cell; the aggregation lands
    # on the reference grid's exact cells because the export
    # requests the output in that grid's crs + crsTransform
    # (section 4.1) -- no explicit reproject() is used (over the
    # whole province it fails with "Reprojection output too
    # large"). round() stores TPI as integer metres (rounding
    # after the mean so averaging keeps full precision; it
    # rounds symmetrically to nearest, correct for negative
    # valley values -- drop it for a continuous raster).
    # toFloat() keeps float32 storage on purpose: clip(aoi) masks
    # the pixels outside Alberta, and Earth Engine writes masked
    # *integer* pixels as 0 with no nodata flag (so downstream
    # tools read the out-of-province background as a valid 0),
    # whereas masked *float* pixels export as NaN, which GDAL /
    # terra read as NA. clip(aoi) crops to the true AOI,
    # discarding the buffer ring used only to keep the focal
    # mean unbiased.
    if EXPORT_TARGET == "reference_grid":
        tpi_1km = (
            tpi
            .reduceResolution(
                reducer=ee.Reducer.mean(),
                maxPixels=AGG_MAX_PIXELS,
            )
            .round()
            .toFloat()
            .clip(aoi)
            .rename(f"tpi_{radius}")
        )

        # 4.1 Export this radius as a 1 km GeoTIFF ----
        # crs_transform pins pixels to the reference grid; scale
        # is intentionally not passed.
        task = export_image_to_drive(
            image=tpi_1km,
            description=f"FABDEM_TPI_Alberta_1km_r{radius}",
            region=aoi,
            folder=DRIVE_FOLDER,
            file_name_prefix=f"fabdem_tpi_alberta_1km_r{radius}",
            crs=GRID_CRS,
            crs_transform=GRID_CRS_TRANSFORM,
            max_pixels=1e13,
            wait=False,
        )
    else:
        # 4.1 Export this radius at FOCAL_BASE_M ----
        # No reduceResolution and no round(): this is the raw
        # input to the aggregation, so it keeps sub-metre TPI
        # values (the 1 km path rounds only after averaging).
        # toFloat() for the same nodata reason as above, and
        # clip(aoi) discards the compute buffer ring.
        task = export_image_to_drive(
            image=(
                tpi.toFloat().clip(aoi).rename(f"tpi_{radius}")
            ),
            description=(
                f"FABDEM_TPI_Alberta_{FOCAL_BASE_M}m_r{radius}"
            ),
            region=aoi,
            folder=DRIVE_FOLDER,
            file_name_prefix=(
                f"fabdem_tpi_alberta_{FOCAL_BASE_M}m_r{radius}"
            ),
            scale=FOCAL_BASE_M,
            crs=GRID_CRS,
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
