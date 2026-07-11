# ---
# title:   FABDEM DEV (Deviation from Mean Elevation) for Alberta
# author:  Brendan Casey
# created: 2026-07-11
# inputs:
#   - FABDEM ImageCollection
#     (projects/sat-io/open-datasets/FABDEM)
#   - FAO GAUL province boundaries
# outputs:
#   - One DEV GeoTIFF per focal radius for Alberta
#     (exported to Google Drive)
# notes:
#   This script calculates DEV (deviation from mean
#   elevation) from the FABDEM bare-earth DEM (30 m, forests
#   and buildings removed), following De Reu et al. (2013):
#
#     DEV = (z - mean_z) / SD_z
#
#   where mean_z and SD_z are the mean and standard deviation
#   of elevation within a focal window. The numerator is the
#   Topographic Position Index (TPI); dividing by
#   SD_z standardizes it by local relief, so DEV is expressed
#   in standard-deviation units rather than metres. A 5 m rise
#   on flat muskeg and a 300 m ridge in the Rockies can then
#   score alike, and values are comparable across radii by
#   construction. DEV is computed at each focal radius in
#   DEV_RADII, with one GeoTIFF exported to Google Drive per
#   radius.
#
#   Data citations:
#   Hawker, L., et al. (2022). A 30 m global map of
#   elevation with forests and buildings removed.
#   Environmental Research Letters, 17(2), 024016.
#   doi:10.1088/1748-9326/ac4d4f
#
#   De Reu, J., et al. (2013). Application of the topographic
#   position index to heterogeneous landscapes. Geomorphology,
#   186, 39-49. doi:10.1016/j.geomorph.2012.12.015
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
DEV_RADII = [250, 1000, 2000]  # one export per radius
DEV_WINDOW_SHAPE = "circle"  # "circle" or "square"
DEV_UNITS = "meters"  # "meters" or "pixels"
SD_EPSILON = 0.001  # floor for SD_z to avoid divide-by-zero
USE_TEST_AOI = True  # True: small test AOI; False: Alberta
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
    "fabdem_dev_alberta",
    out_dir=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "gee_compute_reports",
    ),
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

# 3. Prepare the DEM ----
# This section mosaics the FABDEM collection, pins it to a
# metric projection (EPSG:3402, Alberta 10-TM) so the focal
# radius maps to real ground distance, and clips to the AOI.
# The same elevation image feeds every focal radius below.

elevation = (
    ee.ImageCollection("projects/sat-io/open-datasets/FABDEM")
    .mosaic()
    .setDefaultProjection("EPSG:3402", None, 30)
    .clip(aoi)
    .double()
)

# 4. Compute and export DEV per focal radius ----
# For each radius in DEV_RADII, elevation mean and SD are
# reduced over one shared kernel, DEV = (z - mean) / SD is
# formed and exported to Google Drive as its own GeoTIFF.
# Larger radii use bigger kernels and cost proportionally
# more compute; the per-task batch EECU-seconds in the
# report show where. Set wait=True on the export to block;
# otherwise monitor progress at
# https://code.earthengine.google.com/tasks

tasks = []
for radius in DEV_RADII:
    # One kernel drives both the mean and SD reductions, so
    # numerator and denominator share the exact same window.
    if DEV_WINDOW_SHAPE == "circle":
        kernel = ee.Kernel.circle(radius=radius, units=DEV_UNITS)
    elif DEV_WINDOW_SHAPE == "square":
        kernel = ee.Kernel.square(radius=radius, units=DEV_UNITS)
    else:
        raise ValueError(
            f"Unsupported window shape: {DEV_WINDOW_SHAPE}"
        )

    # Mean and standard deviation of elevation in the window
    mean_z = elevation.reduceNeighborhood(
        reducer=ee.Reducer.mean(),
        kernel=kernel,
    )
    sd_z = elevation.reduceNeighborhood(
        reducer=ee.Reducer.stdDev(),
        kernel=kernel,
    )

    # DEV: TPI (z - mean) standardized by local relief (SD).
    # Floor SD at SD_EPSILON so near-flat windows (SD ~ 0) do
    # not blow up the ratio. DEV is unitless (SD units) and
    # continuous, so it is exported as float (no rounding).
    dev = (
        elevation
        .subtract(mean_z)
        .divide(sd_z.max(SD_EPSILON))
        .rename(f"dev_{radius}")
    )

    # 4.1 Export this radius as a GeoTIFF ----
    task = export_image_to_drive(
        image=dev,
        description=f"FABDEM_DEV_Alberta_{radius}m",
        region=aoi,
        folder=DRIVE_FOLDER,
        file_name_prefix=f"fabdem_dev_alberta_{radius}m",
        scale=EXPORT_SCALE,
        crs=EXPORT_CRS,
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
