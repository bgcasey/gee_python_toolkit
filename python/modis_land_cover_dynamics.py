# ---
# title:   MODIS Annual Land Cover Dynamics
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - MODIS MCD12Q2 phenology collection
#     (MODIS/061/MCD12Q2)
#   - FAO GAUL province boundaries (Alberta)
# outputs:
#   - Annual multiband phenology GeoTIFFs exported to
#     Google Drive at native (500 m) resolution and at
#     focal scales (0/150/250 m) in EPSG:3978
# notes:
#   Python port of modis_land_cover_dynamics.js for the
#   Earth Engine Python API. Extracts all bands from the
#   MODIS MCD12Q2 phenology product, applies band scaling
#   factors, casts to Float32, and exports annual
#   multiband images. Focal analyses (0/150/250 m) are
#   exported separately in EPSG:3978.
#
#   The original Map.addLayer/Map.setCenter calls, vis
#   parameters, and debug print() blocks are omitted.
#
#   Setup (once):
#     pip install earthengine-api
#     earthengine authenticate
#   Then set EE_PROJECT in _gee_config.py and run.
# ---

import os
import sys

import ee

# Make utils importable regardless of the working
# directory VS Code runs the script from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gee_config import DRIVE_FOLDER
from utils.compute_report import ComputeReport
from utils.gee_helpers import export_image_collection, focal_stats
from utils.gee_utils import initialize_ee

# 1. Setup ----

# 1.1 User parameters ----
MODIS_START_DATE = "2024-01-01"  # phenology year start
MODIS_END_DATE = "2024-12-31"  # phenology year end
EXPORT_SCALE = 500  # native MCD12Q2 resolution (m)
EXPORT_CRS = "EPSG:4326"  # native export CRS
FOCAL_SCALE = 990  # focal export scale (m)
FOCAL_CRS = "EPSG:3978"  # focal export CRS
FOCAL_KERNELS = [150, 250]  # focal radii (m), circle
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
    "modis_land_cover_dynamics",
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

# 3. Load MODIS MCD12Q2 dataset ----
# Loads the phenology collection, tags each image with its
# year, and clips it to the AOI.


def add_year_and_clip(image):
    """Tag an image with its year and clip it to the AOI."""
    year = image.date().format("yyyy")
    return image.set("year", year).clip(aoi)


dataset = (
    ee.ImageCollection("MODIS/061/MCD12Q2")
    .filter(ee.Filter.date(MODIS_START_DATE, MODIS_END_DATE))
    .map(add_year_and_clip)
)

# 3.1 Apply scaling factors to selected bands ----
# EVI minima/amplitudes are scaled by 0.0001 and EVI areas
# by 0.1, overwriting the original band values.


def apply_scaling(image):
    """Scale the EVI phenology bands of a MODIS image."""
    scaled = (
        image.select(["EVI_Minimum_1"])
        .multiply(0.0001)
        .rename("EVI_Minimum_1")
        .addBands(
            image.select(["EVI_Minimum_2"])
            .multiply(0.0001)
            .rename("EVI_Minimum_2")
        )
        .addBands(
            image.select(["EVI_Amplitude_1"])
            .multiply(0.0001)
            .rename("EVI_Amplitude_1")
        )
        .addBands(
            image.select(["EVI_Amplitude_2"])
            .multiply(0.0001)
            .rename("EVI_Amplitude_2")
        )
        .addBands(
            image.select(["EVI_Area_1"])
            .multiply(0.1)
            .rename("EVI_Area_1")
        )
        .addBands(
            image.select(["EVI_Area_2"])
            .multiply(0.1)
            .rename("EVI_Area_2")
        )
    )
    return ee.Image(
        image.addBands(scaled, None, True).copyProperties(
            image, image.propertyNames()
        )
    )


dataset = dataset.map(apply_scaling)

# 3.2 Ensure all bands are Float32 ----
dataset = dataset.map(lambda img: img.toFloat())

# 4. Check bands (optional) ----
# Earth Engine is lazy, so the profiler needs an evaluated
# computation to measure per-algorithm EECU usage.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("MODIS band min/max (reduceRegion)"):
        stats = (
            dataset.first()
            .reduceRegion(
                reducer=ee.Reducer.minMax(),
                geometry=aoi,
                scale=EXPORT_SCALE,
                maxPixels=1e13,
                bestEffort=True,
            )
            .getInfo()
        )
    print("MODIS first-image min/max:", stats)

# 5. Export time series to Google Drive ----
# Exports each image in the collection as a multiband
# GeoTIFF. export_image_collection iterates client-side and
# starts one export task per image.


def modis_file_name(img):
    """File name for the native-resolution export."""
    year = img.date().format("yyyy").getInfo()
    return "MODIS_MCD12Q2_" + year


export_image_collection(
    dataset,
    aoi,
    DRIVE_FOLDER,
    EXPORT_SCALE,
    EXPORT_CRS,
    modis_file_name,
)

# 6. Focal analysis ----
# Exports focal (neighbourhood) statistics at 0/150/250 m
# in EPSG:3978. The 0 m case renames bands with a "_0"
# suffix but applies no smoothing.

# 6.1 Zero-metre focal (no smoothing) ----


def rename_zero_focal(img):
    """Append a "_0" suffix to every band name."""
    new_names = img.bandNames().map(
        lambda name: ee.String(name).cat("_0")
    )
    return img.rename(new_names)


modis_0 = dataset.map(rename_zero_focal)


def modis_file_name_0(img):
    """File name for the 0 m focal export."""
    year = img.get("year").getInfo() or "unknown"
    return "MODIS_MCD12Q2__0_" + str(year)


export_image_collection(
    modis_0,
    aoi,
    DRIVE_FOLDER,
    FOCAL_SCALE,
    FOCAL_CRS,
    modis_file_name_0,
)

# 6.2 Circular focal means (150 m, 250 m) ----
for kernel_size in FOCAL_KERNELS:
    modis_focal = dataset.map(
        lambda img, k=kernel_size: focal_stats(
            img, k, "circle", ["year"]
        )
    )

    def make_focal_file_name(k):
        def focal_file_name(img):
            year = img.get("year").getInfo() or "unknown"
            return "MODIS_MCD12Q2__" + str(k) + "_" + str(
                year
            )

        return focal_file_name

    export_image_collection(
        modis_focal,
        aoi,
        DRIVE_FOLDER,
        FOCAL_SCALE,
        FOCAL_CRS,
        make_focal_file_name(kernel_size),
    )

# 7. Compute usage report ----
# Writes the profiled sections to gee_compute_reports/.
# Collection exports start many batch tasks, so per-task
# EECU totals are not logged here; monitor progress at
# https://code.earthengine.google.com/tasks
report.write()

# End of script ----
