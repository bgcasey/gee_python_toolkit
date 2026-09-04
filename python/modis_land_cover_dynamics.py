# ---
# title:   MODIS Annual Land Cover Dynamics
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - MODIS MCD12Q2 phenology collection
#     (MODIS/061/MCD12Q2)
#   - AB2020 provincial boundary (Earth Engine asset;
#     _gee_config.PROVINCIAL_BOUNDARY_ASSET) for the crop
# outputs:
#   - One annual multiband phenology GeoTIFF per year, either
#     aggregated to the ABMI 1 km reference grid or at
#     NATIVE_SCALE_M, selected by EXPORT_TARGET
#   - Focal (neighbourhood) GeoTIFFs at 0/150/250 m in
#     FOCAL_CRS when RUN_FOCAL is True, exported alongside
#     and not affected by EXPORT_TARGET
# notes:
#   Extracts all bands from the MODIS MCD12Q2 phenology
#   product, applies the EVI band scaling factors, and casts
#   to Float32.
#
#   Every band uses 32767 as its fill value. With 
#   MASK_FILL_VALUES True those pixels are
#   masked before any scaling or rebasing, so fill is never
#   averaged into a coarser cell.
#
#   The phenology date bands are stored as days since
#   1970-01-01. With CONVERT_DATES_TO_DOY True they are
#   rebased to days since January 1 of the image's phenology
#   year, so exported values are ordinal days.
#
#   NATIVE_SCALE_M is MCD12Q2's nominal 500 m; native_scale
#   reports 463.31 for the underlying sinusoidal grid if you
#   prefer it exact.
#
#   Setup (once):
#     pip install earthengine-api
#     earthengine authenticate
#   Then set EE_PROJECT in _gee_config.py and run.
# ---

import math
import os
import sys

import ee

# Make utils importable regardless of the working
# directory VS Code runs the script from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _gee_config import COARSE_SCALE, DRIVE_FOLDER, GRID_CRS
from utils.compute_report import ComputeReport
from utils.gee_helpers import export_image_collection, focal_stats
from utils.gee_utils import (
    define_study_area,
    export_collection_to_reference_grid,
    initialize_ee,
)

# 1. Setup ----

# 1.1 User parameters ----
FILE_PREFIX = "MODIS_MCD12Q2"

SOURCE_ASSET = "MODIS/061/MCD12Q2"
MODIS_START_DATE = "2001-01-01"  # phenology year start
MODIS_END_DATE = "2024-12-31"  # phenology year end

# MCD12Q2 is already coarse, so the aggregation base is its
# native resolution.
BASE_SCALE_M = 500
NATIVE_SCALE_M = 500  # resolution of the "native" export

# True: mask the product fill value before scaling and
# rebasing, so it is never read as data.
MASK_FILL_VALUES = True

# True: rebase the date bands from days since 1970-01-01 to
# days since January 1 of the image's phenology year.
CONVERT_DATES_TO_DOY = True

# Separate focal product (section 5), not aggregated to the
# ABMI grid and not affected by EXPORT_TARGET.
RUN_FOCAL = False  # False: skip the focal exports
FOCAL_SCALE = 1000  # focal export scale (m)
FOCAL_CRS = "EPSG:3978"  # focal export CRS
FOCAL_KERNELS = [150, 250]  # focal radii (m), circle

FOCAL_REACH_M = max(FOCAL_KERNELS) if RUN_FOCAL else 0

# Reach of the nearest-neighbour gap fill applied after
# aggregation, in 1 km cells; 0 disables it.
FILL_GAPS_PX = 20

# "native" skips the aggregation and writes NATIVE_SCALE_M
# pixels in the grid CRS, ungridded - useful for inspecting
# the input to the aggregation.
EXPORT_TARGET = "reference_grid"  # "native" or "reference_grid"
    
USE_TEST_AOI = False  # True: small test AOI; False: Alberta
PRINT_STATS = False  # value preview (slow for large AOIs)
COMPUTE_REPORT = False  # write EECU usage report (txt)
# Costs the full export runtime (hours province-wide), so
# turn it on only when profiling the test AOI.
WAIT_FOR_EXPORTS = False

# Derived: keeps test-AOI tasks and files
# distinguishable from full-extent runs.
TEST_SUFFIX = "_test" if USE_TEST_AOI else ""

# 1.2 Validate parameters ----
if EXPORT_TARGET not in ("native", "reference_grid"):
    raise ValueError(
        "Unknown EXPORT_TARGET: "
        f"{EXPORT_TARGET!r} (use 'native' or 'reference_grid')"
    )
if FILL_GAPS_PX < 0:
    raise ValueError(
        f"FILL_GAPS_PX must be >= 0, got {FILL_GAPS_PX!r}"
    )

# 1.3 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.4 Derived scales ----
# AGG_BUFFER_M is 2x the output scale, so a 1 km cell touching
# the aoi at a corner still reads a full diagonal (1414 m)
# beyond it; AGG_MAX_PIXELS is reduceResolution's per-cell
# input budget, +2 covering a cell that straddles a base pixel
# on each side.
AGG_BUFFER_M = 2 * (
    COARSE_SCALE
    if EXPORT_TARGET == "reference_grid"
    else NATIVE_SCALE_M
)
COMPUTE_BUFFER_M = max(FOCAL_REACH_M, AGG_BUFFER_M)
AGG_MAX_PIXELS = (
    math.ceil(COARSE_SCALE / BASE_SCALE_M) + 2
) ** 2

# 1.5 Set up run bookkeeping ----
# report profiles compute usage; export_tasks collects the
# export tasks it logs.
report = ComputeReport(FILE_PREFIX, enabled=COMPUTE_REPORT)
export_tasks = []

# 2. Define study area ----
# aoi is the export / crop boundary; aoi_compute adds the ring
# the source is read over.
aoi, aoi_compute = define_study_area(
    use_test_aoi=USE_TEST_AOI,
    buffer_m=COMPUTE_BUFFER_M,
)

# 3. Build the collection ----
# Reads over aoi_compute, never aoi; clipped back in section 4.


def add_year_and_clip(image):
    """Tag an image with its year and clip it to the ring."""
    year = image.date().format("yyyy")
    return image.set("year", year).clip(aoi_compute)


# Fill value shared by every MCD12Q2 band (C6.1 user guide,
# Table 1).
FILL_VALUE = 32767


def mask_fill(image):
    """Mask the fill value in every band."""
    return ee.Image(
        image.updateMask(image.neq(FILL_VALUE)).copyProperties(
            image, image.propertyNames()
        )
    )


def apply_scaling(image):
    """Rescale the EVI phenology bands to physical units."""
    scale_by = {
        "EVI_Minimum_1": 0.0001,
        "EVI_Minimum_2": 0.0001,
        "EVI_Amplitude_1": 0.0001,
        "EVI_Amplitude_2": 0.0001,
        "EVI_Area_1": 0.1,
        "EVI_Area_2": 0.1,
    }
    scaled = ee.Image.cat(
        [
            image.select([band]).multiply(factor).rename(band)
            for band, factor in scale_by.items()
        ]
    )
    return ee.Image(
        image.addBands(scaled, None, True).copyProperties(
            image, image.propertyNames()
        )
    )


# Bands holding a phenological transition date, stored as
# days since 1970-01-01.
DATE_BANDS = [
    "Greenup_1",
    "Greenup_2",
    "MidGreenup_1",
    "MidGreenup_2",
    "Peak_1",
    "Peak_2",
    "Maturity_1",
    "Maturity_2",
    "Senescence_1",
    "Senescence_2",
    "MidGreendown_1",
    "MidGreendown_2",
    "Dormancy_1",
    "Dormancy_2",
]


def dates_to_day_of_year(image):
    """Rebase the date bands to the image's phenology year."""
    year_start = ee.Date.fromYMD(image.date().get("year"), 1, 1)
    offset = year_start.difference(ee.Date("1970-01-01"), "day")
    date_bands = image.bandNames().filter(
        ee.Filter.inList("item", DATE_BANDS)
    )
    rebased = image.select(date_bands).subtract(offset)
    return ee.Image(
        image.addBands(rebased, None, True).copyProperties(
            image, image.propertyNames()
        )
    )


dataset = (
    ee.ImageCollection(SOURCE_ASSET)
    .filter(ee.Filter.date(MODIS_START_DATE, MODIS_END_DATE))
    .map(add_year_and_clip)
)

if MASK_FILL_VALUES:
    dataset = dataset.map(mask_fill)

dataset = dataset.map(apply_scaling)

if CONVERT_DATES_TO_DOY:
    dataset = dataset.map(dates_to_day_of_year)

dataset = dataset.map(lambda img: img.toFloat())

# 3.1 Check layer values (optional) ----
# Earth Engine is lazy, so the profiler needs an evaluated
# computation to measure EECU usage; this also runs when
# COMPUTE_REPORT is on.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("MODIS band min/max (reduceRegion)"):
        stats = (
            dataset.first()
            .reduceRegion(
                reducer=ee.Reducer.minMax(),
                geometry=aoi,
                scale=NATIVE_SCALE_M,
                maxPixels=1e13,
                bestEffort=True,
            )
            .getInfo()
        )
    print("MODIS first-image min/max:", stats)

# 4. Export the time series ----
# One export task per image. Monitor progress at
# https://code.earthengine.google.com/tasks
target_suffix = (
    "abmi1km" if EXPORT_TARGET == "reference_grid" else "native"
) + TEST_SUFFIX


def modis_file_name(img):
    """File name for one multiband export."""
    year = img.date().format("yyyy").getInfo()
    return f"{FILE_PREFIX}_{year}_{target_suffix}"


if EXPORT_TARGET == "reference_grid":
    export_tasks += export_collection_to_reference_grid(
        dataset,
        aoi,
        modis_file_name,
        folder=DRIVE_FOLDER,
        reducer=ee.Reducer.mean(),
        agg_base_m=BASE_SCALE_M,
        agg_max_pixels=AGG_MAX_PIXELS,
        fill_gaps_px=FILL_GAPS_PX,
    )
else:
    export_tasks += export_image_collection(
        dataset,
        aoi,
        DRIVE_FOLDER,
        NATIVE_SCALE_M,
        GRID_CRS,
        modis_file_name,
    )

# 5. Focal analysis ----
# An optional separate product, controlled by RUN_FOCAL: focal
# (neighbourhood) statistics at 0/150/250 m in FOCAL_CRS. The
# 0 m case appends a "_0" band suffix but applies no smoothing.

if RUN_FOCAL:

    def make_focal_file_name(kernel_size):
        """Build the file-name function for one focal radius."""

        def focal_file_name(img):
            year = img.get("year").getInfo() or "unknown"
            return f"{FILE_PREFIX}_{kernel_size}_{year}"

        return focal_file_name

    def rename_zero_focal(img):
        """Append a "_0" suffix to every band name."""
        new_names = img.bandNames().map(
            lambda name: ee.String(name).cat("_0")
        )
        return img.rename(new_names)

    export_tasks += export_image_collection(
        dataset.map(rename_zero_focal),
        aoi,
        DRIVE_FOLDER,
        FOCAL_SCALE,
        FOCAL_CRS,
        make_focal_file_name(0),
    )

    for kernel_size in FOCAL_KERNELS:
        modis_focal = dataset.map(
            lambda img, k=kernel_size: focal_stats(
                img, k, "circle", ["year"]
            )
        )
        export_tasks += export_image_collection(
            modis_focal,
            aoi,
            DRIVE_FOLDER,
            FOCAL_SCALE,
            FOCAL_CRS,
            make_focal_file_name(kernel_size),
        )

# 6. Compute usage report ----
# Records each task's total EECU-seconds and writes the txt
# report to gee_compute_reports/.

if WAIT_FOR_EXPORTS:
    for task in export_tasks:
        report.log_task(task)
report.write()

# End of script ----
