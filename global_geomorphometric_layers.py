# ---
# title:   Global Geomorphometric Layers (Geomorpho90m)
# author:  Brendan Casey
# created: 2026-07-10
# inputs:
#   - Geomorpho90m ImageCollections
#     (projects/sat-io/open-datasets/Geomorpho90m)
#   - FAO GAUL province boundaries
# outputs:
#   - Multiband Geomorpho90m GeoTIFF for Alberta (exported
#     to Google Drive)
# notes:
#   This script loads multiple geomorphometric variables
#   from the Geomorpho90m dataset, mosaics and clips them
#   to the AOI, and combines them into a single multiband
#   image. Map visualization layers from the original GEE
#   JavaScript are dropped.
#
#   Citation:
#   Amatulli, G., McInerney, D., Sethi, T., Strobl, P.,
#   Domisch, S. (2020). Geomorpho90m, empirical evaluation
#   and accuracy assessment of global high-resolution
#   geomorphometric layers. Scientific Data 7(1), 1-18.
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
EXPORT_SCALE = 90  # meters (Geomorpho90m native ~90 m)
EXPORT_CRS = "EPSG:4326"
PRINT_STATS = True  # min/max check (slow for large AOIs)
USE_TEST_AOI = True  # True: small test AOI; False: Alberta
COMPUTE_REPORT = True  # write EECU usage report (txt);
# blocks until the export task finishes

# Base path and Geomorpho90m collections to combine, in the
# order they are stacked into the multiband export image.
BASE_PATH = "projects/sat-io/open-datasets/Geomorpho90m/"
COLLECTION_NAMES = [
    "aspect",           # Aspect
    "aspect-cosine",    # Aspect-Cosine
    "aspect-sine",      # Aspect-Sine
    "convergence",      # Convergence Index
    "cti",              # Compound Topographic Index (CTI)
    "dev-magnitude",    # Deviation Magnitude
    "dev-scale",        # Deviation Scale
    "eastness",         # Eastness
    "elev-stdev",       # Elevation Standard Deviation
    "northness",        # Northness
    "rough-magnitude",  # Multiscale Roughness Magnitude
    "rough-scale",      # Multiscale Roughness Scale
    "roughness",        # Roughness
    "slope",            # Slope
    "spi",              # Stream Power Index
    "tpi",              # Topographic Position Index (TPI)
    "tri",              # Terrain Ruggedness Index (TRI)
    "vrm",              # Vector Ruggedness Measure (VRM)
]

# 1.2 Initialize Earth Engine ----
# Project ID is read from _gee_config.py
initialize_ee()

# 1.3 Set up compute usage report ----
# Profiles EECU usage per section and per export task.
# Best used with USE_TEST_AOI = True to find choke
# points cheaply before a full-province run.
report = ComputeReport(
    "global_geomorphometric_layers",
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

# 3. Geomorpho90m processing ----
# This section loads, mosaics, clips, and renames each
# Geomorpho90m collection, then combines them into a single
# multiband image.


def load_and_process(collection_name, aoi):
    """Load, mosaic, clip, and rename a collection.

    Parameters
    ----------
    collection_name : str
        Geomorpho90m collection short name.
    aoi : ee.Geometry
        Area of interest to clip to.

    Returns
    -------
    ee.Image
        Single-band image renamed to ``collection_name``.
    """
    return (
        ee.ImageCollection(BASE_PATH + collection_name)
        .mosaic()
        .clip(aoi)
        .rename(collection_name)
    )


geomorpho90m = load_and_process(COLLECTION_NAMES[0], aoi)
for name in COLLECTION_NAMES[1:]:
    geomorpho90m = geomorpho90m.addBands(
        load_and_process(name, aoi)
    )

# 3.1 Check min and max values (optional) ----
# Also runs when COMPUTE_REPORT is on: Earth Engine is
# lazy, so the profiler needs an evaluated computation
# (getInfo) to measure per-algorithm EECU usage.
if PRINT_STATS or COMPUTE_REPORT:
    with report.section("Geomorpho90m min/max (reduceRegion)"):
        stats = geomorpho90m.reduceRegion(
            reducer=ee.Reducer.minMax(),
            geometry=aoi,
            scale=EXPORT_SCALE,
            maxPixels=1e13,
            bestEffort=True,
        ).getInfo()
    print("Geomorpho90m min and max values:", stats)

# 4. Export data ----
# This section exports the multiband Geomorpho90m image to
# Google Drive as a GeoTIFF. Set wait=True to block until
# the task finishes; otherwise monitor progress at
# https://code.earthengine.google.com/tasks

task = export_image_to_drive(
    image=geomorpho90m,
    description="Geomorpho90m_Export",
    region=aoi,
    folder=DRIVE_FOLDER,
    file_name_prefix="global_geomorphometric_layers",
    scale=EXPORT_SCALE,
    crs=EXPORT_CRS,
    max_pixels=1e13,
    wait=False,
)

# 5. Compute usage report ----
# This section waits for the export to finish, records
# its total EECU-seconds, and writes the txt report to
# gee_compute_reports/. Note: a full-province export can
# take hours; for a quick profile use the test AOI.

report.log_task(task)
report.write()

# End of script ----
