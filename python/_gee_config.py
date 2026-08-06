# ---
# title:   GEE Configuration
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Shared configuration for all GEE Python scripts in
#   this folder. Set EE_PROJECT to your Google Cloud
#   project ID registered for Earth Engine. Find it in
#   the Code Editor (code.earthengine.google.com, profile
#   icon, top right) or register a project at
#   code.earthengine.google.com/register.
#
#   Scripts pick this up automatically via
#   utils.gee_utils.initialize_ee().
# ---

# Google Cloud project ID registered for Earth Engine
EE_PROJECT = "ee-bgcasey-abmi"

# Default Google Drive folder for exports
DRIVE_FOLDER = "gee_exports"

# --- ABMI 1 km reference grid --------------------------------
# Align every export to this grid so all rasters share the same
# CRS, cell size, and origin and stack without resampling. Used
# by the reference-grid helpers in utils/gee_utils.py.
GRID_CRS = "EPSG:3400"  # NAD83 / Alberta 10-TM (Forest)
# Affine transform [xScale, xShear, originX, yShear, yScale,
# originY]; origin is the grid-aligned top-left corner covering
# Alberta (west 170616.1822, north 6659532.4311).
GRID_CRS_TRANSFORM = [1000, 0, 170616.1822, 0, -1000, 6659532.4311]

# Alberta provincial boundary (Earth Engine table asset) used
# to crop full-province exports to the province outline.
PROVINCIAL_BOUNDARY_ASSET = (
    "projects/ee-bgcasey-abmi/assets/AB2020_provincial_boundary"
)

# Small test AOI (lon/lat ring) for cheap test runs.
TEST_AOI_COORDS = [
    [-113.5, 55.5],  # Top-left corner
    [-113.5, 55.0],  # Bottom-left corner
    [-112.8, 55.0],  # Bottom-right corner
    [-112.8, 55.5],  # Top-right corner
]
