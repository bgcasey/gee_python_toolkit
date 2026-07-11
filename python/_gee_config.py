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
