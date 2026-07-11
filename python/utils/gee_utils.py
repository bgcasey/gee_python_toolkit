# ---
# title:   GEE Utility Functions
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Helper functions for running Google Earth Engine from
#   Python (e.g., in VS Code). Includes authentication /
#   initialization and a Drive export wrapper with optional
#   task monitoring.
# ---

import time

import ee


def initialize_ee(project=None):
    """Authenticate and initialize the Earth Engine API.

    Reads the project ID from _gee_config.py unless one is
    passed explicitly. Tries to initialize with existing
    credentials first. If that fails, runs the interactive
    authentication flow (opens a browser) and initializes
    again.

    Args:
        project (str): Google Cloud project ID registered
            for Earth Engine. Overrides the EE_PROJECT
            value in _gee_config.py.
    """
    if project is None:
        from _gee_config import EE_PROJECT
        project = EE_PROJECT

    if not project or project == "ee-your-project-id":
        raise ValueError(
            "Set EE_PROJECT in _gee_config.py to your "
            "registered Earth Engine cloud project ID "
            "(see code.earthengine.google.com, profile "
            "icon, top right)."
        )

    try:
        ee.Initialize(project=project)
    except Exception:
        ee.Authenticate()
        ee.Initialize(project=project)
    print("Earth Engine initialized.")


def export_image_to_drive(
    image,
    description,
    region,
    folder="gee_exports",
    file_name_prefix=None,
    scale=30,
    crs="EPSG:4326",
    max_pixels=1e13,
    wait=False,
):
    """Export an ee.Image to Google Drive as a GeoTIFF.

    Args:
        image (ee.Image): Image to export.
        description (str): Task name shown in the Task list.
        region (ee.Geometry): Export region.
        folder (str): Google Drive folder name.
        file_name_prefix (str): Output file name. Defaults
            to the task description.
        scale (float): Pixel resolution in meters.
        crs (str): Output coordinate reference system.
        max_pixels (float): Maximum allowable pixel count.
        wait (bool): If True, block and print task status
            until the export finishes.

    Returns:
        ee.batch.Task: The started export task.
    """
    task = ee.batch.Export.image.toDrive(
        image=image,
        description=description,
        folder=folder,
        fileNamePrefix=file_name_prefix or description,
        region=region,
        scale=scale,
        crs=crs,
        maxPixels=max_pixels,
    )
    task.start()
    print(f"Started export task: {description}")

    if wait:
        monitor_task(task)

    return task


def monitor_task(task, poll_interval=30):
    """Poll an export task and print status until it ends.

    Args:
        task (ee.batch.Task): Task to monitor.
        poll_interval (int): Seconds between status checks.
    """
    while task.active():
        status = task.status()
        print(f"  Task {status['state']}...")
        time.sleep(poll_interval)

    status = task.status()
    print(f"  Task finished with state: {status['state']}")
    if status["state"] == "FAILED":
        print(f"  Error: {status.get('error_message')}")
