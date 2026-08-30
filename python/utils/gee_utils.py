# ---
# title:   GEE Utility Functions
# author:  Brendan Casey
# created: 2026-07-10
# notes:
#   Helper functions for running Google Earth Engine from
#   Python (e.g., in VS Code). Includes authentication /
#   initialization, a Drive export wrapper with optional task
#   monitoring, and reference-grid helpers so every script can
#   produce ABMI 1 km grid-aligned exports with one call
#   instead of copy-pasting the grid parameters and aggregation
#   logic. Grid constants live in _gee_config.py.
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
    crs_transform=None,
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
        scale (float): Pixel resolution in meters. Ignored
            when crs_transform is given.
        crs (str): Output coordinate reference system.
        crs_transform (list): Optional 6-element affine
            transform [xScale, xShear, xTranslate, yShear,
            yScale, yTranslate] in the output CRS. When set,
            pixels are pinned to this exact grid (scale is
            not used), so outputs align to a reference grid
            and stack without resampling. region still bounds
            the export; GEE snaps coverage to the transform.
        max_pixels (float): Maximum allowable pixel count.
        wait (bool): If True, block and print task status
            until the export finishes.

    Returns:
        ee.batch.Task: The started export task.
    """
    export_kwargs = dict(
        image=image,
        description=description,
        folder=folder,
        fileNamePrefix=file_name_prefix or description,
        region=region,
        crs=crs,
        maxPixels=max_pixels,
    )
    # crsTransform and scale are mutually exclusive; passing
    # both makes Earth Engine reject the task.
    if crs_transform is not None:
        export_kwargs["crsTransform"] = crs_transform
    else:
        export_kwargs["scale"] = scale

    task = ee.batch.Export.image.toDrive(**export_kwargs)
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


# --- Reference-grid helpers ----------------------------------
# These wrap the shared ABMI 1 km grid workflow so scripts do
# not repeat the grid parameters or the aggregation / export
# logic. Grid constants come from _gee_config.py.


def define_study_area(use_test_aoi=False, buffer_m=0):
    """Return (aoi, aoi_compute) for an aligned export.

    aoi is the export / crop boundary: the small test polygon
    (TEST_AOI_COORDS) when use_test_aoi is True, otherwise the
    AB2020 provincial boundary asset (PROVINCIAL_BOUNDARY_ASSET).
    aoi_compute is aoi grown by buffer_m metres -- a ring so
    neighborhood operations (focal mean, slope's 3x3, ...) stay
    unbiased at the true AOI edge. Clip the DEM to aoi_compute,
    and the final output back to aoi. buffer_m=0 returns
    aoi_compute == aoi.

    Args:
        use_test_aoi (bool): Use the small test polygon.
        buffer_m (float): Compute-ring width in metres; set it
            >= the largest neighborhood reach (e.g. the biggest
            focal radius, or one pixel for slope).

    Returns:
        tuple(ee.Geometry, ee.Geometry): (aoi, aoi_compute).
    """
    from _gee_config import PROVINCIAL_BOUNDARY_ASSET, TEST_AOI_COORDS

    if use_test_aoi:
        aoi = ee.Geometry.Polygon(TEST_AOI_COORDS)
    else:
        aoi = ee.FeatureCollection(
            PROVINCIAL_BOUNDARY_ASSET
        ).geometry()

    if buffer_m:
        # maxError keeps buffering the detailed boundary cheap;
        # the ring is discarded so a coarse approximation is fine.
        aoi_compute = aoi.buffer(buffer_m, max(1.0, buffer_m * 0.05))
    else:
        aoi_compute = aoi
    return aoi, aoi_compute


def native_scale(source, floor_m=0):
    """Nominal scale (m) of a source's native projection.

    Reads the projection from an ee.Image, or from the first
    image of an ee.ImageCollection. Never pass a mosaic:
    mosaic() drops the fixed projection and reports 1 degree
    (~111320 m).

    The value is the true grid spacing, so it is rarely round
    (FABDEM 30.92, MERIT Hydro 92.77, MODIS 463.31). Pass it
    to setDefaultProjection and reduceResolution as-is.

    floor_m guards a computed focal layer: aggregating a 30 m
    base to 1 km over the full province exceeds Earth Engine's
    per-tile reprojection limit, so those layers need a floor
    of ~50 m. Leave it 0 to take the native scale unchanged.

    Blocks on a getInfo call.

    Args:
        source (ee.Image or ee.ImageCollection): Unmosaicked
            source to read the projection from.
        floor_m (float): Lower bound on the returned scale.

    Returns:
        float: Scale in metres, at least floor_m.
    """
    if isinstance(source, ee.ImageCollection):
        source = source.first()
    scale = ee.Image(source).projection().nominalScale().getInfo()
    return max(scale, floor_m)


def fabdem_elevation(aoi_compute, base_m=50):
    """Mosaicked FABDEM elevation pinned to the grid CRS.

    Returns the FABDEM bare-earth DEM (30 m, forests and
    buildings removed) mosaicked, given a fixed base projection
    (GRID_CRS at base_m metres) so neighborhood radii map to
    real ground distance and FABDEM is served from its pyramids,
    then clipped to aoi_compute. base_m is the resolution the
    science is computed at: keep it >= ~50 m for full-province
    1 km aggregation (see to_reference_grid).

    Args:
        aoi_compute (ee.Geometry): Buffered compute AOI.
        base_m (float): Base resolution in metres.

    Returns:
        ee.Image: FABDEM elevation (single band, double).
    """
    from _gee_config import GRID_CRS

    return (
        ee.ImageCollection("projects/sat-io/open-datasets/FABDEM")
        .mosaic()
        .setDefaultProjection(GRID_CRS, None, base_m)
        .clip(aoi_compute)
        .double()
    )


def reference_grid_mask():
    """The ABMI 1 km reference grid as a mask image.

    Reads GRID_ASSET -- the uploaded grid_1km.tif, 1 over the
    grid's Alberta cells and masked elsewhere -- and pins it to
    the grid projection. Use it with updateMask() in place of
    clipping to the provincial boundary: clip() renders a
    fractional coverage mask and drops cells that only clip the
    province edge, leaving them NA in the export even though
    ab_grid() holds them.

    updateMask() also leaves the image footprint intact, where
    clip() shrinks it to the geometry, so a later
    fill_gaps_nearest() can still write into masked cells.

    Returns:
        ee.Image: 1 over the grid's cells, masked elsewhere.
    """
    from _gee_config import GRID_ASSET, GRID_CRS, GRID_CRS_TRANSFORM

    if not GRID_ASSET:
        raise ValueError(
            "GRID_ASSET is not set in _gee_config.py. Upload "
            "sciSpatialR's inst/extdata/grid_1km.tif as an "
            "Earth Engine image asset (no-data value 0) and "
            "set GRID_ASSET to its id."
        )
    return (
        ee.Image(GRID_ASSET)
        .gt(0)
        .reproject(crs=GRID_CRS, crsTransform=GRID_CRS_TRANSFORM)
    )


def to_reference_grid(
    image,
    aoi,
    reducer=None,
    agg_max_pixels=1024,
    round_values=False,
    grid_mask=True,
):
    """Aggregate a computed image onto the ABMI 1 km grid.

    reduceResolution aggregates the image to the 1 km reference
    cells and the result is cast to float32. The export
    (export_to_reference_grid, or export_image_to_drive with the
    grid crs + crsTransform) is what pins the cells to the grid.

    grid_mask=True restricts the result to the reference grid's
    own cells (reference_grid_mask); aoi is then only the export
    region. Set it False to fall back on clipping to aoi, which
    drops cells that merely clip the province edge.

    IMPORTANT: the input image must already be pinned to a
    metric base projection FINER than 1 km via
    setDefaultProjection (fabdem_elevation does this). A 30 m
    base fails on a full-province run with "Reprojection output
    too large" -- filling one 1 km tile forces the base compute
    over ~256 km, ~8600 px at 30 m (over EE's cap). A ~50 m base
    drops that to ~5200 px.

    float32 is deliberate: pixels masked outside aoi export as
    NaN -> NA, whereas an integer output writes them as 0 with
    no nodata flag (read as valid 0 downstream).

    Args:
        image (ee.Image): Image at a metric base projection.
        aoi (ee.Geometry): Boundary to crop the result to when
            grid_mask is False; ignored otherwise.
        reducer (ee.Reducer): Aggregation reducer; defaults to
            ee.Reducer.mean() (area-weighted mean, the right
            downsample for a continuous surface).
        agg_max_pixels (int): Max base pixels aggregated per
            1 km cell (~(1000 / base)^2; 1024 covers a base of
            ~31 m or coarser -- raise it for a finer base).
        round_values (bool): Round to integer before storing
            (e.g. TPI in whole metres). Leave False for
            continuous data (slope degrees, ratios, ...).
        grid_mask (bool): Mask to the reference grid's cells
            instead of clipping to aoi.

    Returns:
        ee.Image: 1 km, float32, on the reference grid's cells.
    """
    reducer = reducer if reducer is not None else ee.Reducer.mean()
    out = image.reduceResolution(
        reducer=reducer, maxPixels=agg_max_pixels
    )
    if round_values:
        out = out.round()
    out = out.toFloat()
    if grid_mask:
        return out.updateMask(reference_grid_mask())
    return out.clip(aoi)


def fill_gaps_nearest(
    image,
    max_px=0,
    reducer=None,
    aoi=None,
):
    """Fill masked cells from their nearest valid neighbours.

    Gap fill for an image already on the ABMI 1 km grid (what
    to_reference_grid returns). Each pass gives a masked cell
    the reducer's result over its 3x3 neighbourhood and leaves
    cells that already hold a value untouched, so the fill grows
    one ring per pass and a gap takes its value from the nearest
    data. Under the default mean a deep fill comes out smoothed,
    each ring averaging the one outside it; ee.Reducer.mode()
    returns exact class codes at any depth. Cells further than
    max_px from any valid cell stay NA.

    Args:
        image (ee.Image): Image on the 1 km reference grid.
        max_px (int): Maximum fill distance in grid cells, one
            neighbourhood pass per cell, closing a hole up to
            2 * max_px cells wide. 0 returns image unchanged.
        reducer (ee.Reducer): Neighbourhood reducer; defaults to
            ee.Reducer.mean(). Pass the reducer the image was
            aggregated with (ee.Reducer.mode() for class codes).
        aoi (ee.Geometry): Boundary to re-clip to, so the fill
            cannot grow past the AOI edge. Skipped when None.

    Returns:
        ee.Image: Image with gaps up to max_px cells filled.
    """
    from _gee_config import GRID_CRS, GRID_CRS_TRANSFORM

    max_px = int(max_px or 0)
    if max_px <= 0:
        return image

    reducer = reducer if reducer is not None else ee.Reducer.mean()
    band_names = image.bandNames()
    kernel = ee.Kernel.square(radius=1, units="pixels")
    filled = image.reproject(
        crs=GRID_CRS, crsTransform=GRID_CRS_TRANSFORM
    )

    # reduceNeighborhood measures the kernel in pixels, so the
    # grid is pinned before each pass.
    for _ in range(max_px):
        # skipMasked=False reduces at masked cells too; the
        # reducer still ignores masked inputs.
        neighbours = filled.reduceNeighborhood(
            reducer=reducer, kernel=kernel, skipMasked=False
        ).rename(band_names)
        # unmask writes only where the mask is exactly zero.
        filled = filled.unmask(neighbours).reproject(
            crs=GRID_CRS, crsTransform=GRID_CRS_TRANSFORM
        )

    if aoi is not None:
        filled = filled.clip(aoi)
    return filled.rename(band_names)


def export_to_reference_grid(
    image,
    aoi,
    description,
    folder="gee_exports",
    file_name_prefix=None,
    aggregate=True,
    reducer=None,
    agg_max_pixels=1024,
    round_values=False,
    fill_gaps_px=0,
    grid_mask=True,
    max_pixels=1e13,
    wait=False,
):
    """Aggregate (optional) and export an image on the ABMI grid.

    Wraps to_reference_grid() and export_image_to_drive() with
    the reference grid's crs + crsTransform, so the GeoTIFF
    lands on the exact 1 km cells and stacks with every other
    grid export. Pass aggregate=False when the image is already
    on the grid (it is then only clipped to aoi).

    Args:
        image (ee.Image): Image to aggregate and export.
        aoi (ee.Geometry): Export region and crop boundary.
        description (str): Task name / default file name.
        folder (str): Google Drive folder.
        file_name_prefix (str): Output file name.
        aggregate (bool): Aggregate to 1 km first (True) or
            export the image as-is on the grid (False).
        reducer (ee.Reducer): Passed to to_reference_grid.
        agg_max_pixels (int): Passed to to_reference_grid.
        round_values (bool): Passed to to_reference_grid.
        fill_gaps_px (int): Nearest-neighbour gap fill applied
            to the gridded image before export, in grid cells
            (see fill_gaps_nearest). 0 disables it.
        grid_mask (bool): Mask to the reference grid's cells
            (reference_grid_mask) instead of clipping to aoi, so
            the export carries exactly the cells ab_grid() holds.
            aoi then only bounds the export region.
        max_pixels (float): Export pixel budget.
        wait (bool): Block until the task finishes.

    Returns:
        ee.batch.Task: The started export task.
    """
    from _gee_config import GRID_CRS, GRID_CRS_TRANSFORM

    grid = reference_grid_mask() if grid_mask else None

    if aggregate:
        out = to_reference_grid(
            image,
            aoi,
            reducer,
            agg_max_pixels,
            round_values,
            grid_mask=grid_mask,
        )
    elif grid is not None:
        out = image.updateMask(grid)
    else:
        out = image.clip(aoi)

    if fill_gaps_px:
        # The grid mask is re-applied after the fill so filled
        # cells cannot spread past the grid's own cells.
        out = fill_gaps_nearest(
            out,
            max_px=fill_gaps_px,
            reducer=reducer,
            aoi=None if grid is not None else aoi,
        )
        if grid is not None:
            out = out.updateMask(grid)

    # A polygon region re-imposes the coverage mask the grid is
    # meant to replace, so the grid path exports over the bounds.
    region = aoi.bounds(1, GRID_CRS) if grid is not None else aoi

    return export_image_to_drive(
        image=out,
        description=description,
        region=region,
        folder=folder,
        file_name_prefix=file_name_prefix,
        crs=GRID_CRS,
        crs_transform=GRID_CRS_TRANSFORM,
        max_pixels=max_pixels,
        wait=wait,
    )


def export_collection_to_reference_grid(
    collection,
    aoi,
    file_name_fn,
    folder="gee_exports",
    reducer=None,
    agg_base_m=None,
    agg_max_pixels=1024,
    round_values=False,
    fill_gaps_px=0,
    grid_mask=True,
    max_pixels=1e13,
):
    """Aggregate each image in a collection onto the ABMI grid.

    The collection version of export_to_reference_grid: it
    iterates the collection client-side and, for each image,
    aggregates it (area mean by default) onto the ABMI 1 km
    reference grid and exports it on the grid's crs +
    crsTransform, so every per-image surface (e.g. one per year)
    stacks with the other 1 km covariates. Mirrors
    export_image_collection (gee_helpers) but snaps to the grid.

    agg_base_m pins each image to that base projection
    (setDefaultProjection at GRID_CRS) before aggregating. Set
    it coarse enough that ~256 km / agg_base_m stays under Earth
    Engine's per-tile reprojection limit (~8600 px) -- e.g.
    ~50 m for a 10-30 m source; it can equal the native scale
    for coarse sources (e.g. 500 m for MODIS). reduceResolution
    needs it to know the input resolution.

    Args:
        collection (ee.ImageCollection): Images to aggregate.
        aoi (ee.Geometry): Crop boundary and export region.
        file_name_fn (callable): Returns a file name per image.
        folder (str): Google Drive folder.
        reducer (ee.Reducer): Aggregation reducer; defaults to
            ee.Reducer.mean(). Use ee.Reducer.mode() for
            categorical layers (e.g. land cover class).
        agg_base_m (float): Base resolution (m) each image is
            pinned to before aggregating (see above).
        agg_max_pixels (int): reduceResolution maxPixels.
        round_values (bool): Round to integer before storing.
        fill_gaps_px (int): Nearest-neighbour gap fill applied
            to each gridded image before export, in grid cells
            (see fill_gaps_nearest). 0 disables it.
        grid_mask (bool): Mask each image to the reference grid's
            cells instead of clipping to aoi (see
            export_to_reference_grid).
        max_pixels (float): Export pixel budget per image.

    Returns:
        list: The started ee.batch.Task export tasks.
    """
    from _gee_config import GRID_CRS, GRID_CRS_TRANSFORM

    grid = reference_grid_mask() if grid_mask else None
    region = aoi.bounds(1, GRID_CRS) if grid is not None else aoi

    col_list = collection.toList(collection.size())
    size = collection.size().getInfo()

    tasks = []
    for i in range(size):
        try:
            img = ee.Image(col_list.get(i))
            file_name = file_name_fn(img)
            if not file_name or not isinstance(file_name, str):
                raise ValueError("Invalid file name generated.")

            if agg_base_m is not None:
                img = img.setDefaultProjection(
                    crs=GRID_CRS, scale=agg_base_m
                )
            gridded = to_reference_grid(
                img,
                aoi,
                reducer,
                agg_max_pixels,
                round_values,
                grid_mask=grid_mask,
            )
            if fill_gaps_px:
                gridded = fill_gaps_nearest(
                    gridded,
                    max_px=fill_gaps_px,
                    reducer=reducer,
                    aoi=None if grid is not None else aoi,
                )
                if grid is not None:
                    gridded = gridded.updateMask(grid)

            task = ee.batch.Export.image.toDrive(
                image=gridded,
                description=file_name,
                folder=folder,
                fileNamePrefix=file_name,
                region=region,
                crs=GRID_CRS,
                crsTransform=GRID_CRS_TRANSFORM,
                maxPixels=max_pixels,
            )
            task.start()
            print(f"Started export task: {file_name}")
            tasks.append(task)
        except Exception as err:
            print(f"Error processing image: {err}")
            continue

    return tasks
