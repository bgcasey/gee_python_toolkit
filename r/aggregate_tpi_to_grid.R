# ---
# title:   Aggregate 30 m TPI tiles to the ABMI 1 km grid
# author:  Brendan Casey
# created: 2026-07-17
# inputs:
#   - A directory of 30 m TPI GeoTIFF tiles exported by
#     python/fabdem_tpi_alberta.py and downloaded from Google
#     Drive (Earth Engine shards large exports into several
#     tiles named <prefix>-<row>-<col>.tif).
#   - The ABMI 1 km reference grid (a .gdb fishnet or a .tif
#     template) that defines the target CRS, cell size, and
#     lattice registration.
#   - The AB2020 provincial boundary (a vector) to crop to.
# outputs:
#   - One 1 km GeoTIFF aligned to the ABMI reference grid,
#     cropped to Alberta, written to out_path.
# notes:
#   Companion to python/fabdem_tpi_alberta.py. That script now
#   exports TPI at native 30 m because aggregating a computed
#   30 m layer to 1 km over the whole province exceeds Earth
#   Engine's per-tile reprojection limit. This function does
#   the aggregation here instead, where terra streams the data
#   from disk and there is no such limit:
#     1. stitch the exported 30 m tiles into one virtual mosaic
#        (terra::vrt -- no data copy),
#     2. aggregate 30 m -> 1 km by area-weighted mean onto the
#        reference grid template (resample method = "average"),
#        which lands every value on the grid's exact cells and
#        handles the non-integer 1000/30 ratio,
#     3. crop / mask to the provincial boundary,
#     4. round to integer metres and write.
#   Rounding to integer metres happens here, after the mean, so
#   the averaging keeps full 30 m precision (matching the
#   original "round after aggregate" intent).
# ---

library(terra) # raster I/O, mosaicking, resampling

# Build a 1 km raster template from the ABMI reference grid.
# Accepts a raster template (.tif) directly, or derives one
# from the vector fishnet: full (unclipped) cells are the
# largest, so any one of them gives a corner on the grid
# lattice and the exact cell size, which the Alberta-clipped
# edge cells do not. The template's extent is the grid's own
# extent expanded outward to that lattice.
build_reference_template <- function(reference_grid,
                                     reference_grid_layer =
                                       "Grid_1KM_revAB2020") {
  is_raster_ref <- tolower(tools::file_ext(reference_grid)) %in%
    c("tif", "tiff")
  if (is_raster_ref) {
    return(rast(reference_grid))
  }

  grid_v <- if (is.na(reference_grid_layer)) {
    vect(reference_grid)
  } else {
    vect(reference_grid, layer = reference_grid_layer)
  }
  cell_area <- if ("SHAPE_Area" %in% names(grid_v)) {
    grid_v$SHAPE_Area
  } else {
    expanse(grid_v)
  }
  fe <- ext(grid_v[which.max(cell_area)])
  cres <- c(xmax(fe) - xmin(fe), ymax(fe) - ymin(fe))
  ge <- ext(grid_v)
  xmn <- xmin(fe) + floor((xmin(ge) - xmin(fe)) / cres[1]) * cres[1]
  xmx <- xmin(fe) + ceiling((xmax(ge) - xmin(fe)) / cres[1]) * cres[1]
  ymn <- ymin(fe) + floor((ymin(ge) - ymin(fe)) / cres[2]) * cres[2]
  ymx <- ymin(fe) + ceiling((ymax(ge) - ymin(fe)) / cres[2]) * cres[2]
  rast(
    xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx,
    resolution = cres, crs = crs(grid_v)
  )
}

# Stitch 30 m TPI tiles, aggregate to the reference grid by
# area mean, and crop to the AOI.
#
# Args:
#   tiles         Directory holding the 30 m tiles, a vector of
#                 tile paths, or a single GeoTIFF.
#   out_path      Output 1 km GeoTIFF path.
#   reference_grid  ABMI grid .gdb (with reference_grid_layer)
#                 or a .tif template.
#   reference_grid_layer  Layer name in the .gdb; NA for a .tif.
#   aoi_boundary  Provincial boundary vector to crop to; NA to
#                 skip cropping.
#   tile_pattern  Regex to pick tiles when `tiles` is a folder.
#   agg_method    resample() method; "average" = area-weighted
#                 mean (the correct downsample for a continuous
#                 surface).
#   round_values  Round the 1 km means to integer metres.
#   datatype      Output data type; defaults to INT2S when
#                 rounding, else FLT4S.
#   overwrite     Overwrite out_path if it exists.
#
# Returns (invisibly) the 1 km SpatRaster.
aggregate_tpi_to_grid <- function(tiles,
                                  out_path,
                                  reference_grid,
                                  reference_grid_layer =
                                    "Grid_1KM_revAB2020",
                                  aoi_boundary = NA,
                                  tile_pattern = "\\.tif$",
                                  agg_method = "average",
                                  round_values = TRUE,
                                  datatype = NULL,
                                  overwrite = TRUE) {
  # 1. Gather the tiles -------------------------------------
  if (length(tiles) == 1 && dir.exists(tiles)) {
    tile_files <- list.files(
      tiles,
      pattern = tile_pattern,
      full.names = TRUE,
      ignore.case = TRUE
    )
  } else {
    tile_files <- tiles
  }
  tile_files <- tile_files[file.exists(tile_files)]
  if (length(tile_files) == 0) {
    stop("No tiles found to aggregate.")
  }
  message("Aggregating ", length(tile_files), " tile(s).")

  # 2. Virtual mosaic (no data copy) ------------------------
  mosaic <- if (length(tile_files) == 1) {
    rast(tile_files)
  } else {
    vrt(tile_files, overwrite = TRUE)
  }
  if (nlyr(mosaic) > 1) {
    mosaic <- mosaic[[1]]
  }

  # 3. Reference-grid template ------------------------------
  template <- build_reference_template(
    reference_grid, reference_grid_layer
  )
  # resample() needs a shared CRS. The export is already in the
  # grid's CRS, so this normally does nothing; project only if
  # they somehow differ.
  if (!same.crs(mosaic, template)) {
    warning("Tiles and grid differ in CRS; reprojecting tiles.")
    mosaic <- project(mosaic, crs(template), method = "bilinear")
  }

  # 4. Aggregate 30 m -> 1 km by area-weighted mean ---------
  # resample onto the template lands values on the grid's exact
  # cells; "average" weights each 30 m pixel by its overlap.
  message("Aggregating to ", paste(res(template), collapse = " x "),
          " m grid by '", agg_method, "' ...")
  agg <- resample(mosaic, template, method = agg_method)

  # 5. Crop / mask to the AOI -------------------------------
  if (length(aoi_boundary) == 1 && !is.na(aoi_boundary)) {
    aoi <- vect(aoi_boundary)
    if (!same.crs(aoi, agg)) {
      aoi <- project(aoi, crs(agg))
    }
    agg <- crop(mask(agg, aoi), aoi)
    message("Cropped to AOI: ", basename(aoi_boundary))
  }

  # 6. Round and write --------------------------------------
  if (round_values) {
    agg <- round(agg)
  }
  if (is.null(datatype)) {
    datatype <- if (round_values) "INT2S" else "FLT4S"
  }
  writeRaster(
    agg,
    out_path,
    overwrite = overwrite,
    datatype = datatype,
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")
  )
  message(
    "Wrote ", out_path, " (",
    paste(dim(agg)[2:1], collapse = " x "), " cells)."
  )

  invisible(agg)
}

# Example (edit paths, then source this file and run):
#
# source("r/aggregate_tpi_to_grid.R")
# aggregate_tpi_to_grid(
#   tiles = "H:/My Drive/gee_exports",
#   out_path = "tmp/fabdem_tpi_alberta_1km_r1000.tif",
#   reference_grid = paste0(
#     "//ABMI-DATA2/science/spatial_data/temp/",
#     "GRID1SQKM_AB2020.gdb"
#   ),
#   reference_grid_layer = "Grid_1KM_revAB2020",
#   aoi_boundary = paste0(
#     "//ABMI-DATA2/science/spatial_data/temp/",
#     "AB2020_provincial_boundary.shp"
#   ),
#   tile_pattern = "fabdem_tpi_alberta_30m_r1000.*\\.tif$"
# )
#
# End of script ----
