# ---
# title:   Check Raster Exports
# author:  Brendan Casey
# created: 2026-07-11
# inputs:
#   - A directory of GeoTIFFs or a single GeoTIFF (the rasters
#     exported by the GEE scripts and downloaded from Google
#     Drive). Set raster_input below.
# outputs:
#   - Console table of per-raster alignment with the ABMI 1 km
#     reference grid: CRS, extent, resolution, and origin.
#   - Console table of per-layer pixel statistics: counts, % NA,
#     min, max, mean, sd, and the stretch quantiles.
#   - One map PNG and one histogram PNG per band, written to
#     out_dir/<file>/map/ and out_dir/<file>/hist/.
# notes:
#   A quick look at exported rasters before they are used. The
#   inspection itself is four calls into sciSpatialR
#   (https://github.com/ABbiodiversity/sciSpatialR):
#
#     check_alignment()  test the geometry against ab_grid()
#     raster_stats()     summarise raster cell values
#     plot_raster()      map with a stretched colour scale and
#                        the provincial boundary
#     plot_hist()        histogram of the cell values
#
#   Install the package with:
#     remotes::install_github("ABbiodiversity/sciSpatialR")
#
#   Alignment is tested against ab_grid() cropped to the raster,
#   so an AOI subset is compared with the grid cells it should
#   cover rather than with the whole province and can match on
#   extent as well. All four properties TRUE means the raster
#   stacks on the grid with no resampling; send it through
#   resample_to_grid() when they are not.
#
#   raster_stats() reports every band of every file. Bands are
#   plotted one at a time rather than as facets of one figure, so
#   each gets its own colour scale and its own PNG, filed under a
#   folder named for the file it came from.
# ---

# 1. Setup ----

## 1.1 Load packages ----
# sciSpatialR - checks, statistics, figures (version: 0.1.0)
# terra       - read each raster for the alignment check
#               (version: 1.9.34)
# ggplot2     - ggsave for the returned plots (version: 4.0.3)
library(sciSpatialR) # check_alignment, stats, figures
library(terra) # raster I/O for the alignment check
library(ggplot2) # saving the returned ggplots

## 1.2 User parameters ----

# Directory holding the GeoTIFFs to check, or a single GeoTIFF.
raster_input <- "H:/My Drive/gee_exports"

# Directory for the PNG figures; created if it does not exist.
out_dir <- "tmp/raster_inspection"

# Empty out_dir before plotting, so figures left by an earlier
# run of a different set of rasters cannot be mistaken for this
# one's. Set to FALSE to add to what is already there.
clear_out_dir <- TRUE

# Boundary drawn on the maps. "alberta" is the packaged
# provincial boundary and sets the map's extent, so a small AOI
# raster shows as a speck within the province; use NULL to let a
# test-AOI raster fill the frame.
map_boundary <- "alberta"

hist_bins <- 50 # number of histogram bins
stretch <- c(0.02, 0.98) # quantiles the map scale is clamped to
png_width <- 9 # PNG width in inches
png_height <- 6 # PNG height in inches
png_dpi <- 200 # PNG resolution

## 1.3 List input rasters ----
# Accept either a directory of rasters or a single .tif file.
if (dir.exists(raster_input)) {
  raster_files <- list.files(
    raster_input,
    pattern = "\\.tif$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(raster_files) == 0) {
    stop("No .tif rasters found in directory: ", raster_input)
  }
} else if (file.exists(raster_input)) {
  raster_files <- raster_input
} else {
  stop("raster_input does not exist: ", raster_input)
}

message("Found ", length(raster_files), " raster(s) to check.")

## 1.4 Clear previous figures ----
# Delete the contents of out_dir, not the folder itself, and only
# when it already exists.
if (clear_out_dir && dir.exists(out_dir)) {
  old_figures <- list.files(out_dir, full.names = TRUE)
  unlink(old_figures, recursive = TRUE)
  message(
    "Cleared ",
    length(old_figures),
    " item(s) from ",
    out_dir
  )
}

# 2. Reference-grid alignment ----
# check_alignment() reports CRS, extent, resolution, and origin;
# all TRUE means the raster stacks on the ABMI 1 km grid without
# resampling. Each raster is compared with the grid cells it
# covers rather than with the whole province: crop() snaps to the
# grid's own cells, so an off-lattice raster still comes back
# with a different extent and is caught.
ref_grid <- ab_grid()
fit_list <- list()

for (raster_path in raster_files) {
  layer_name <- tools::file_path_sans_ext(basename(raster_path))
  r <- rast(raster_path)
  ref_sub <- if (same.crs(r, ref_grid)) {
    crop(ref_grid, ext(r))
  } else {
    ref_grid
  }
  fit_list[[layer_name]] <- check_alignment(
    r,
    ref = ref_sub,
    verbose = FALSE
  )
}

fit_df <- data.frame(
  raster = names(fit_list),
  do.call(rbind, fit_list),
  row.names = NULL
)
print(fit_df)

# 3. Pixel statistics ----
# raster_stats() streams each file from disk, so the whole folder
# is summarised in one call. The stretch quantiles are requested
# as columns so the table shows the range the maps are painted
# over.
stats_df <- raster_stats(
  raster_files,
  quantiles = stretch,
  verbose = FALSE
)
print(stats_df)

# 4. Figures ----
# A map and a histogram for every band, each on its own colour
# scale rather than facetted with the rest of the file. Figures
# are filed under a folder named for the source file, split into
# map/ and hist/, and named for the band. Bands with no valid
# cells are skipped: there is nothing to stretch a colour scale
# or bin. n_valid comes from the table above, whose rows follow
# band order within each file.
for (raster_path in raster_files) {
  file_name <- tools::file_path_sans_ext(basename(raster_path))
  message("Plotting: ", file_name)

  r <- rast(raster_path)
  n_valid <- stats_df$n_valid[
    stats_df$source == basename(raster_path)
  ]

  map_dir <- file.path(out_dir, file_name, "map")
  hist_dir <- file.path(out_dir, file_name, "hist")
  dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(hist_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nlyr(r))) {
    band <- r[[i]]
    # terra names an unnamed band "<file>_<n>", which is the name
    # the figures are filed under.
    layer_name <- names(band)
    if (!nzchar(layer_name)) {
      layer_name <- paste0(file_name, "_", i)
    }

    if (n_valid[i] == 0) {
      message("  Skipping ", layer_name, ": no valid pixels.")
      next
    }

    ggsave(
      file.path(map_dir, paste0(layer_name, "_map.png")),
      plot = plot_raster(
        band,
        stretch = stretch,
        boundary = map_boundary,
        main = layer_name
      ),
      width = png_width,
      height = png_height,
      dpi = png_dpi
    )

    ggsave(
      file.path(hist_dir, paste0(layer_name, "_hist.png")),
      plot = plot_hist(
        band,
        bins = hist_bins,
        main = layer_name
      ),
      width = png_width,
      height = png_height,
      dpi = png_dpi
    )

    message("  Wrote ", layer_name, " map and histogram.")
  }
}

# End of script ----
