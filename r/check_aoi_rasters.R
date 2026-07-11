# ---
# title:   Check small AOI test rasters
# author:  Brendan Casey
# created: 2026-07-11
# inputs:
#   - A directory of GeoTIFFs or a single GeoTIFF (the
#     small test-AOI rasters exported by the FABDEM
#     scripts and downloaded from Google Drive). Set
#     RASTER_INPUT below.
# outputs:
#   - Optional console print of per-raster pixel statistics
#     (counts, % NA, min, max, mean, sd).
#   - Optional PNG raster plot per raster, coloured with the
#     MetBrewer "Hiroshige" palette, written to out_dir.
#   - Optional PNG NA-mask plot per raster, with missing
#     cells highlighted, written to out_dir.
#   - Optional PNG histogram per raster, bars coloured with
#     the MetBrewer "Hiroshige" palette, written to out_dir.
# notes:
#   A quick QA pass over the small AOI rasters generated for
#   testing before committing to a full-province export. Each
#   check is toggled independently in the setup section so you
#   can run only what you need (e.g. stats without plots).
#
#   The script reads every *.tif in RASTER_INPUT when it is a
#   directory, or checks that single file when it is a
#   GeoTIFF path. For each raster it optionally: (1)
#   computes summary pixel statistics, (2) writes a coloured
#   PNG raster plot, and (3) computes and optionally writes a
#   coloured PNG histogram. NA / masked pixels are dropped
#   before any statistic or histogram is computed.
#
#   Assumptions: rasters are single-band. For multi-band
#   inputs only the first band is checked (noted at runtime).
# ---

# 1. Setup ----

## 1.1 Load packages ----
# terra   - read rasters, summary stats (version: 1.7-78)
# ggplot2 - histogram plots (version: 3.5.1)
# MetBrewer - "Hiroshige" colour palette (version: 0.2.0)
library(terra) # raster I/O and statistics
library(ggplot2) # histogram figures
library(MetBrewer) # Hiroshige palette

## 1.2 User parameters ----
# Toggle each check with a TRUE / FALSE flag, set the input
# and output directories, and tune the histogram resolution.

# Directory holding the small AOI GeoTIFFs to check, or a
# single GeoTIFF file.
raster_input <- "H:/My Drive/gee_exports"

# Directory for raster/histogram PNGs. Created if it does
# not exist and a plot output is enabled.
out_dir <- "tmp/aoi_raster_checks"

# Check toggles ----
compute_stats <- TRUE # counts, % NA, min, max, mean, sd
compute_hist <- TRUE # compute the value histogram
plot_raster <- TRUE # write coloured raster PNGs
plot_na_mask <- TRUE # write NA-mask PNGs
plot_hist <- TRUE # write coloured histogram PNGs
# plot_hist implies the histogram is computed; it does not
# require compute_hist to be TRUE on its own.

# Histogram / plot settings ----
hist_bins <- 50 # number of histogram bins
raster_stretch_quantiles <- c(0.02, 0.98) # clamp plot scale
png_width <- 1800 # PNG width in pixels
png_height <- 1200 # PNG height in pixels
png_dpi <- 200 # PNG resolution

## 1.3 Prepare output directory ----
# Ensure out_dir exists before anything is written.
if (plot_raster || plot_na_mask || plot_hist) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

## 1.4 List input rasters ----
# Gather every GeoTIFF in raster_input. Accept either a
# directory of rasters or a single .tif file.
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
  if (tolower(tools::file_ext(raster_input)) != "tif") {
    stop("raster_input is a file, but not a .tif: ", raster_input)
  }
  raster_files <- raster_input
} else {
  stop("raster_input does not exist: ", raster_input)
}

message("Found ", length(raster_files), " raster(s) to check.")

# 2. Check each raster ----
# Loop over the rasters. For each one, pull the valid (non-NA)
# pixel values once and reuse them for every enabled check so
# the raster is only read a single time.

# Collector for per-raster stats when compute_stats is TRUE.
stats_list <- list()

for (raster_path in raster_files) {
  layer_name <- tools::file_path_sans_ext(basename(raster_path))
  message("Checking: ", layer_name)

  # Read the raster; use only the first band if multi-band.
  r <- rast(raster_path)
  if (nlyr(r) > 1) {
    message(
      "  Note: ",
      nlyr(r),
      " bands found; checking band 1 only."
    )
    r <- r[[1]]
  }

  # Pull values once so missingness and valid-value summaries
  # can both be derived from the same vector.
  vals_all <- values(r, mat = FALSE)
  n_total <- length(vals_all)
  n_na <- sum(is.na(vals_all))
  n_valid <- n_total - n_na
  pct_na <- if (n_total > 0) 100 * n_na / n_total else NA_real_
  vals <- vals_all[!is.na(vals_all)]

  message(
    "  NA pixels: ",
    n_na,
    " / ",
    n_total,
    " (",
    sprintf("%.2f", pct_na),
    "%)"
  )

  if (length(vals) == 0) {
    message("  Skipping: no valid (non-NA) pixels.")
    next
  }

  ## 2.1 Pixel statistics ----
  # Summary statistics over the valid pixels, collected into a
  # single data frame and printed after the loop.
  if (compute_stats) {
    stats_list[[layer_name]] <- data.frame(
      raster = layer_name,
      n_total = n_total,
      n_na = n_na,
      pct_na = pct_na,
      n_valid = n_valid,
      min = min(vals),
      max = max(vals),
      mean = mean(vals),
      sd = sd(vals)
    )
  }

  ## 2.2 Histogram ----
  # Compute the value histogram (used for the console summary
  # and, when enabled, the PNG). plot = FALSE returns the bin
  # breaks and counts without drawing base graphics.
  if (compute_hist || plot_hist) {
    h <- hist(vals, breaks = hist_bins, plot = FALSE)
    if (compute_hist) {
      message(
        "  Histogram: ",
        length(h$counts),
        " bins, peak count ",
        max(h$counts)
      )
    }
  }

  ## 2.3 Coloured raster PNG ----
  # Save a quick raster preview using the same palette as the
  # histogram so the spatial pattern can be checked visually.
  # Constrain the display range to central quantiles so a
  # few extreme values do not wash out the main spatial
  # pattern, while keeping the raster values continuous.
  if (plot_raster) {
    hiroshige <- met.brewer("Hiroshige", type = "continuous")
    stretch_limits <- stats::quantile(
      vals,
      probs = raster_stretch_quantiles,
      na.rm = TRUE,
      names = FALSE
    )
    if (diff(stretch_limits) > 0) {
      message(
        "  Raster stretch: ",
        signif(stretch_limits[1], 5),
        " to ",
        signif(stretch_limits[2], 5)
      )
    } else {
      stretch_limits <- range(vals)
      message(
        "  Raster stretch skipped: plot values are flat."
      )
    }
    raster_png_path <- file.path(
      out_dir,
      paste0(layer_name, "_raster.png")
    )
    grDevices::png(
      filename = raster_png_path,
      width = png_width,
      height = png_height,
      res = png_dpi
    )
    terra::plot(
      r,
      col = grDevices::colorRampPalette(hiroshige)(256),
      range = stretch_limits,
      main = layer_name,
      axes = FALSE,
      box = FALSE,
      legend = TRUE
    )
    grDevices::dev.off()
    message("  Wrote raster plot: ", raster_png_path)
  }

  ## 2.4 NA-mask PNG ----
  # Plot missing cells directly so white areas in the raster
  # preview can be confirmed as masked/NA rather than a
  # display-range artefact.
  if (plot_na_mask) {
    na_mask <- terra::ifel(is.na(r), 1, 0)
    na_mask_png_path <- file.path(
      out_dir,
      paste0(layer_name, "_na_mask.png")
    )
    grDevices::png(
      filename = na_mask_png_path,
      width = png_width,
      height = png_height,
      res = png_dpi
    )
    terra::plot(
      na_mask,
      col = c("grey95", "firebrick2"),
      main = paste(layer_name, "NA mask"),
      axes = FALSE,
      box = FALSE,
      legend = FALSE
    )
    grDevices::dev.off()
    message("  Wrote NA mask: ", na_mask_png_path)
  }

  ## 2.5 Coloured histogram PNG ----
  # ggplot histogram with bars coloured along the "Hiroshige"
  # palette by pixel value (a continuous gradient), saved as a
  # PNG named after the raster.
  if (plot_hist) {
    hiroshige <- met.brewer("Hiroshige", type = "continuous")

    df <- data.frame(value = vals)
    p <- ggplot(df, aes(x = value, fill = after_stat(x))) +
      geom_histogram(bins = hist_bins, colour = NA) +
      scale_fill_gradientn(colours = hiroshige) +
      labs(
        title = layer_name,
        x = "Pixel value",
        y = "Count",
        fill = "Value"
      ) +
      theme_minimal(base_size = 14)

    png_path <- file.path(
      out_dir,
      paste0(layer_name, "_hist.png")
    )
    ggsave(
      png_path,
      plot = p,
      width = png_width / png_dpi,
      height = png_height / png_dpi,
      dpi = png_dpi
    )
    message("  Wrote histogram: ", png_path)
  }
}

# 3. Print stats table ----
# Combine the per-raster stats into one table and print it to
# the console. Only runs when compute_stats is TRUE and at
# least one raster produced statistics.
if (compute_stats && length(stats_list) > 0) {
  stats_df <- do.call(rbind, stats_list)
  rownames(stats_df) <- NULL
  message("Pixel statistics:")
  print(stats_df)
}

# End of script ----
