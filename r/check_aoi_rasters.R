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
#   - Optional console print of a per-raster reference-grid
#     fit check (CRS, resolution, lattice alignment, whether
#     the raster stacks onto the reference grid, and whether
#     all valid pixels fall inside the AOI boundary).
#   - Optional PNG raster plot per raster, coloured with the
#     MetBrewer "Hiroshige" palette, with the AOI boundary
#     overlaid, written to out_dir.
#   - Optional PNG NA-mask plot per raster, with missing
#     cells highlighted and the AOI boundary overlaid, written
#     to out_dir.
#   - Optional PNG histogram per raster, bars coloured with
#     the MetBrewer "Hiroshige" palette, written to out_dir.
#   - Optional PNG grid-stack plot per raster: the raster
#     beside the dummy reference grid it stacks with (a
#     checkerboard of 1 km cells cropped to the AOI), written
#     to out_dir.
# notes:
#   A quick QA pass over the small AOI rasters generated for
#   testing before committing to a full-province export. Each
#   check is toggled independently in the setup section so you
#   can run only what you need (e.g. stats without plots).
#
#   The script reads every *.tif in RASTER_INPUT when it is a
#   directory, or checks that single file when it is a
#   GeoTIFF path. For each raster it optionally: (1) checks
#   that it fits the reference grid, (2) computes summary
#   pixel statistics, (3) writes a coloured PNG raster plot,
#   and (4) computes and optionally writes a coloured PNG
#   histogram. NA / masked pixels are dropped before any
#   statistic or histogram is computed.
#
#   The reference-grid fit check reads the grid the exports
#   are aligned to (the ABMI 1 km fishnet by default), derives
#   a raster template from its CRS, cell size, and lattice
#   registration, and confirms each raster shares them and
#   stacks onto the template without resampling. A raster that
#   stacks can be combined directly with other grid-aligned
#   layers. The stack test crops that template to the raster,
#   fills it with a dummy checkerboard, and crops it to the
#   AOI boundary, so the alignment proof doubles as a
#   visualization. It also reports whether any valid pixels
#   fall outside the AOI boundary (i.e. whether the raster is
#   cropped to it).
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
check_grid_fit <- TRUE # verify raster aligns to and stacks
#                        on the reference grid, and is cropped
#                        to the AOI boundary
compute_stats <- TRUE # counts, % NA, min, max, mean, sd
compute_hist <- TRUE # compute the value histogram
plot_raster <- TRUE # write coloured raster PNGs
plot_na_mask <- TRUE # write NA-mask PNGs
plot_hist <- TRUE # write coloured histogram PNGs
plot_grid_stack <- TRUE # write raster + dummy-grid stack PNGs
#                         (needs check_grid_fit = TRUE)
overlay_boundary <- TRUE # draw the AOI boundary on the raster
#                          and NA-mask PNGs
# plot_hist implies the histogram is computed; it does not
# require compute_hist to be TRUE on its own.

# Histogram / plot settings ----
hist_bins <- 50 # number of histogram bins
raster_stretch_quantiles <- c(0.02, 0.98) # clamp plot scale
png_width <- 1800 # PNG width in pixels
png_height <- 1200 # PNG height in pixels
png_dpi <- 200 # PNG resolution

# Reference-grid fit settings ----
# The grid every export is aligned to. Either a vector fishnet
# (a .gdb with reference_grid_layer set) or a raster template
# (a .tif, with reference_grid_layer = NA). The check derives
# the target CRS, cell size, and lattice registration from it.
reference_grid <- paste0(
  "//ABMI-DATA2/science/spatial_data/temp/",
  "GRID1SQKM_AB2020.gdb"
)
reference_grid_layer <- "Grid_1KM_revAB2020" # NA for a .tif
grid_fit_tolerance <- 1e-4 # metres of slack for alignment

# AOI boundary settings ----
# The provincial boundary the exports are cropped to (the
# "AOI"). Used to overlay on plots, check each raster's valid
# pixels fall inside it, and crop the dummy reference grid in
# the stack visualization. Any vector terra can read (.shp,
# .gpkg, a .gdb layer via a second argument, etc.).
aoi_boundary <- paste0(
  "//ABMI-DATA2/science/spatial_data/temp/",
  "AB2020_provincial_boundary.shp"
)

## 1.3 Prepare output directory ----
# Ensure out_dir exists before anything is written.
if (plot_raster || plot_na_mask || plot_hist || plot_grid_stack) {
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

## 1.5 Load reference grid template ----
# When check_grid_fit is on, load the grid the exports are
# aligned to and derive a raster template from it: its CRS,
# cell size, and lattice registration. Each raster is later
# tested against this template. Reading the fishnet is the slow
# step, so it happens once here. A read/derivation failure
# disables the check with a warning rather than aborting the
# whole QA run.
ref_template <- NULL
ref_res_xy <- NULL
ref_anchor <- NULL

if (check_grid_fit) {
  ref_template <- tryCatch(
    {
      is_raster_ref <- tolower(tools::file_ext(reference_grid)) %in%
        c("tif", "tiff")
      if (is_raster_ref) {
        # A raster already defines the target grid directly.
        rast(reference_grid)
      } else {
        # Vector fishnet: read it, then derive a template from
        # its CRS, cell size, and registration. Full (unclipped)
        # cells have the largest area; any one of them gives a
        # corner on the grid lattice and the exact cell size,
        # which edge cells (clipped to the AOI) do not.
        message("Reading reference grid: ", reference_grid)
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
        full_cell <- grid_v[which.max(cell_area)]
        fe <- ext(full_cell)
        cres <- c(xmax(fe) - xmin(fe), ymax(fe) - ymin(fe))
        # Expand the grid's own extent outward to the lattice
        # anchored on the full-cell corner, so the template
        # covers every cell with the correct registration.
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
    },
    error = function(e) {
      warning(
        "Could not load reference grid; skipping grid-fit ",
        "check. Reason: ", conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(ref_template)) {
    check_grid_fit <- FALSE
  } else {
    ref_res_xy <- res(ref_template)
    ref_anchor <- c(xmin(ref_template), ymin(ref_template))
    ref_name <- tryCatch(
      crs(ref_template, describe = TRUE)$name,
      error = function(e) NA_character_
    )
    message(
      "Reference grid: ",
      if (is.na(ref_name)) "custom CRS" else ref_name,
      ", res ", paste(signif(ref_res_xy, 6), collapse = " x "),
      ", lattice origin (",
      format(ref_anchor[1], nsmall = 2), ", ",
      format(ref_anchor[2], nsmall = 2), ")"
    )
  }
}

## 1.6 Load AOI boundary ----
# Read the provincial boundary once for plot overlays, the
# crop-coverage check, and cropping the dummy reference grid.
# Loaded only when a feature that needs it is enabled. A read
# failure disables those features with a warning rather than
# aborting the run.
aoi_vect <- NULL
if (check_grid_fit || overlay_boundary || plot_grid_stack) {
  aoi_vect <- tryCatch(
    vect(aoi_boundary),
    error = function(e) {
      warning(
        "Could not read AOI boundary; overlay, crop check, ",
        "and grid stack are disabled. Reason: ",
        conditionMessage(e)
      )
      NULL
    }
  )
  if (!is.null(aoi_vect)) {
    message("Loaded AOI boundary: ", basename(aoi_boundary))
  }
}

# Return the AOI in a target object's CRS, projecting only
# when they differ, so overlays and masks line up.
aoi_in_crs <- function(target) {
  if (is.null(aoi_vect)) {
    return(NULL)
  }
  if (terra::same.crs(aoi_vect, target)) {
    aoi_vect
  } else {
    terra::project(aoi_vect, crs(target))
  }
}

# 2. Check each raster ----
# Loop over the rasters. For each one, pull the valid (non-NA)
# pixel values once and reuse them for every enabled check so
# the raster is only read a single time.

# Collectors for per-raster stats and grid-fit results.
stats_list <- list()
grid_fit_list <- list()

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

  # Raster + dummy reference grid, built by the stack test
  # below and reused by the grid-stack plot. NULL when the
  # stack could not be built.
  grid_stack <- NULL

  ## 2.1 Reference-grid fit ----
  # Confirm the raster shares the reference grid's CRS, cell
  # size, and lattice registration, then prove it by cropping
  # the reference template to this raster and stacking them.
  # A raster that stacks needs no resampling before it is
  # combined with other grid-aligned layers. This check uses
  # only geometry, so it runs even for all-NA rasters.
  if (check_grid_fit) {
    crs_ok <- terra::same.crs(r, ref_template)
    res_ok <- isTRUE(all.equal(
      res(r), ref_res_xy,
      tolerance = grid_fit_tolerance, scale = 1
    ))
    # Distance from the raster's lower-left corner to the
    # nearest reference grid line; aligned rasters sit on the
    # lattice, so both offsets are ~0.
    dx <- (xmin(r) - ref_anchor[1]) %% ref_res_xy[1]
    dy <- (ymin(r) - ref_anchor[2]) %% ref_res_xy[2]
    off_x <- min(dx, ref_res_xy[1] - dx)
    off_y <- min(dy, ref_res_xy[2] - dy)
    align_ok <- crs_ok &&
      off_x <= grid_fit_tolerance &&
      off_y <= grid_fit_tolerance

    # Definitive test: crop the reference to this raster's
    # extent, fill it with a dummy checkerboard so the 1 km
    # cells are visible, crop it to the AOI boundary, and
    # stack. c() errors unless the two geometries match
    # exactly, so a clean stack proves the raster fits and
    # yields a raster + grid stack for the grid-stack plot.
    stackable <- FALSE
    stack_msg <- ""
    if (crs_ok) {
      stackable <- tryCatch(
        {
          ref_crop <- crop(ref_template, ext(r))
          ref_dummy <- init(ref_crop, "chess")
          aoi_grid <- aoi_in_crs(ref_dummy)
          if (!is.null(aoi_grid)) {
            ref_dummy <- mask(ref_dummy, aoi_grid)
          }
          names(ref_dummy) <- "reference_grid"
          grid_stack <- c(r, ref_dummy)
          TRUE
        },
        error = function(e) {
          stack_msg <<- conditionMessage(e)
          FALSE
        }
      )
    } else {
      stack_msg <- "CRS differs from reference grid"
    }

    # Crop coverage: count valid pixels that fall outside the
    # AOI boundary. A raster cropped to the AOI has none.
    n_outside <- NA_integer_
    in_boundary <- NA
    aoi_r <- aoi_in_crs(r)
    if (!is.null(aoi_r)) {
      outside <- mask(r, aoi_r, inverse = TRUE)
      n_outside <- sum(!is.na(values(outside, mat = FALSE)))
      in_boundary <- n_outside == 0
    }

    message(
      "  Grid fit: crs ", if (crs_ok) "OK" else "MISMATCH",
      ", res ", if (res_ok) "OK" else "MISMATCH",
      ", aligned ", if (align_ok) "OK" else "OFF-GRID",
      ", stackable ", if (stackable) "YES" else "NO"
    )
    if (!stackable && nzchar(stack_msg)) {
      message("    Stack blocked: ", stack_msg)
    }
    if (!is.na(in_boundary)) {
      message(
        "  Crop to AOI: ",
        if (in_boundary) {
          "OK (no valid pixels outside boundary)"
        } else {
          paste0(n_outside, " valid pixels outside boundary")
        }
      )
    }

    grid_fit_list[[layer_name]] <- data.frame(
      raster = layer_name,
      crs_ok = crs_ok,
      res_ok = res_ok,
      aligned = align_ok,
      offset_x = off_x,
      offset_y = off_y,
      stackable = stackable,
      n_outside = n_outside,
      in_boundary = in_boundary
    )
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

  ## 2.2 Pixel statistics ----
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

  ## 2.3 Histogram ----
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

  ## 2.4 Coloured raster PNG ----
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
    if (overlay_boundary) {
      aoi_ov <- aoi_in_crs(r)
      if (!is.null(aoi_ov)) {
        terra::plot(aoi_ov, add = TRUE, border = "grey20", lwd = 1.2)
      }
    }
    grDevices::dev.off()
    message("  Wrote raster plot: ", raster_png_path)
  }

  ## 2.5 NA-mask PNG ----
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
    if (overlay_boundary) {
      aoi_ov <- aoi_in_crs(na_mask)
      if (!is.null(aoi_ov)) {
        terra::plot(aoi_ov, add = TRUE, border = "grey20", lwd = 1.2)
      }
    }
    grDevices::dev.off()
    message("  Wrote NA mask: ", na_mask_png_path)
  }

  ## 2.6 Coloured histogram PNG ----
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

  ## 2.7 Grid-stack visualization PNG ----
  # Plot the tested raster beside the dummy reference grid it
  # stacks with: a checkerboard of 1 km cells, cropped to the
  # AOI, sharing the raster's exact cells. Seeing the two
  # panels register confirms the raster sits on the grid.
  # Written only when the stack was built (needs check_grid_fit
  # and a raster that stacks).
  if (plot_grid_stack && !is.null(grid_stack)) {
    grid_stack_png_path <- file.path(
      out_dir,
      paste0(layer_name, "_grid_stack.png")
    )
    grDevices::png(
      filename = grid_stack_png_path,
      width = png_width,
      height = png_height,
      res = png_dpi
    )
    terra::plot(
      grid_stack,
      main = c(layer_name, "reference grid (dummy)"),
      axes = FALSE,
      box = FALSE,
      legend = TRUE
    )
    grDevices::dev.off()
    message("  Wrote grid stack: ", grid_stack_png_path)
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

# 4. Print grid-fit table ----
# Combine the per-raster grid-fit results into one table and
# print it. Only runs when check_grid_fit succeeded and at
# least one raster was tested. Any FALSE in crs_ok / res_ok /
# aligned / stackable flags a raster that needs resampling
# before it can be stacked on the reference grid; in_boundary
# FALSE flags a raster with valid pixels outside the AOI.
if (length(grid_fit_list) > 0) {
  grid_fit_df <- do.call(rbind, grid_fit_list)
  rownames(grid_fit_df) <- NULL
  message("Reference-grid fit:")
  print(grid_fit_df)
  if (all(grid_fit_df$stackable)) {
    message("All rasters stack on the reference grid.")
  } else {
    message(
      "WARNING: ",
      sum(!grid_fit_df$stackable),
      " raster(s) do NOT stack on the reference grid."
    )
  }
  n_out <- sum(!grid_fit_df$in_boundary, na.rm = TRUE)
  if (n_out > 0) {
    message(
      "WARNING: ",
      n_out,
      " raster(s) have valid pixels outside the AOI boundary."
    )
  }
}

# End of script ----
