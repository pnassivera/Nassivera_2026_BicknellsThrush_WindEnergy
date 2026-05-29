################################################################################
# Wind energy development alters breeding habitat of the Bicknell’s Thrush,
# a threatened range-restricted boreal songbird
#
# Nassivera et al. (2026)
#
# Objective: fully reproducible workflow for the paper titled above.
#            
# Includes:
# - North-sector pre-development reconstruction by IDW
# - injection of reconstructed values into the full Seigneurie-scale SDM
# - creation of pre-/post-development categorical habitat rasters
# - class-level landscape metrics before vs after development
# - park-by-park contribution decomposition
# - habitat patch-size distribution analysis
# - constraint-space modelling and spatial budget selection
# - iterative development scenarios and their landscape-metric impacts
#
# Last update: 2026-05-29
################################################################################

# ---- 0) LIBRARIES ------------------------------------------------------------

library(rstudioapi)
library(sf)
library(terra)
library(gstat)
library(sp)
library(landscapemetrics)

library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(tidytext)
library(ggplot2)
library(units)
library(stringr)

library(raster)
library(mapview)
library(leafsync)
library(patchwork)
library(mapedit)
library(viridisLite)
library(ggspatial)

library(plotly)
library(htmlwidgets)

# ---- 1) PATHS & SOME UTILITIES -----------------------------------------------

root_dir <- "C:/Users/Dell/Desktop/Bicknells_Thrush/Nassivera et al. 2026"

dir_data    <- file.path(root_dir, "02_Data")
dir_results <- file.path(root_dir, "03_Results")

# Output architecture
dir_ras <- file.path(dir_results, "Rasters")
dir_tab <- file.path(dir_results, "Tables")
dir_fig <- file.path(dir_results, "Figures")
dir_vec <- file.path(dir_results, "Vectors")
dir_int <- file.path(dir_results, "Interactive")

invisible(lapply(
  list(dir_results, dir_ras, dir_tab, dir_fig, dir_vec, dir_int),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

if (!dir.exists(dir_data)) {
  stop("Missing project data directory: ", dir_data)
}

# Geodatabase included in the reproducible project

gdb_path <- file.path(
  root_dir,
  "00_Géodatabase_GRBI",
  "GRBI_Projet_Thèse.gdb"
)

if (!dir.exists(gdb_path)) {
  stop("Missing geodatabase: ", gdb_path)
}

# Utilities

write_csv_bom <- function(df, path) {
  con <- file(path, open = "wb")
  writeBin(charToRaw("\xEF\xBB\xBF"), con)
  close(con)
  
  write.table(
    as.data.frame(df),
    file = path,
    sep = ",",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE,
    append = TRUE,
    fileEncoding = "UTF-8",
    qmethod = "double"
  )
}

calc_metrics_class <- function(r) {
  list(
    ca       = landscapemetrics::lsm_c_ca(r),
    np       = landscapemetrics::lsm_c_np(r),
    area_mn  = landscapemetrics::lsm_c_area_mn(r),
    lpi      = landscapemetrics::lsm_c_lpi(r),
    ed       = landscapemetrics::lsm_c_ed(r, count_boundary = TRUE),
    ai       = landscapemetrics::lsm_c_ai(r),
    tca      = landscapemetrics::lsm_c_tca(r, directions = 8, consider_boundary = FALSE, edge_depth = 3),
    cohesion = landscapemetrics::lsm_c_cohesion(r, directions = 8),
    mesh     = landscapemetrics::lsm_c_mesh(r),
    enn_mn   = landscapemetrics::lsm_c_enn_mn(r),
    area_cv  = landscapemetrics::lsm_c_area_cv(r),
    shape_cv = landscapemetrics::lsm_c_shape_cv(r, directions = 8)
  )
}

norm_lsm <- function(x, state_label) {
  x <- as.data.frame(x)
  x$state <- state_label
  
  keep <- intersect(names(x), c("layer", "level", "class", "metric", "value", "state"))
  x <- x[, keep, drop = FALSE]
  
  x$layer <- as.integer(x$layer)
  x$class <- as.integer(x$class)
  x$value <- as.numeric(x$value)
  x
}

# ---- 2) NORTH SECTOR REINTERPOLATION (IDW) -----------------------------------
# Goal: reconstruct pre-development occurrence probability within the existing
#       North-sector footprint using IDW (gstat), constrained to elevation >= 600 m.
# Notes:
# - "mitigation" and "protection" are the only thresholds used for reporting.
# - A lower internal threshold is used only to define reference points for IDW.

# 2.0 Inputs -------------------------------------------------------------------

prediction_path <- file.path(dir_data, "prediction_2024_beaupre.tif")
altitude_path   <- file.path(dir_data, "elevation_beaupre.tif")

layer_roads   <- "Chemins_20m_nord"
layer_turbine <- "Eoliennes_polygones_nord"

if (!file.exists(prediction_path)) {
  stop("Missing prediction raster: ", prediction_path)
}
if (!file.exists(altitude_path)) {
  stop("Missing elevation raster: ", altitude_path)
}

# 2.1 Parameters ---------------------------------------------------------------

# Reporting thresholds from Tremblay et al. (2025)
mitigation_threshold <- 0.103
protection_threshold <- 0.194

# Internal threshold for IDW reference points ONLY (eq. to average habitat quality)
reference_min_threshold <- 0.04

# Elevational constraint
elevation_min_m <- 600

# Footprint buffer (meters)
footprint_buffer_m <- 10

# LOOCV sample size
n_cv_points <- 1000

# IDW tuning
run_cv    <- TRUE   # TRUE = LOOCV; FALSE = use known best values
beta_best <- 1.85   # used only if run_cv == FALSE
n_best    <- 14     # used only if run_cv == FALSE

# 2.2 Load rasters -------------------------------------------------------------

r_pred <- terra::rast(prediction_path)
r_alt  <- terra::rast(altitude_path)

# 2.3 Build North footprint ----------------------------------------------------

merge_footprint <- function(gdb, layers) {
  sf_list <- lapply(layers, function(lyr) {
    sf::st_read(gdb, layer = lyr, quiet = TRUE) |>
      sf::st_geometry() |>
      sf::st_sf(geometry = _)
  })
  do.call(rbind, sf_list) |>
    sf::st_make_valid() |>
    sf::st_union() |>
    sf::st_sf(geometry = _)
}

footprint_sf <- merge_footprint(gdb_path, c(layer_roads, layer_turbine))

crs_pred_wkt <- terra::crs(r_pred, proj = TRUE)
if (is.na(sf::st_crs(footprint_sf))) stop("Missing CRS on footprint layers (GDB).")

footprint_sf <- sf::st_transform(footprint_sf, crs_pred_wkt)

footprint_buf_sf <- footprint_sf |>
  sf::st_make_valid() |>
  sf::st_buffer(dist = footprint_buffer_m) |>
  sf::st_make_valid() |>
  sf::st_union() |>
  sf::st_sf(geometry = _)

r_footprint <- terra::rasterize(
  x          = terra::vect(footprint_buf_sf),
  y          = r_pred,
  field      = 1,
  background = NA,
  touches    = TRUE
)

# 2.4 Define interpolation domain + reference points ---------------------------

# Restrict to the convex hull of the buffered footprint to reduce edge effects
aoi_sf   <- sf::st_sf(geometry = sf::st_convex_hull(sf::st_union(footprint_buf_sf)))
aoi_vect <- terra::vect(aoi_sf)

r_pred_aoi <- terra::mask(terra::crop(r_pred, aoi_vect), aoi_vect)
r_alt_aoi  <- terra::mask(terra::crop(r_alt,  aoi_vect), aoi_vect)
r_alt_aoi  <- terra::resample(r_alt_aoi, r_pred_aoi, method = "bilinear")

r_fp_aoi <- terra::mask(terra::crop(r_footprint, aoi_vect), aoi_vect)
fp_mask  <- !is.na(r_fp_aoi)

# Cells to fill: inside buffered footprint AND elevation >= threshold
fill_mask <- fp_mask & (r_alt_aoi >= elevation_min_m)
r_fill    <- terra::ifel(fill_mask, 1, NA)

# Masked raster (NA only where we will interpolate)
r_pred_masked <- r_pred_aoi
r_pred_masked[fill_mask] <- NA

# Reference points: outside footprint, elevation >= threshold, and >= internal minimum
ref_mask <- is.na(r_fp_aoi) &
  (r_alt_aoi >= elevation_min_m) &
  (r_pred_aoi >= reference_min_threshold)

r_pred_ref <- terra::mask(r_pred_aoi, ref_mask)

terra::writeRaster(r_fill,        file.path(dir_ras, "north_fill_domain.tif"), overwrite = TRUE)
terra::writeRaster(r_pred_masked, file.path(dir_ras, "north_pred_masked.tif"), overwrite = TRUE)

# 2.5 Utilities ----------------------------------------------------------------

rmse <- function(x) sqrt(mean(x^2, na.rm = TRUE))

sample_reference_points <- function(r_ref, n, min_thr, seed = 123) {
  set.seed(seed)
  r_ok <- terra::mask(r_ref, r_ref >= min_thr)
  
  pts <- terra::spatSample(
    x      = r_ok,
    size   = n,
    method = "random",
    na.rm  = TRUE,
    xy     = TRUE,
    values = TRUE
  )
  
  pts <- as.data.frame(pts)
  colnames(pts) <- c("x", "y", "value")
  pts <- pts[is.finite(pts$value) & !is.na(pts$value), , drop = FALSE]
  if (nrow(pts) < 50) stop("Too few reference points after filtering.")
  pts
}

spdf_from_xyv <- function(df, crs_wkt) {
  sp::SpatialPointsDataFrame(
    coords = df[, c("x", "y")],
    data   = data.frame(value = df$value),
    proj4string = sp::CRS(SRS_string = crs_wkt)
  )
}

cv_idw_loocv <- function(spatialDF, seq_neighbors, seq_beta, verbose = TRUE) {
  n <- nrow(spatialDF)
  if (n < 10) stop("Too few points for LOOCV.")
  
  seq_neighbors <- unique(seq_neighbors)
  seq_neighbors <- seq_neighbors[seq_neighbors <= (n - 1)]
  seq_neighbors <- unique(pmax(seq_neighbors, 3))
  
  seq_beta <- unique(pmax(seq_beta, 0.01))
  
  grid <- expand.grid(beta = seq_beta, neighbors = seq_neighbors)
  grid$rmse <- NA_real_
  
  for (i in seq_len(nrow(grid))) {
    idw_i <- gstat::gstat(
      formula = value ~ 1,
      data    = spatialDF,
      nmax    = grid$neighbors[i],
      set     = list(idp = grid$beta[i])
    )
    cv <- gstat::gstat.cv(idw_i, nfold = n)
    grid$rmse[i] <- rmse(cv$residual)
    
    if (isTRUE(verbose)) {
      message(sprintf(
        "CV %d/%d | beta=%.3f | n=%d | RMSE=%.6f",
        i, nrow(grid), grid$beta[i], grid$neighbors[i], grid$rmse[i]
      ))
    }
  }
  
  best_i <- which.min(grid$rmse)
  list(
    best_beta      = grid$beta[best_i],
    best_neighbors = grid$neighbors[best_i],
    best_rmse      = grid$rmse[best_i],
    grid           = grid
  )
}

# 2.6 LOOCV --------------------------------------------------------------------

# Runs LOOCV only if run_cv == TRUE. Otherwise uses beta_best / n_best from 2.1.

crs_wkt <- terra::crs(r_pred_aoi, proj = TRUE)

if (run_cv) {
  
  # Build CV sample
  df_cv <- sample_reference_points(
    r_ref   = r_pred_ref,
    n       = n_cv_points,
    min_thr = reference_min_threshold,
    seed    = 123
  )
  write.csv(df_cv, file.path(dir_tab, "cv_points_sample.csv"), row.names = FALSE)
  sp_cv <- spdf_from_xyv(df_cv, crs_wkt)
  
  # Coarse grid
  beta_gross <- c(0.1, seq(0.5, 3, 0.5))
  n_gross    <- unique(pmax(round(seq(3, nrow(sp_cv) - 1, length.out = 10)), 3))
  
  res_gross <- cv_idw_loocv(sp_cv, n_gross, beta_gross, verbose = TRUE)
  write.csv(res_gross$grid, file.path(dir_tab, "cv_grid_coarse.csv"), row.names = FALSE)
  
  best_gross <- data.frame(
    beta      = res_gross$best_beta,
    neighbors = res_gross$best_neighbors,
    rmse      = res_gross$best_rmse
  )
  write.csv(best_gross, file.path(dir_tab, "cv_best_coarse.csv"), row.names = FALSE)
  
  p_coarse <- ggplot(res_gross$grid, aes(x = neighbors, y = beta, fill = rmse)) +
    geom_tile(color = "black") +
    scale_fill_viridis_c() +
    theme_bw() +
    labs(
      title = "IDW LOOCV (coarse grid)",
      x = "Neighbors (nmax)",
      y = "Power (beta / idp)",
      fill = "RMSE"
    )
  ggsave(file.path(dir_fig, "IDW_LOOCV_coarse.png"), p_coarse, width = 8, height = 5, dpi = 300)
  
  # Fine grid around best
  b0 <- res_gross$best_beta
  n0 <- res_gross$best_neighbors
  nmax_allowed <- nrow(sp_cv) - 1
  
  beta_fine <- seq(max(0.05, b0 - 0.8), b0 + 0.8, by = 0.05)
  n_fine    <- seq(max(3, n0 - 100), min(nmax_allowed, n0 + 100), by = 10)
  
  res_fine <- cv_idw_loocv(sp_cv, n_fine, beta_fine, verbose = TRUE)
  write.csv(res_fine$grid, file.path(dir_tab, "cv_grid_fine.csv"), row.names = FALSE)
  
  best_fine <- data.frame(
    beta      = res_fine$best_beta,
    neighbors = res_fine$best_neighbors,
    rmse      = res_fine$best_rmse
  )
  write.csv(best_fine, file.path(dir_tab, "cv_best_fine.csv"), row.names = FALSE)
  
  best_fine$label <- sprintf(
    "Best\nβ = %.2f\nnmax = %d\nRMSE = %.4f",
    best_fine$beta,
    best_fine$neighbors,
    best_fine$rmse
  )
  
  p_fine <- ggplot(res_fine$grid, aes(x = neighbors, y = beta, fill = rmse)) +
    geom_tile(color = "black") +
    scale_fill_viridis_c() +
    geom_point(
      data = best_fine,
      aes(x = neighbors, y = beta),
      inherit.aes = FALSE,
      shape = 21, size = 3.5,
      fill = "white", color = "white"
    ) +
    geom_text(
      data = best_fine,
      aes(x = neighbors, y = beta, label = label),
      inherit.aes = FALSE,
      nudge_x = 8, nudge_y = 0.1,
      hjust = 0, vjust = 1,
      size = 4,
      color = "white"
    ) +
    theme_bw() +
    labs(
      title = "IDW LOOCV (fine grid)",
      x = "Neighbors (nmax)",
      y = "Power (beta / idp)",
      fill = "RMSE"
    )
  
  ggsave(file.path(dir_fig, "IDW_LOOCV_fine.png"), p_fine, width = 8, height = 5, dpi = 300)
  
  # Set best parameters for final interpolation
  beta_best <- res_fine$best_beta
  n_best    <- res_fine$best_neighbors
}

cat("\n--- IDW final parameters ---\n")

if (run_cv) {
  cat(sprintf(
    "beta (idp) = %.3f\nnmax = %d\nRMSE (LOOCV) = %.6f\n",
    beta_best, n_best, res_fine$best_rmse
  ))
} else {
  cat(sprintf("beta (idp) = %.3f\nnmax = %d\n", beta_best, n_best))
}

cat("----------------------------\n")

# 2.7 Final IDW interpolation --------------------------------------------------

cells_fill <- which(!is.na(terra::values(r_fill)) & terra::values(r_fill) == 1)
if (!length(cells_fill)) stop("No cells to fill (check footprint and elevation constraint).")

df_train <- as.data.frame(r_pred_ref, xy = TRUE, na.rm = TRUE)
colnames(df_train) <- c("x", "y", "value")

train_cells <- terra::cellFromXY(r_pred_aoi, df_train[, c("x", "y")])
df_train <- df_train[!(train_cells %in% cells_fill), , drop = FALSE]
if (nrow(df_train) < 50) stop("Too few training points after excluding fill cells.")

sp_train <- spdf_from_xyv(df_train, crs_wkt)

idw_model <- gstat::gstat(
  formula = value ~ 1,
  data    = sp_train,
  nmax    = n_best,
  set     = list(idp = beta_best)
)

xy_fill <- terra::xyFromCell(r_pred_aoi, cells_fill)
sp_new  <- sp::SpatialPoints(coords = xy_fill, proj4string = sp::CRS(SRS_string = crs_wkt))

pred <- predict(idw_model, newdata = sp_new, debug.level = 0)
pred_fill <- pred$var1.pred

r_predev_reconstructed <- r_pred_masked
r_predev_reconstructed[cells_fill] <- pred_fill

out_raster <- file.path(
  dir_ras,
  sprintf("predev_reconstructed_IDW_beta_%s_n_%s.tif", beta_best, n_best)
)
terra::writeRaster(r_predev_reconstructed, out_raster, overwrite = TRUE)

# 2.8 Summary tables -----------------------------------------------------------

px_ha <- (terra::res(r_predev_reconstructed)[1] * terra::res(r_predev_reconstructed)[2]) / 10000

v_fill <- terra::values(r_predev_reconstructed)[cells_fill]
v_fill <- v_fill[is.finite(v_fill) & !is.na(v_fill)]

summary_fill <- data.frame(
  n_cells = length(v_fill),
  area_ha = length(v_fill) * px_ha,
  mean    = mean(v_fill),
  median  = median(v_fill),
  sd      = sd(v_fill),
  p05     = as.numeric(stats::quantile(v_fill, 0.05, na.rm = TRUE)),
  p95     = as.numeric(stats::quantile(v_fill, 0.95, na.rm = TRUE)),
  min     = min(v_fill),
  max     = max(v_fill),
  beta    = beta_best,
  nmax    = n_best
)
write.csv(summary_fill, file.path(dir_tab, "reinterpolation_summary.csv"), row.names = FALSE)

v_class <- pmax(0, pmin(v_fill, 1))
v_class <- v_class[v_class >= mitigation_threshold]

breaks <- c(mitigation_threshold, protection_threshold, 1.0000001)
labels <- c("mitigation", "protection")

cls <- cut(v_class, breaks = breaks, labels = labels, right = FALSE, include.lowest = TRUE)
tab <- as.data.frame(table(cls), stringsAsFactors = FALSE)
colnames(tab) <- c("class", "n_pixels")
tab$area_ha <- tab$n_pixels * px_ha
tab$pct     <- 100 * tab$n_pixels / sum(tab$n_pixels)

write.csv(tab, file.path(dir_tab, "reinterpolation_class_composition.csv"), row.names = FALSE)

# 2.9 Figures -----------------------------------------------------------------

# Clamp interpolated values to [0,1]
v <- pmax(0, pmin(v_fill, 1))

# Cumulative distribution
ecdf_v <- stats::ecdf(v)
x_seq  <- seq(0, 1, length.out = 500)

cdf_df <- data.frame(
  suitability = x_seq,
  pct_ge      = 100 * (1 - ecdf_v(x_seq))
)

# Threshold palette
cols_thr <- c(
  Mitigation = "#e8ac31",
  Protection = "#a86466"
)

thr_df <- data.frame(
  threshold = c("Mitigation", "Protection"),
  x         = c(mitigation_threshold, protection_threshold),
  stringsAsFactors = FALSE
)

thr_df$pct_ge <- 100 * (1 - ecdf_v(thr_df$x))

# Add area (ha) to labels using the class table `tab`
area_map <- setNames(tab$area_ha, tab$class)

thr_df$class   <- tolower(thr_df$threshold)
thr_df$area_ha <- unname(area_map[thr_df$class])

thr_df$label <- sprintf(
  "≥ %.3f: %.1f%% (%.0f ha)",
  thr_df$x,
  thr_df$pct_ge,
  thr_df$area_ha
)

# Plot
p_cdf <- ggplot(cdf_df, aes(x = suitability, y = pct_ge)) +
  geom_line(linewidth = 1.1, color = "black") +
  geom_vline(
    xintercept = mitigation_threshold,
    linetype   = "dotted",
    linewidth  = 0.9,
    color      = cols_thr["Mitigation"]
  ) +
  geom_vline(
    xintercept = protection_threshold,
    linetype   = "dotted",
    linewidth  = 0.9,
    color      = cols_thr["Protection"]
  ) +
  geom_point(
    data = thr_df,
    aes(x = x, y = pct_ge, fill = threshold),
    shape = 21,
    size  = 3.5,
    color = "black"
  ) +
  geom_text(
    data = thr_df,
    aes(x = x, y = pct_ge, label = label),
    color = cols_thr[thr_df$threshold],
    hjust = -0.05,
    vjust = -0.6,
    size  = 4
  ) +
  scale_fill_manual(
    values = cols_thr,
    name   = NULL
  ) +
  theme_bw() +
  labs(
    title = "Cumulative profile of interpolated values",
    x     = "Occurrence probability (IDW)",
    y     = "Pixels ≥ X (%)"
  )

ggsave(
  file.path(dir_fig, "Reinterpolation_cumulative_profile.png"),
  p_cdf,
  width  = 7.5,
  height = 5,
  dpi    = 300
)

# 2.10 Zoom: before / after reinterpolation -----------------------------------

# Local objects
r_ref <- r_pred_aoi
r_int <- r_predev_reconstructed

# Convert to RasterLayer for mapview
rl_ref <- raster::raster(r_ref)
rl_int <- raster::raster(r_int)

# Shared visual contrast
zlim <- c(0, 0.6)
rl_ref_v <- raster::clamp(rl_ref, zlim[1], zlim[2])
rl_int_v <- raster::clamp(rl_int, zlim[1], zlim[2])

# Avoid pixel downsampling artefacts
mapviewOptions(mapview.maxpixels = max(ncell(rl_ref_v), ncell(rl_int_v)))

# Build mapview objects
mv_ref <- mapview(
  rl_ref_v,
  layer.name = "Post-development (before interpolation)",
  legend = TRUE
)

mv_int <- mapview(
  rl_int_v,
  layer.name = "Reconstructed pre-development (after interpolation)",
  legend = TRUE
)

# Extract exactly one leaflet widget per map
l_ref <- if (inherits(mv_ref@map, "leaflet")) mv_ref@map else mv_ref@map[[1]]
l_int <- if (inherits(mv_int@map, "leaflet")) mv_int@map else mv_int@map[[1]]

leafsync::sync(
  list(l_ref, l_int),
  ncol = 2
)

# 2.11 ROI selection and figure -----------------------------------------------

# Interactive ROI selection (draw rectangle on map)
drawn  <- mapedit::editMap(l_ref, editor = "leafpm")
roi_sf <- drawn$finished

# Ensure ROI has CRS (mapedit usually returns WGS84)
if (is.na(sf::st_crs(roi_sf))) {
  sf::st_crs(roi_sf) <- 4326
}

# Reproject ROI to raster CRS
roi_sf_r <- sf::st_transform(roi_sf, terra::crs(r_ref))
roi_v    <- terra::vect(roi_sf_r)

# Crop rasters to ROI
r_before_crop <- terra::crop(r_ref, roi_v)
r_before_crop <- terra::mask(r_before_crop, roi_v)

r_after_crop  <- terra::crop(r_int, roi_v)
r_after_crop  <- terra::mask(r_after_crop, roi_v)

# Convert to data frame
df_before <- as.data.frame(r_before_crop, xy = TRUE, na.rm = FALSE)
df_after  <- as.data.frame(r_after_crop,  xy = TRUE, na.rm = FALSE)

colnames(df_before)[3] <- "value"
colnames(df_after)[3]  <- "value"

# Side by side figure
zlim <- c(0, 0.6)

p_before <- ggplot(df_before, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c(limits = zlim, na.value = "white") +
  theme_void() +
  labs(title = "Post-development (before interpolation)")

p_after <- ggplot(df_after, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c(limits = zlim, na.value = "white") +
  theme_void() +
  labs(title = "Reconstructed pre-development (after interpolation)") +
  theme(legend.position = "none")

p_final <- p_before + p_after +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

p_final

ggsave(
  filename = file.path(dir_fig, "ROI_before_after_interpolation.png"),
  plot     = p_final,
  width    = 10,
  height   = 5,
  dpi      = 300
)

# ---- 3) INJECT NORTH REINTERPOLATION INTO THE FULL SDM -----------------------
# Goal: build a "pre-development" SDM at the Seigneurie scale by replacing
#       post-development values within the North-sector AOI with the reconstructed
#       values (IDW) computed in Step 2.

# Full pre-development raster initialized from the post-development SDM
r_predev_full <- r_pred
names(r_predev_full) <- "predev_SDM"

idx <- which(is.finite(terra::values(r_predev_reconstructed)))

if (!length(idx)) stop("No non-NA values found in r_predev_reconstructed.")

# Map AOI cells -> full raster cells using coordinates
xy_aoi     <- terra::xyFromCell(r_predev_reconstructed, idx)
cells_full <- terra::cellFromXY(r_predev_full, xy_aoi)

if (anyNA(cells_full)) {
  stop("Some reconstructed cells could not be matched to the full SDM grid.")
}

# Overwrite full raster at those locations
r_predev_full[cells_full] <- terra::values(r_predev_reconstructed)[idx]

# Quick sanity check: how many cells were replaced?
cat(sprintf("\nInjected %d cells into the full SDM.\n", length(cells_full)))

# Visual: ORIGINAL (post-dev) vs INJECTED PRE-DEV (North sector)

buffer_m <- 200

# Outline of the injected domain = buffered North footprint
fp_outline_sf <- sf::st_boundary(footprint_buf_sf)

# Buffer around the footprint for local crop (same window for both rasters)
aoi_buf_sf <- sf::st_buffer(footprint_buf_sf, dist = buffer_m)
aoi_buf_v  <- terra::vect(aoi_buf_sf)

# Local crop (same extent)
r_post_crop <- terra::crop(r_pred,        aoi_buf_v)
r_pre_crop  <- terra::crop(r_predev_full, aoi_buf_v)

# Shared clamp for visualization
zlim <- c(0, 0.6)
r_post_viz <- terra::clamp(r_post_crop, zlim[1], zlim[2])
r_pre_viz  <- terra::clamp(r_pre_crop,  zlim[1], zlim[2])

# Mask + trim (remove bbox frame)
r_post_viz <- terra::trim(terra::mask(r_post_viz, aoi_buf_v))
r_pre_viz  <- terra::trim(terra::mask(r_pre_viz,  aoi_buf_v))

# Mapview (synced) + footprint outline overlay
mapviewOptions(mapview.maxpixels = max(
  terra::ncell(r_post_viz),
  terra::ncell(r_pre_viz)
))

mv_post <- mapview(
  r_post_viz,
  layer.name = "Post-development SDM (original) — North +200 m",
  legend     = TRUE,
  na.color   = "transparent"
) + mapview(fp_outline_sf, layer.name = "North footprint", color = "black", lwd = 2)

mv_pre <- mapview(
  r_pre_viz,
  layer.name = "Pre-development SDM (injected) — North +200 m",
  legend     = TRUE,
  na.color   = "transparent"
) + mapview(fp_outline_sf, layer.name = "North footprint", color = "black", lwd = 2)

# Extract exactly one leaflet widget per map
l_post <- if (inherits(mv_post@map, "leaflet")) mv_post@map else mv_post@map[[1]]
l_pre  <- if (inherits(mv_pre@map,  "leaflet")) mv_pre@map  else mv_pre@map[[1]]

leafsync::sync(list(l_post, l_pre), ncol = 2)

# Save
out_predev_full <- file.path(
  dir_ras,
  sprintf("SDM_predev_injected_North_IDW_beta_%s_n_%s.tif", beta_best, n_best)
)

terra::writeRaster(r_predev_full, out_predev_full, overwrite = TRUE)

# ---- 4) LANDSCAPEMETRICS INPUT RASTERS (PRE/POST DEVELOPMENT) ----------------
# Goal: produce integer categorical habitat rasters (0–2) for landscapemetrics:
#       (i) pre-development (North injected) and (ii) post-development.
#       Then enforce the wind footprint as class 0 on the post-development raster.

# Classes:
# 0 = non_habitat      : < mitigation_threshold
# 1 = mitigation_class : [mitigation_threshold, protection_threshold)
# 2 = protection_class : [protection_threshold, 1]

# Wind farm footprints (all sectors) + buffered footprint

# Public-data reconstruction assumptions:
# - turbine pad area = 1.2 ha
# - road corridor total width = 40 m
# - final additional buffer = footprint_buffer_m

# Target surface per turbine: 1.2 ha = 12,000 m²
turbine_pad_area_m2 <- 12000

# Equivalent disk radius: A = pi * r²
turbine_pad_radius_m <- sqrt(turbine_pad_area_m2 / pi)

# Desired total road width = 40 m
# therefore buffer of 20 m on each side
road_buffer_m <- 20

# North layers

layers_north <- c(
  "Chemins_20m_nord",
  "Eoliennes_polygones_nord"
)

north_list <- lapply(layers_north, function(lyr) {
  x <- sf::st_read(gdb_path, layer = lyr, quiet = TRUE)
  if (nrow(x) == 0) return(NULL)
  x
})
north_list <- Filter(Negate(is.null), north_list)

# Public layers reconstructed for East / West / South

layers_points <- c(
  "Eoliennes_points_est",
  "Eoliennes_points_ouest",
  "Eoliennes_points_sud"
)

layers_lines <- c(
  "Chemins_lignes_est",
  "Chemins_lignes_ouest",
  "Chemins_lignes_sud"
)

pts_list <- lapply(layers_points, function(lyr) {
  x <- sf::st_read(gdb_path, layer = lyr, quiet = TRUE)
  if (nrow(x) == 0) return(NULL)
  x
})
pts_list <- Filter(Negate(is.null), pts_list)

lines_list <- lapply(layers_lines, function(lyr) {
  x <- sf::st_read(gdb_path, layer = lyr, quiet = TRUE)
  if (nrow(x) == 0) return(NULL)
  x
})
lines_list <- Filter(Negate(is.null), lines_list)

# Shared preparation of all input layers

prep_sf <- function(x) {
  x <- sf::st_transform(x, terra::crs(r_pred))
  x <- sf::st_zm(x, drop = TRUE, what = "ZM")
  x <- sf::st_make_valid(x)
  x
}

north_list <- lapply(north_list, prep_sf)
pts_list   <- lapply(pts_list,   prep_sf)
lines_list <- lapply(lines_list, prep_sf)

if (!length(north_list) && !length(pts_list) && !length(lines_list)) {
  stop("No wind-footprint features were found in the geodatabase.")
}

# Keep North geometries as they are

geom_north <- if (length(north_list)) {
  do.call(
    c,
    lapply(north_list, function(x) sf::st_geometry(x))
  )
} else {
  NULL
}

# Reconstruct turbine footprints (1.2 ha disks)

geom_pts <- if (length(pts_list)) {
  do.call(
    c,
    lapply(pts_list, function(x) {
      sf::st_buffer(sf::st_geometry(x), dist = turbine_pad_radius_m)
    })
  )
} else {
  NULL
}

# Reconstruct road corridors (40 m total width)

geom_lines <- if (length(lines_list)) {
  do.call(
    c,
    lapply(lines_list, function(x) {
      sf::st_buffer(sf::st_geometry(x), dist = road_buffer_m)
    })
  )
} else {
  NULL
}

# Merge North + reconstructed East/West/South

geom_all <- c(geom_north, geom_pts, geom_lines)

footprint_all_buf <- sf::st_union(geom_all)
footprint_all_buf <- sf::st_collection_extract(footprint_all_buf, "POLYGON")

# Apply final buffer

footprint_all_buf <- sf::st_buffer(footprint_all_buf, dist = footprint_buffer_m)
footprint_all_buf <- sf::st_make_valid(footprint_all_buf)
footprint_all_buf <- sf::st_collection_extract(footprint_all_buf, "POLYGON")

footprint_all_sf <- sf::st_sf(geometry = footprint_all_buf)
footprint_all_v  <- terra::vect(footprint_all_sf)

out_fp_gpkg <- file.path(
  dir_vec,
  sprintf("wind_footprint_total_reconstructed_buffer%sm.gpkg", footprint_buffer_m)
)

sf::st_write(footprint_all_sf, out_fp_gpkg, delete_dsn = TRUE, quiet = TRUE)

# Waterbodies mask 

water_mask_path <- file.path(dir_data, "water_mask_SeigneurieBeaupre_COTQ2022.tif")

r_water <- terra::rast(water_mask_path)

# Force exact CRS match
if (!terra::same.crs(r_water, r_predev_full)) {
  r_water <- terra::project(r_water, r_predev_full, method = "near")
}

# Then force exact grid match
if (!terra::compareGeom(r_water, r_predev_full, stopOnError = FALSE)) {
  r_water <- terra::resample(r_water, r_predev_full, method = "near")
}

# Masking
r_predev_full[!is.na(r_water)] <- NA
r_pred[!is.na(r_water)]        <- NA

# Quick visual check of water mask

cols_water <- c("#2b83ba")

mv_water <- mapview(
  raster::raster(r_water),
  layer.name    = "Water mask",
  col.regions   = cols_water,
  alpha.regions = 0.7,
  na.color      = "transparent",
  legend        = TRUE
)

mv_pre_watercheck <- mapview(
  raster::raster(terra::clamp(r_predev_full, 0, 0.6)),
  layer.name = "Pre-development SDM masked",
  na.color   = "transparent",
  legend     = TRUE
) + mv_water

mv_post_watercheck <- mapview(
  raster::raster(terra::clamp(r_pred, 0, 0.6)),
  layer.name = "Post-development SDM masked",
  na.color   = "transparent",
  legend     = TRUE
) + mv_water

l1 <- if (inherits(mv_pre_watercheck@map, "leaflet")) mv_pre_watercheck@map else mv_pre_watercheck@map[[1]]
l2 <- if (inherits(mv_post_watercheck@map, "leaflet")) mv_post_watercheck@map else mv_post_watercheck@map[[1]]

leafsync::sync(list(l1, l2), ncol = 2)

cat(
  "Non-NA habitat pixels still present inside water mask (r_predev_full):",
  terra::global(!is.na(r_predev_full) & !is.na(r_water), "sum", na.rm = TRUE)[1,1],
  "\n"
)

cat(
  "Non-NA habitat pixels still present inside water mask (r_pred):",
  terra::global(!is.na(r_pred) & !is.na(r_water), "sum", na.rm = TRUE)[1,1],
  "\n"
)

# Classification (2 thresholds)

thr_mit <- mitigation_threshold
thr_con <- protection_threshold

rcl <- matrix(c(
  -Inf,    thr_mit, 0,
  thr_mit, thr_con, 1,
  thr_con, Inf,     2
), ncol = 3, byrow = TRUE)

r_class_pre <- terra::classify(r_predev_full, rcl, include.lowest = TRUE, right = FALSE)
r_class_pre <- terra::as.int(r_class_pre)
names(r_class_pre) <- "habitat_class_predev"

r_class_post <- terra::classify(r_pred, rcl, include.lowest = TRUE, right = FALSE)
r_class_post <- terra::as.int(r_class_post)
names(r_class_post) <- "habitat_class_postdev"

lvl_df_raster <- data.frame(
  value = 0:2,
  class = c("non_habitat", "mitigation_class", "protection_class")
)

levels(r_class_pre)  <- list(lvl_df_raster)
levels(r_class_post) <- list(lvl_df_raster)

# Rasterize the full wind footprint explicitly on the habitat grid so that
# footprint application is fully consistent with the sector-by-sector logic
# used elsewhere in the workflow.

r_footprint_all <- terra::rasterize(
  x          = footprint_all_v,
  y          = r_class_post,
  field      = 1,
  background = NA,
  touches    = TRUE
)

r_footprint_all <- terra::ifel(r_footprint_all == 1, 1L, NA)
names(r_footprint_all) <- "wind_footprint_all"

# Apply footprint to post-development class raster (force class 0)
r_class_post_imp <- r_class_post
r_class_post_imp[!is.na(r_footprint_all)] <- 0L
r_class_post_imp <- terra::as.int(r_class_post_imp)
names(r_class_post_imp) <- "habitat_class_postdev_impacted"
levels(r_class_post_imp) <- list(lvl_df_raster)

print(landscapemetrics::check_landscape(r_class_pre))
print(landscapemetrics::check_landscape(r_class_post_imp))

out_pre_tif  <- file.path(dir_ras, "habitat_class_2thresholds_predev.tif")
out_post_tif <- file.path(dir_ras, "habitat_class_2thresholds_postdev_impacted.tif")
out_fp_tif   <- file.path(
  dir_ras,
  sprintf("wind_footprint_total_reconstructed_buffer%sm_raster.tif", footprint_buffer_m)
)

out_predev_full_masked_tif <- file.path(dir_ras, "SDM_predev_injected_North_masked_water.tif")
out_pred_masked_tif        <- file.path(dir_ras, "prediction_2024_beaupre_masked_water.tif")
out_water_tif              <- file.path(dir_ras, "water_mask_SeigneurieBeaupre_COTQ2022.tif")

terra::writeRaster(r_class_pre,      out_pre_tif,               overwrite = TRUE)
terra::writeRaster(r_class_post_imp, out_post_tif,              overwrite = TRUE)
terra::writeRaster(r_footprint_all,  out_fp_tif,                overwrite = TRUE)
terra::writeRaster(r_predev_full,    out_predev_full_masked_tif, overwrite = TRUE)
terra::writeRaster(r_pred,           out_pred_masked_tif,        overwrite = TRUE)
terra::writeRaster(r_water,          out_water_tif,              overwrite = TRUE)

# Mapview (pre-dev vs post-dev impacted)

cols_2thr_map <- c(
  "#BEBEBE",
  "#e8ac31",
  "#a86466"
)

r_pre_viz  <- r_class_pre
r_post_viz <- r_class_post_imp

max_cells <- 2e7
if (terra::ncell(r_pre_viz) > max_cells) {
  fact <- ceiling(sqrt(terra::ncell(r_pre_viz) / max_cells))
  r_pre_viz  <- terra::aggregate(r_pre_viz,  fact = fact, fun = "modal", na.rm = TRUE)
  r_post_viz <- terra::aggregate(r_post_viz, fact = fact, fun = "modal", na.rm = TRUE)
}

rl_pre  <- raster::ratify(raster::raster(r_pre_viz))
rl_post <- raster::ratify(raster::raster(r_post_viz))

rat_pre <- raster::levels(rl_pre)[[1]]
rat_pre$class <- c("Non habitat", "Mitigation", "Protection")
levels(rl_pre) <- rat_pre

rat_post <- raster::levels(rl_post)[[1]]
rat_post$class <- c("Non habitat", "Mitigation", "Protection")
levels(rl_post) <- rat_post

mapviewOptions(mapview.maxpixels = max(raster::ncell(rl_pre), raster::ncell(rl_post)))

mv_pre <- mapview(
  rl_pre,
  layer.name  = "Pre-development — classified (2 thresholds)",
  col.regions = cols_2thr_map,
  na.color    = "transparent",
  legend      = TRUE
)

mv_post <- mapview(
  rl_post,
  layer.name  = "Post-development — classified (footprint forced to 0)",
  col.regions = cols_2thr_map,
  na.color    = "transparent",
  legend      = TRUE
)

l_pre  <- if (inherits(mv_pre@map,  "leaflet")) mv_pre@map  else mv_pre@map[[1]]
l_post <- if (inherits(mv_post@map, "leaflet")) mv_post@map else mv_post@map[[1]]

leafsync::sync(list(l_pre, l_post), ncol = 2)

# Exploratory class statistics
# Quantify the area (ha) and percentage of the Seigneurie occupied by the
# mitigation and protection habitat classes under pre-development conditions.

calc_class_stats_seigneurie <- function(r_cls, total_ha) {
  
  pix_ha <- prod(terra::res(r_cls)) / 10000
  
  n_mit <- terra::global(r_cls == 1, "sum", na.rm = TRUE)[1, 1]
  n_con <- terra::global(r_cls == 2, "sum", na.rm = TRUE)[1, 1]
  
  out <- tibble::tibble(
    class = c(1L, 2L),
    class_label = c("Mitigation", "Protection"),
    n_cells = c(n_mit, n_con),
    area_ha = n_cells * pix_ha,
    pct_seigneurie = 100 * area_ha / total_ha
  )
  
  return(out)
}

seigneurie_total_ha <- 159221

stats_pre <- calc_class_stats_seigneurie(
  r_cls    = r_class_pre,
  total_ha = seigneurie_total_ha
)

print(stats_pre, n = Inf)

# ---- 5) LANDSCAPE METRICS (CLASS LEVEL) --------------------------------------
# Goal: quantify how wind farm development altered habitat composition and
#       configuration by comparing pre-development vs post-development class maps

r_before <- r_class_pre
names(r_before) <- "habitat_predev"

r_after <- r_class_post_imp
names(r_after) <- "habitat_postdev"

print(landscapemetrics::check_landscape(r_before))
print(landscapemetrics::check_landscape(r_after))

# 5.1 Class-level metrics ------------------------------------------------------

metrics_before <- calc_metrics_class(r_before)
metrics_after  <- calc_metrics_class(r_after)

# 5.2 Dictionaries: category / units / labels ---------------------------------

metric_dict <- data.frame(
  metric = c("ca","np","area_mn","lpi","ed","ai","tca","cohesion","mesh","enn_mn","area_cv","shape_cv"),
  category = c(
    "Amount",
    "Fragmentation",
    "Fragmentation",
    "Fragmentation",
    "Fragmentation",
    "Fragmentation",
    "Fragmentation",
    "Connectivity",
    "Connectivity",
    "Connectivity",
    "Heterogeneity",
    "Heterogeneity"
  ),
  stringsAsFactors = FALSE
)

metric_units <- data.frame(
  metric = c("ca","np","area_mn","lpi","ed","ai","tca","cohesion","mesh","enn_mn","area_cv","shape_cv"),
  unit = c("ha","count","ha","%","m/ha","%","ha","%","ha","m","-","-"),
  stringsAsFactors = FALSE
)

metric_labels <- c(
  ca        = "CA (ha)",
  np        = "NP (count)",
  area_mn   = "AREA_MN (ha)",
  lpi       = "LPI (%)",
  ed        = "ED (m/ha)",
  ai        = "AI (%)",
  tca       = "TCA (ha)",
  cohesion  = "COHESION (%)",
  mesh      = "MESH (ha)",
  enn_mn    = "ENN_MN (m)",
  area_cv   = "AREA_CV (-)",
  shape_cv  = "SHAPE_CV (-)"
)

metric_labels_df <- tibble(
  metric = names(metric_labels),
  metric_label = unname(metric_labels)
)

# 5.3 Normalize landscapemetrics outputs --------------------------------------

before_df <- do.call(rbind, lapply(metrics_before, norm_lsm, state_label = "before"))
after_df  <- do.call(rbind, lapply(metrics_after,  norm_lsm, state_label = "after"))
all_df    <- rbind(before_df, after_df)

all_df <- merge(all_df, metric_dict,  by = "metric", all.x = TRUE)
all_df <- merge(all_df, metric_units, by = "metric", all.x = TRUE)

if (any(is.na(all_df$category))) {
  stop("Some metrics could not be matched to a category.")
}
if (any(is.na(all_df$unit))) {
  stop("Some metrics could not be matched to a unit.")
}

all_df <- all_df %>% dplyr::filter(is.na(class) | class != 0)

# 5.4 Wide table + delta -------------------------------------------------------

wide <- reshape(
  all_df[, c("category","metric","unit","level","class","state","value")],
  timevar   = "state",
  idvar     = c("category","metric","unit","level","class"),
  direction = "wide"
)

wide <- wide[, c("category","metric","unit","level","class","value.before","value.after")]
wide$delta <- wide$value.after - wide$value.before

cat_order <- c("Amount","Fragmentation","Connectivity","Heterogeneity")
wide$category <- factor(wide$category, levels = cat_order, ordered = TRUE)
wide <- wide[order(wide$category, wide$metric, wide$class), ]

wide_tbl <- as_tibble(wide)

options(scipen = 999, pillar.sigfig = 4)
print(wide_tbl, n = Inf)

# 5.5 Export CSV ---------------------------------------------------------------

out_csv <- file.path(dir_tab, "table_all_metrics_2thr.csv")
write_csv_bom(wide_tbl, out_csv)

# 5.6 Forest plot (Δ% per metric × class) -------------------------------------

pal_class <- c(
  "Mitigation" = "#e8ac31",
  "Protection" = "#a86466"
)

forest_df <- wide_tbl %>%
  dplyr::filter(level == "class", class %in% c(1, 2)) %>%
  dplyr::left_join(metric_labels_df, by = "metric") %>%
  dplyr::mutate(
    class = factor(
      as.character(class),
      levels = c("1", "2"),
      labels = c("Mitigation", "Protection")
    ),
    delta_pct = dplyr::if_else(
      is.finite(value.before) & !is.na(value.before) & value.before > 0,
      (value.after - value.before) / value.before * 100,
      NA_real_
    ),
    out_5 = !is.na(delta_pct) & abs(delta_pct) > 5,
    force = abs(delta_pct)
  )

metric_order <- forest_df %>%
  dplyr::filter(!is.na(delta_pct)) %>%
  dplyr::group_by(category, metric_label) %>%
  dplyr::summarise(score_ord = max(abs(delta_pct), na.rm = TRUE), .groups = "drop")

forest_df <- forest_df %>%
  dplyr::left_join(metric_order, by = c("category", "metric_label")) %>%
  dplyr::mutate(
    score_ord = dplyr::coalesce(score_ord, -Inf),
    metric_lab = tidytext::reorder_within(metric_label, -score_ord, category)
  )

p_forest <- ggplot(
  forest_df,
  aes(x = delta_pct, y = metric_lab, color = class)
) +
  annotate(
    "rect",
    xmin = -5, xmax = 5,
    ymin = -Inf, ymax = Inf,
    alpha = 0.30,
    fill  = "grey60"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_point(
    data = dplyr::filter(forest_df, !is.na(delta_pct)),
    aes(alpha = out_5, size = force)
  ) +
  scale_alpha_manual(values = c(`FALSE` = 0.65, `TRUE` = 1), guide = "none") +
  scale_size_continuous(
    range = c(2.5, 10),
    name  = "Change magnitude (%)",
    guide = guide_legend(
      override.aes = list(color = "grey30", alpha = 1, size = 6)
    )
  ) +
  facet_wrap(~ category, scales = "free_y", ncol = 1) +
  tidytext::scale_y_reordered() +
  scale_color_manual(
    values = pal_class,
    name   = "Habitat class",
    drop   = FALSE,
    guide  = guide_legend(
      override.aes = list(size = 6, alpha = 1)
    )
  ) +
  labs(
    x = "Δ (%)",
    y = NULL
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4)
  )

print(p_forest)

ggsave(
  filename = file.path(dir_fig, "forest_plot_metrics_2thr.png"),
  plot     = p_forest,
  width    = 44,
  height   = 32,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

out_csv_forest <- file.path(dir_tab, "table_forest_plot_metrics_2thr.csv")
write_csv_bom(forest_df, out_csv_forest)

# ---- 6) PARK-BY-PARK MAIN EFFECTS (CONTRIBUTION DECOMPOSITION) ---------------

sector_specs <- list(
  North = list(
    mode = "raw",
    raw_layers = c("Chemins_20m_nord", "Eoliennes_polygones_nord")
  ),
  West = list(
    mode = "reconstruct",
    point_layer = "Eoliennes_points_ouest",
    line_layer  = "Chemins_lignes_ouest"
  ),
  South = list(
    mode = "reconstruct",
    point_layer = "Eoliennes_points_sud",
    line_layer  = "Chemins_lignes_sud"
  ),
  East = list(
    mode = "reconstruct",
    point_layer = "Eoliennes_points_est",
    line_layer  = "Chemins_lignes_est"
  )
)

read_layer_sf <- function(gdb, layer_name) {
  x <- sf::st_read(gdb, layer = layer_name, quiet = TRUE)
  if (nrow(x) == 0) return(NULL)
  x <- sf::st_zm(x, drop = TRUE, what = "ZM")
  x <- sf::st_make_valid(x)
  x
}

prep_to_raster_crs <- function(x, r_ref) {
  x <- sf::st_transform(x, terra::crs(r_ref))
  x <- sf::st_zm(x, drop = TRUE, what = "ZM")
  x <- sf::st_make_valid(x)
  x
}

build_sector_footprint_raster <- function(gdb, spec, r_ref,
                                          footprint_buffer_m,
                                          turbine_pad_area_m2 = 12000,
                                          road_total_width_m = 40) {
  
  turbine_pad_radius_m <- sqrt(turbine_pad_area_m2 / pi)
  road_buffer_m <- road_total_width_m / 2
  
  if (identical(spec$mode, "raw")) {
    
    sf_list <- lapply(spec$raw_layers, function(lyr) read_layer_sf(gdb, lyr))
    sf_list <- Filter(Negate(is.null), sf_list)
    
    if (!length(sf_list)) {
      stop("No valid raw features were found for this sector.")
    }
    
    sf_list <- lapply(sf_list, prep_to_raster_crs, r_ref = r_ref)
    geom_all <- do.call(c, lapply(sf_list, sf::st_geometry))
    
  } else if (identical(spec$mode, "reconstruct")) {
    
    pts <- read_layer_sf(gdb, spec$point_layer)
    lns <- read_layer_sf(gdb, spec$line_layer)
    
    geom_pts <- NULL
    geom_lns <- NULL
    
    if (!is.null(pts)) {
      pts <- prep_to_raster_crs(pts, r_ref = r_ref)
      geom_pts <- sf::st_buffer(sf::st_geometry(pts), dist = turbine_pad_radius_m)
    }
    
    if (!is.null(lns)) {
      lns <- prep_to_raster_crs(lns, r_ref = r_ref)
      geom_lns <- sf::st_buffer(sf::st_geometry(lns), dist = road_buffer_m)
    }
    
    geom_all <- c(geom_pts, geom_lns)
    
    if (!length(geom_all)) {
      stop("No reconstructed features were found for this sector.")
    }
    
  } else {
    stop("Unknown sector mode.")
  }
  
  geom_u <- sf::st_union(geom_all)
  geom_u <- sf::st_make_valid(geom_u)
  geom_u <- sf::st_collection_extract(geom_u, "POLYGON")
  
  geom_u <- sf::st_buffer(geom_u, dist = footprint_buffer_m)
  geom_u <- sf::st_make_valid(geom_u)
  geom_u <- sf::st_collection_extract(geom_u, "POLYGON")
  
  fp_sf <- sf::st_sf(geometry = geom_u)
  fp_v  <- terra::vect(fp_sf)
  
  fp_r <- terra::rasterize(fp_v, r_ref, field = 1, background = NA, touches = TRUE)
  terra::ifel(fp_r == 1, 1, NA)
}

apply_footprint_as_removed <- function(r0, fp_r) {
  r <- r0
  r[!is.na(fp_r)] <- 0L
  r
}

fp_rasters <- lapply(names(sector_specs), function(nm) {
  build_sector_footprint_raster(
    gdb                 = gdb_path,
    spec                = sector_specs[[nm]],
    r_ref               = r_before,
    footprint_buffer_m  = footprint_buffer_m,
    turbine_pad_area_m2 = 12000,
    road_total_width_m  = 40
  )
})
names(fp_rasters) <- names(sector_specs)

states_main <- c(
  list(none = r_before),
  lapply(names(fp_rasters), function(nm) {
    apply_footprint_as_removed(r_before, fp_rasters[[nm]])
  })
)
names(states_main) <- c("none", names(fp_rasters))

metrics_by_state <- lapply(names(states_main), function(st) {
  met <- calc_metrics_class(states_main[[st]])
  do.call(rbind, lapply(met, norm_lsm, state_label = st))
})
names(metrics_by_state) <- names(states_main)

states_df <- do.call(rbind, metrics_by_state)

states_df <- merge(states_df, metric_dict,  by = "metric", all.x = TRUE)
states_df <- merge(states_df, metric_units, by = "metric", all.x = TRUE)

if (any(is.na(states_df$category))) {
  stop("Some metrics could not be matched to a category.")
}
if (any(is.na(states_df$unit))) {
  stop("Some metrics could not be matched to a unit.")
}

states_df <- states_df %>% filter(is.na(class) | class != 0)

wide_states <- reshape(
  states_df[, c("category","metric","unit","level","class","state","value")],
  timevar   = "state",
  idvar     = c("category","metric","unit","level","class"),
  direction = "wide"
)

wide_base <- wide_states %>%
  filter(level == "class", class %in% c(1, 2))

impact_one <- function(sector_name) {
  before <- wide_base$value.none
  after  <- wide_base[[paste0("value.", sector_name)]]
  
  ok <- is.finite(before) &
    !is.na(before) &
    is.finite(after) &
    !is.na(after) &
    before > 0
  
  if (!any(ok)) {
    return(NA_real_)
  }
  
  sum(abs((after[ok] - before[ok]) / before[ok] * 100), na.rm = TRUE)
}

impact_tbl <- tibble(
  sector = names(fp_rasters),
  impact = sapply(names(fp_rasters), impact_one)
) %>%
  mutate(
    pct = if (all(is.na(impact))) NA_real_ else impact / sum(impact, na.rm = TRUE) * 100,
    sector = factor(sector, levels = c("North", "West", "South", "East"))
  )

pal_sector <- c(
  "North" = "#4E79A7",
  "West"  = "#F28E2B",
  "South" = "#59A14F",
  "East"  = "#E15759"
)

p_stack <- ggplot(impact_tbl, aes(x = "Total", y = pct, fill = sector)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = paste0(round(pct, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 6
  ) +
  scale_fill_manual(values = pal_sector, name = "Sector") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Relative contribution (%)",
    title = "Relative contributions of each sector"
  ) +
  theme_bw(base_size = 22) +
  theme(
    legend.position    = "top",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p_stack)

ggsave(
  filename = file.path(dir_fig, "effects_stacked_barplot.png"),
  plot     = p_stack,
  width    = 26,
  height   = 24,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

write_csv_bom(
  impact_tbl,
  file.path(dir_tab, "main_effects_contributions.csv")
)

# ---- 7) HABITAT PATCH-SIZE DISTRIBUTIONS (AREA-WEIGHTED) ---------------------

# Goal: characterize patch-size structure (mitigation vs protection) under
#       pre- vs post-development, using area-weighted mean patch size: μw = Σ a² / Σ a

class_labs <- c(
  "1" = "Mitigation",
  "2" = "Protection"
)

pal_thr <- c(
  "Mitigation" = "#e8ac31",
  "Protection" = "#a86466"
)

time_labs <- c(
  before = "Pre-development",
  after  = "Post-development"
)

extract_patch_area_tbl <- function(r_cat, time_label, class_labs) {
  a <- lsm_p_area(r_cat)
  
  a %>%
    filter(level == "patch") %>%
    transmute(
      time_raw  = time_label,
      hab_class = factor(
        as.character(class),
        levels = names(class_labs),
        labels = unname(class_labs)
      ),
      patch_id  = as.integer(id),
      area_ha   = as.numeric(value)
    ) %>%
    filter(!is.na(hab_class), is.finite(area_ha), area_ha > 0)
}

patch_tbl <- bind_rows(
  extract_patch_area_tbl(r_before, "before", class_labs),
  extract_patch_area_tbl(r_after,  "after",  class_labs)
) %>%
  mutate(
    time_raw = factor(time_raw, levels = c("before", "after")),
    time     = factor(time_raw, levels = c("before", "after"), labels = unname(time_labs))
  )

write_csv_bom(patch_tbl, file.path(dir_tab, "patch_areas_by_class_time.csv"))

mu_w <- function(a) {
  a <- a[is.finite(a) & a > 0]
  if (!length(a)) return(NA_real_)
  sum(a^2) / sum(a)
}

mu_tbl_long <- patch_tbl %>%
  group_by(hab_class, time_raw) %>%
  summarise(
    n_patches  = n(),
    total_area = sum(area_ha),
    mu_w_ha    = mu_w(area_ha),
    .groups    = "drop"
  )

mu_tbl <- mu_tbl_long %>%
  pivot_wider(
    names_from  = time_raw,
    values_from = c(n_patches, total_area, mu_w_ha),
    names_sep   = "_"
  ) %>%
  mutate(
    delta_mu_w = mu_w_ha_after - mu_w_ha_before,
    pct_mu_w   = 100 * (mu_w_ha_after / mu_w_ha_before - 1)
  )

write_csv_bom(mu_tbl, file.path(dir_tab, "mu_w_by_class_before_after.csv"))

x_map <- c("Pre-development" = 1, "Post-development" = 2)

build_weighted_violin <- function(df, width = 0.38, n = 700) {
  df %>%
    filter(is.finite(area_ha), area_ha > 0, !is.na(time), !is.na(hab_class)) %>%
    mutate(
      x0   = unname(x_map[as.character(time)]),
      ylog = log10(area_ha),
      w    = area_ha
    ) %>%
    group_by(hab_class, time) %>%
    group_map(~{
      ylog <- .x$ylog
      w    <- .x$w
      x0   <- .x$x0[1]
      
      d <- stats::density(ylog, weights = w / sum(w), n = n)
      
      dens <- d$y
      dens <- dens / max(dens, na.rm = TRUE) * width
      
      y <- 10^d$x
      
      tibble(
        hab_class = .y$hab_class,
        time      = .y$time,
        x         = c(x0 - dens, rev(x0 + dens)),
        y         = c(y, rev(y))
      )
    }) %>%
    bind_rows()
}

vio_poly <- build_weighted_violin(patch_tbl, width = 0.38, n = 700)

mu_plot <- mu_tbl_long %>%
  mutate(
    time = factor(time_raw, levels = c("before","after"), labels = unname(time_labs)),
    x    = unname(x_map[as.character(time)])
  ) %>%
  select(hab_class, time, x, mu_w_ha)

ann_mu <- mu_tbl %>%
  mutate(
    label = sprintf(
      "μ: %.1f → %.1f ha\n(%+.0f%%)",
      mu_w_ha_before, mu_w_ha_after, pct_mu_w
    )
  ) %>%
  transmute(hab_class, label) %>%
  distinct(hab_class, .keep_all = TRUE)

y_max <- vio_poly %>%
  group_by(hab_class) %>%
  summarise(y = max(y, na.rm = TRUE) * 1.35, .groups = "drop")

ann_mu <- left_join(ann_mu, y_max, by = "hab_class") %>%
  mutate(x = 1.5)

p_violin <- ggplot() +
  geom_polygon(
    data = vio_poly,
    aes(x = x, y = y, group = interaction(hab_class, time), fill = hab_class),
    color = "grey20",
    alpha = 0.55,
    linewidth = 0.4
  ) +
  facet_wrap(~ hab_class, nrow = 1) +
  scale_fill_manual(values = pal_thr, guide = "none") +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Pre-development", "Post-development"),
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  scale_y_log10(
    limits = c(0.01, 5000),
    breaks = c(0.01, 0.1, 1, 10, 50, 100, 500, 1000, 5000),
    labels = c("0.01","0.1","1","10","50","100","500","1000","5000"),
    expand = expansion(mult = c(0.03, 0.01))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Patch area (ha, log scale)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position    = "none",
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background   = element_rect(fill = "grey95"),
    strip.text         = element_text(face = "bold"),
    plot.margin        = margin(t = 12, r = 8, b = 8, l = 8)
  ) +
  geom_text(
    data = ann_mu,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.6,
    lineheight = 0.95
  ) +
  geom_point(
    data = mu_plot,
    aes(x = x, y = mu_w_ha),
    inherit.aes = FALSE,
    shape = 21,
    size = 3.2,
    stroke = 0.6,
    fill = "black",
    color = "black"
  ) +
  geom_text(
    data = mu_plot,
    aes(x = x + 0.04, y = mu_w_ha, label = "μ"),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold",
    color = "black"
  )

print(p_violin)

ggsave(
  filename = file.path(dir_fig, "patch_violin_area_weighted.png"),
  plot     = p_violin,
  width    = 13,
  height   = 4.8,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

# ---- 8) CONSTRAINT SPACE + WEIGHT-DEPENDENT SCENARIOS (WIND 150 m) ----------
# Goal: evaluate spatial alternatives to minimize habitat impacts by projecting
#       each pixel x -> z(x) = (habitat priority, wind, habitat aggregation),
#       then ranking pixels by a composite score and selecting a fixed area
#       (spatial budget).

# 8.1 Parameters ---------------------------------------------------------------

terra::terraOptions(memfrac = 0.7, todisk = FALSE, progress = 1)
gc()

# Habitat priority axis (sigmoid) — center at mitigation threshold
H_hp_center <- mitigation_threshold
H_hp_scale  <- 0.03

# Local aggregation window (female home-range proxy; 37x37 cells on 10-m grid)
agg_win <- 37

# Wind eligibility + continuous opportunity ramp (150 m layer)
W_min <- 8.0
W_eps <- 0.05

# Upper reference for wind opportunity: top 10% of wind speeds within AOI (>= W_min)
wind_top_q <- 0.90

# AOI layer name in the GDB
layer_aoi <- "Seigneurie_de_Beaupre"

# Iterative allocation parameters
tol_ha   <- 0.5
max_iter <- 200

# Weight scenarios
weight_scenarios <- tibble::tribble(
  ~scenario,       ~w_wind, ~w_hp, ~w_agg,
  "reference",      0.50,    0.35,  0.15,
  "development",    0.70,    0.20,  0.10,
  "conservation",   0.30,    0.50,  0.20
)

if (any(!is.finite(as.matrix(weight_scenarios[, c("w_wind", "w_hp", "w_agg")])))) {
  stop("Non-finite values found in weight scenarios.")
}
if (any(as.matrix(weight_scenarios[, c("w_wind", "w_hp", "w_agg")]) < 0)) {
  stop("Negative weights are not allowed.")
}

# 8.2 Inputs: habitat (pre-development) + wind (150 m) -------------------------

r_hab <- r_predev_full
names(r_hab) <- "habitat_predev_occprob"

wind_path_150m <- file.path(dir_data, "wind_speed_150m_SdB.tif")
if (!file.exists(wind_path_150m)) {
  stop("Missing wind raster: ", wind_path_150m)
}

r_wind <- terra::rast(wind_path_150m)
names(r_wind) <- "wind_150m"

aoi_sf <- sf::st_read(dsn = gdb_path, layer = layer_aoi, quiet = TRUE) |>
  sf::st_make_valid() |>
  sf::st_transform(terra::crs(r_hab))

aoi_v <- terra::vect(sf::st_union(aoi_sf))

r_hab  <- terra::mask(terra::crop(r_hab,  aoi_v), aoi_v)
r_wind <- terra::mask(terra::crop(r_wind, aoi_v), aoi_v)

if (!terra::compareGeom(r_hab, r_wind, stopOnError = FALSE)) {
  r_wind <- terra::resample(r_wind, r_hab, method = "near")
  names(r_wind) <- "wind_150m"
}

# 8.3 Utilities ----------------------------------------------------------------

clamp01 <- function(r) terra::clamp(r, 0, 1)

logistic_r <- function(r, t, s) {
  1 / (1 + exp(-(r - t) / s))
}

ramp01 <- function(r, a, b) {
  clamp01((r - a) / (b - a))
}

ramp01_floor <- function(r, a, b, eps = 0.05) {
  u <- clamp01((r - a) / (b - a))
  eps + (1 - eps) * u
}

hab_class_2thr <- function(r_hab, thr_mit, thr_con) {
  terra::ifel(
    is.na(r_hab), NA,
    terra::ifel(
      r_hab < thr_mit, 0L,
      terra::ifel(r_hab < thr_con, 1L, 2L)
    )
  )
}

ca_by_class_ha <- function(r_cls_int) {
  x <- landscapemetrics::lsm_c_ca(r_cls_int)
  x <- as.data.frame(x)
  x$class <- as.integer(x$class)
  x$value <- as.numeric(x$value)
  x <- x[x$class %in% c(1L, 2L), c("class", "value")]
  names(x) <- c("class", "ca_ha")
  x <- x[order(x$class), ]
  x
}

make_Adev_from_weights <- function(z_wind, z_hp, z_agg, w_wind, w_hp, w_agg) {
  
  S_raw_i <- w_wind * z_wind - w_hp * z_hp - w_agg * z_agg
  names(S_raw_i) <- "S_raw"
  
  S_min_i <- terra::global(S_raw_i, "min", na.rm = TRUE)[1, 1]
  S_max_i <- terra::global(S_raw_i, "max", na.rm = TRUE)[1, 1]
  
  if (!is.finite(S_min_i) || !is.finite(S_max_i) || (S_max_i - S_min_i) <= 0) {
    stop("Invalid score range during weight-dependent analysis.")
  }
  
  A_i <- (S_raw_i - S_min_i) / (S_max_i - S_min_i)
  A_i <- clamp01(A_i)
  names(A_i) <- "A_dev"
  
  list(
    S_raw = S_raw_i,
    A_dev = A_i,
    S_min = S_min_i,
    S_max = S_max_i
  )
}

select_top_ha <- function(A_dev, ha_target, eligible = NULL, seed = 1L) {
  if (!inherits(A_dev, "SpatRaster")) {
    stop("A_dev must be a SpatRaster.")
  }
  if (!is.finite(ha_target) || ha_target <= 0) {
    stop("ha_target must be a positive finite number.")
  }
  
  adm <- if (!is.null(eligible)) eligible else terra::ifel(is.na(A_dev), NA, 1)
  A_adm <- terra::mask(A_dev, adm)
  
  pix_ha <- prod(terra::res(A_adm)) / 10000
  n_req  <- ceiling(ha_target / pix_ha)
  
  df <- terra::as.data.frame(A_adm, cells = TRUE, na.rm = TRUE)
  if (!nrow(df)) {
    stop("No selectable cells are available in the admissible domain.")
  }
  names(df) <- c("cell", "A")
  
  n_req <- min(n_req, nrow(df))
  df <- df[order(df$A, decreasing = TRUE), , drop = FALSE]
  tau <- df$A[n_req]
  
  idx_strict <- which(df$A > tau)
  idx_ties   <- which(df$A == tau)
  n_need     <- n_req - length(idx_strict)
  
  if (n_need < 0) {
    stop("Internal error: negative number of tie cells requested.")
  }
  
  if (n_need > 0) {
    if (!is.null(seed)) set.seed(seed)
    idx_pick <- sample(idx_ties, size = n_need, replace = FALSE)
    sel_cells <- c(df$cell[idx_strict], df$cell[idx_pick])
  } else {
    sel_cells <- df$cell[idx_strict]
  }
  
  A_bin <- terra::ifel(is.na(A_adm), NA, NA)
  A_bin[sel_cells] <- 1L
  names(A_bin) <- "selected"
  
  list(
    tau    = tau,
    A_bin  = A_bin,
    pixels = length(sel_cells),
    ha     = length(sel_cells) * pix_ha,
    pix_ha = pix_ha,
    n_req  = n_req
  )
}

extract_loss_from_selection <- function(r_sel, r_hab, thr_mit, thr_con) {
  
  r_cls_before <- terra::ifel(
    is.na(r_hab), NA,
    terra::ifel(
      r_hab < thr_mit, 0L,
      terra::ifel(r_hab < thr_con, 1L, 2L)
    )
  )
  r_cls_before <- terra::as.int(r_cls_before)
  
  r_cls_after <- r_cls_before
  r_cls_after[!is.na(r_sel)] <- 0L
  r_cls_after <- terra::as.int(r_cls_after)
  
  ca_before <- ca_by_class_ha(r_cls_before)
  ca_after  <- ca_by_class_ha(r_cls_after)
  
  before_vec <- setNames(rep(0, 2), c("1", "2"))
  after_vec  <- setNames(rep(0, 2), c("1", "2"))
  
  before_vec[as.character(ca_before$class)] <- ca_before$ca_ha
  after_vec[as.character(ca_after$class)]   <- ca_after$ca_ha
  
  tibble::tibble(
    loss_mitigation_ha = before_vec["1"] - after_vec["1"],
    loss_protection_ha = before_vec["2"] - after_vec["2"]
  )
}

run_fixed_budget_selection <- function(A_dev, area_target_ha, eligible_wind,
                                       r_hab, thr_mit, thr_con, seed = 1L) {
  
  out_sel <- select_top_ha(
    A_dev     = A_dev,
    ha_target = area_target_ha,
    eligible  = eligible_wind,
    seed      = seed
  )
  
  r_sel <- terra::ifel(out_sel$A_bin == 1L, 1L, NA)
  names(r_sel) <- "selection"
  
  loss_tbl <- extract_loss_from_selection(
    r_sel   = r_sel,
    r_hab   = r_hab,
    thr_mit = thr_mit,
    thr_con = thr_con
  )
  
  list(
    tau              = out_sel$tau,
    selected_area_ha = out_sel$ha,
    selection_raster = r_sel,
    loss_summary     = loss_tbl
  )
}

run_iterative_loss_trajectory <- function(A_dev, r_hab, eligible_wind,
                                          thr_mit, thr_con,
                                          block_ha, tol_ha = 0.5,
                                          max_iter = 200,
                                          save_iters = NULL) {
  
  if (!terra::compareGeom(A_dev, r_hab, stopOnError = FALSE)) {
    stop("Geometry mismatch between A_dev and r_hab.")
  }
  
  pix_ha <- prod(terra::res(A_dev)) / 10000
  
  adm_base <- terra::ifel(is.na(A_dev), NA, 1L)
  adm_base <- terra::mask(adm_base, eligible_wind)
  adm_base <- terra::as.int(adm_base)
  names(adm_base) <- "adm_base"
  
  r_cls0 <- hab_class_2thr(r_hab, thr_mit, thr_con)
  r_cls0 <- terra::mask(r_cls0, adm_base)
  r_cls0 <- terra::as.int(r_cls0)
  names(r_cls0) <- "hab_class"
  
  before_ca <- ca_by_class_ha(r_cls0)
  before_vec <- setNames(rep(0, 2), c("1", "2"))
  before_vec[as.character(before_ca$class)] <- before_ca$ca_ha
  target_loss <- before_vec
  
  sel_cum <- terra::ifel(is.na(adm_base), NA, 0L)
  sel_cum <- terra::as.int(sel_cum)
  names(sel_cum) <- "sel_cum"
  
  r_imp <- r_cls0
  names(r_imp) <- "hab_class_imp"
  
  log_tbl <- list()
  sel_snapshots <- list()
  
  for (iter in seq_len(max_iter)) {
    
    eligible_i <- terra::ifel(
      is.na(adm_base), NA,
      terra::ifel(sel_cum == 1L, NA, 1L)
    )
    eligible_i <- terra::as.int(eligible_i)
    names(eligible_i) <- "eligible_i"
    
    n_avail <- terra::global(!is.na(eligible_i), "sum", na.rm = TRUE)[1, 1]
    if (!is.finite(n_avail) || n_avail == 0) break
    
    out_i <- select_top_ha(
      A_dev     = A_dev,
      ha_target = block_ha,
      eligible  = eligible_i,
      seed      = 1L + iter
    )
    
    sel_i <- terra::ifel(out_i$A_bin == 1L, 1L, NA)
    sel_i <- terra::as.int(sel_i)
    names(sel_i) <- "sel_i"
    
    sel_cum <- terra::ifel(
      is.na(sel_cum), NA,
      terra::ifel(sel_i == 1L, 1L, sel_cum)
    )
    sel_cum <- terra::as.int(sel_cum)
    names(sel_cum) <- "sel_cum"
    
    r_imp[!is.na(sel_i)] <- 0L
    r_imp <- terra::as.int(r_imp)
    names(r_imp) <- "hab_class_imp"
    
    after_ca <- ca_by_class_ha(r_imp)
    after_vec <- setNames(rep(0, 2), c("1", "2"))
    after_vec[as.character(after_ca$class)] <- after_ca$ca_ha
    
    loss_vec <- before_vec - after_vec
    ok <- (loss_vec >= (target_loss - tol_ha))
    
    n_sel  <- terra::global(sel_cum == 1L, "sum", na.rm = TRUE)[1, 1]
    ha_cum <- n_sel * pix_ha
    
    log_tbl[[iter]] <- tibble::tibble(
      iteration          = iter,
      tau                = out_i$tau,
      ha_block           = out_i$ha,
      ha_cum             = ha_cum,
      loss_mitigation_ha = loss_vec["1"],
      loss_protection_ha = loss_vec["2"],
      hit_mitigation     = isTRUE(ok["1"]),
      hit_protection     = isTRUE(ok["2"])
    )
    
    if (!is.null(save_iters) && iter %in% save_iters) {
      sel_snapshots[[paste0("iter_", iter)]] <- sel_cum
    }
    
    stop_now <- all(ok) ||
      (is.finite(n_avail) && n_avail * pix_ha < (block_ha * 0.5))
    
    rm(eligible_i, sel_i, out_i, after_ca, after_vec, loss_vec, ok, n_sel, ha_cum)
    gc()
    
    if (stop_now) break
  }
  
  log_df <- dplyr::bind_rows(log_tbl)
  
  list(
    adm_base      = adm_base,
    r_cls0        = r_cls0,
    before_vec    = before_vec,
    log_df        = log_df,
    loss_table    = log_df %>%
      dplyr::transmute(
        iteration,
        ha_cum = round(ha_cum, 0),
        loss_mitigation_ha = round(loss_mitigation_ha, 0),
        loss_protection_ha = round(loss_protection_ha, 0)
      ),
    sel_cum       = sel_cum,
    sel_snapshots = sel_snapshots
  )
}

run_full_weight_scenario <- function(scenario_name, w_wind, w_hp, w_agg,
                                     z2_wind, z1_hp, z3_agg,
                                     eligible_wind, r_hab,
                                     mitigation_threshold, protection_threshold,
                                     area_target_ha, block_ha,
                                     tol_ha = 0.5, max_iter = 200,
                                     save_iters = c(8L, 15L),
                                     seed_fixed = 1L) {
  
  score_obj <- make_Adev_from_weights(
    z_wind = z2_wind,
    z_hp   = z1_hp,
    z_agg  = z3_agg,
    w_wind = w_wind,
    w_hp   = w_hp,
    w_agg  = w_agg
  )
  
  fixed_obj <- run_fixed_budget_selection(
    A_dev          = score_obj$A_dev,
    area_target_ha = area_target_ha,
    eligible_wind  = eligible_wind,
    r_hab          = r_hab,
    thr_mit        = mitigation_threshold,
    thr_con        = protection_threshold,
    seed           = seed_fixed
  )
  
  iter_obj <- run_iterative_loss_trajectory(
    A_dev         = score_obj$A_dev,
    r_hab         = r_hab,
    eligible_wind = eligible_wind,
    thr_mit       = mitigation_threshold,
    thr_con       = protection_threshold,
    block_ha      = block_ha,
    tol_ha        = tol_ha,
    max_iter      = max_iter,
    save_iters    = save_iters
  )
  
  list(
    scenario = scenario_name,
    weights = tibble::tibble(
      scenario = scenario_name,
      w_wind = w_wind,
      w_hp   = w_hp,
      w_agg  = w_agg
    ),
    score = score_obj,
    fixed = fixed_obj,
    iterative = iter_obj
  )
}

selection_overlap_pct <- function(r_sel_a, r_sel_b) {
  a <- terra::values(r_sel_a, mat = FALSE)
  b <- terra::values(r_sel_b, mat = FALSE)
  
  a <- !is.na(a) & a == 1
  b <- !is.na(b) & b == 1
  
  n_union <- sum(a | b, na.rm = TRUE)
  if (n_union == 0) return(NA_real_)
  
  100 * sum(a & b, na.rm = TRUE) / n_union
}

topk_overlap_pct <- function(A_ref, A_cmp, probs = c(0.05, 0.10, 0.20)) {
  
  v_ref <- terra::values(A_ref, mat = FALSE)
  v_cmp <- terra::values(A_cmp, mat = FALSE)
  
  ok <- is.finite(v_ref) & !is.na(v_ref) & is.finite(v_cmp) & !is.na(v_cmp)
  v_ref <- v_ref[ok]
  v_cmp <- v_cmp[ok]
  
  n <- length(v_ref)
  if (n == 0) {
    return(tibble::tibble(
      top_fraction = probs,
      overlap_pct  = NA_real_
    ))
  }
  
  out <- lapply(probs, function(p) {
    k <- max(1, ceiling(n * p))
    idx_ref <- order(v_ref, decreasing = TRUE)[seq_len(k)]
    idx_cmp <- order(v_cmp, decreasing = TRUE)[seq_len(k)]
    
    tibble::tibble(
      top_fraction = p,
      overlap_pct  = 100 * length(intersect(idx_ref, idx_cmp)) / k
    )
  })
  
  dplyr::bind_rows(out)
}

make_forest_plot <- function(df_plot, plot_title = NULL, reordered = TRUE,
                             x_lim_common, x_breaks_common,
                             force_lim_common, force_breaks_common,
                             pal_class) {
  
  p <- ggplot(
    df_plot,
    aes(x = delta_pct, y = metric_lab, color = class)
  ) +
    annotate(
      "rect",
      xmin = -5, xmax = 5,
      ymin = -Inf, ymax = Inf,
      alpha = 0.30,
      fill  = "grey60"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
    geom_point(aes(alpha = out_5, size = force)) +
    scale_alpha_manual(
      values = c(`FALSE` = 0.65, `TRUE` = 1),
      guide = "none"
    ) +
    scale_size_continuous(
      limits = force_lim_common,
      breaks = force_breaks_common,
      range  = c(2.5, 10),
      name   = "Change magnitude (%)",
      guide  = guide_legend(
        override.aes = list(color = "grey30", alpha = 1)
      )
    ) +
    facet_wrap(~ category, scales = "free_y", ncol = 1) +
    scale_color_manual(
      values = pal_class,
      name   = "Habitat class",
      drop   = FALSE,
      guide  = guide_legend(
        override.aes = list(size = 6, alpha = 1)
      )
    ) +
    scale_x_continuous(
      limits = x_lim_common,
      breaks = x_breaks_common
    ) +
    labs(
      x = "Δ (%)",
      y = NULL,
      title = plot_title
    ) +
    theme_bw(base_size = 20) +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 20, face = "bold"),
      legend.text  = element_text(size = 18),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 22, face = "bold"),
      axis.text.x  = element_text(size = 18),
      axis.text.y  = element_text(size = 18),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4)
    )
  
  if (reordered) {
    p <- p + tidytext::scale_y_reordered()
  } else {
    p <- p + scale_y_discrete()
  }
  
  return(p)
}

# 8.4 Exclude non-exploitable wind ---------------------------------------------
# Keep full habitat information over the AOI for habitat-priority calculations.
# Build the wind-eligible domain separately, and apply it only after all
# habitat-priority surfaces have been computed.

eligible_wind <- terra::ifel(
  is.na(r_wind),
  NA,
  terra::ifel(r_wind >= W_min, 1, NA)
)
names(eligible_wind) <- "eligible_wind"

# Preserve full-AOI versions for habitat-priority calculations
r_hab_full  <- r_hab
r_wind_full <- r_wind

# Wind-masked version used only where the admissible development domain is required
r_wind_elig <- terra::mask(r_wind_full, eligible_wind)
names(r_wind_elig) <- "wind_150m_eligible"

# 8.5 z(x) = (habitat priority, wind, habitat aggregation) ---------------------

# Habitat priority from full habitat surface
z1_hp <- clamp01(logistic_r(r_hab_full, t = H_hp_center, s = H_hp_scale))
names(z1_hp) <- "z1_hp"

# Wind opportunity from eligible wind values only
wind_vals <- terra::values(r_wind_elig, mat = FALSE, na.rm = TRUE)
if (!length(wind_vals)) {
  stop("No valid wind values found after masking to the AOI and eligibility domain.")
}

W_best <- as.numeric(stats::quantile(wind_vals, probs = wind_top_q, na.rm = TRUE))
W_best <- max(W_best, W_min + 1e-6)

z2_wind <- terra::ifel(
  is.na(r_wind_full), NA,
  terra::ifel(
    r_wind_full < W_min, 0,
    terra::ifel(
      r_wind_full <= W_best,
      ramp01_floor(r_wind_full, a = W_min, b = W_best, eps = W_eps),
      1
    )
  )
)
names(z2_wind) <- "z2_wind"

# Local aggregation must be computed on the full habitat landscape in the AOI,
# not on the wind-masked domain
k <- matrix(1, nrow = agg_win, ncol = agg_win)

hab_ge_mit <- terra::ifel(
  is.na(r_hab_full), NA,
  terra::ifel(r_hab_full >= mitigation_threshold, 1, 0)
)
names(hab_ge_mit) <- "hab_ge_mitigation"

z3_agg <- terra::focal(
  hab_ge_mit,
  w      = k,
  fun    = mean,
  na.rm  = TRUE,
  expand = TRUE
)
z3_agg <- clamp01(z3_agg)
names(z3_agg) <- "z3_agg"

# Final admissible domain for the constraint-space analysis
z1_hp   <- terra::mask(z1_hp,   eligible_wind)
z2_wind <- terra::mask(z2_wind, eligible_wind)
z3_agg  <- terra::mask(z3_agg,  eligible_wind)

# Keep the analysis rasters restricted to the same admissible domain afterward
r_hab  <- terra::mask(r_hab_full,  eligible_wind)
r_wind <- terra::mask(r_wind_full, eligible_wind)

gc()

# 8.6 Spatial budget -----------------------------------------------------------

if (!exists("footprint_all_sf") && !exists("footprint_all_v")) {
  stop("A total wind-footprint object is required: footprint_all_sf or footprint_all_v.")
}

if (exists("footprint_all_sf")) {
  fp_sf <- footprint_all_sf
  if (is.na(sf::st_crs(fp_sf))) stop("footprint_all_sf has missing CRS.")
  fp_sf <- sf::st_transform(fp_sf, sf::st_crs(terra::crs(z1_hp)))
  fp_area_m2 <- as.numeric(sf::st_area(fp_sf))
  fp_area_ha <- sum(fp_area_m2, na.rm = TRUE) / 10000
} else {
  fp_v <- footprint_all_v
  fp_area_ha <- sum(terra::expanse(fp_v, unit = "ha", transform = FALSE), na.rm = TRUE)
}

if (!is.finite(fp_area_ha) || fp_area_ha <= 0) {
  stop("Invalid total wind-footprint area.")
}

area_target_ha <- 2 * fp_area_ha
block_ha <- fp_area_ha

cat("\nSpatial budget definition\n")
cat("Total wind footprint area (ha):", round(fp_area_ha, 2), "\n")
cat("Budget (2× footprint) (ha)     :", round(area_target_ha, 2), "\n")

# 8.7 Run all weight scenarios -------------------------------------------------

baseline_name <- "reference"
scenario_results <- list()

for (i in seq_len(nrow(weight_scenarios))) {
  
  sc <- weight_scenarios[i, ]
  
  cat("\nRunning scenario:", sc$scenario, "\n")
  
  keep_detail <- identical(sc$scenario, baseline_name)
  
  scenario_results[[sc$scenario]] <- run_full_weight_scenario(
    scenario_name = sc$scenario,
    w_wind = sc$w_wind,
    w_hp   = sc$w_hp,
    w_agg  = sc$w_agg,
    z2_wind = z2_wind,
    z1_hp   = z1_hp,
    z3_agg  = z3_agg,
    eligible_wind = eligible_wind,
    r_hab = r_hab,
    mitigation_threshold = mitigation_threshold,
    protection_threshold = protection_threshold,
    area_target_ha = area_target_ha,
    block_ha = block_ha,
    tol_ha = tol_ha,
    max_iter = max_iter,
    save_iters = if (keep_detail) c(8L, 15L) else NULL,
    seed_fixed = 1L
  )
  
  write_csv_bom(
    scenario_results[[sc$scenario]]$iterative$loss_table,
    file.path(dir_tab, paste0("loss_table_", sc$scenario, ".csv"))
  )
  
  if (!keep_detail) {
    scenario_results[[sc$scenario]]$score$S_raw <- NULL
    
    scenario_results[[sc$scenario]]$iterative$adm_base      <- NULL
    scenario_results[[sc$scenario]]$iterative$r_cls0        <- NULL
    scenario_results[[sc$scenario]]$iterative$before_vec    <- NULL
    scenario_results[[sc$scenario]]$iterative$log_df        <- NULL
    scenario_results[[sc$scenario]]$iterative$sel_cum       <- NULL
    scenario_results[[sc$scenario]]$iterative$sel_snapshots <- NULL
  }
  
  gc()
}

# 8.8 Selected scenario objects + scenario exports -----------------------------

baseline_name <- "reference"

if (!(baseline_name %in% names(scenario_results))) {
  stop("Reference scenario not found in scenario_results.")
}

ref <- scenario_results[[baseline_name]]

A_dev <- ref$score$A_dev
S_raw <- ref$score$S_raw
S_min <- ref$score$S_min
S_max <- ref$score$S_max

tau        <- ref$fixed$tau
A_bin      <- ref$fixed$selection_raster
sel_out_ha <- ref$fixed$selected_area_ha

loss_table    <- ref$iterative$loss_table
sel_snapshots <- ref$iterative$sel_snapshots
adm_base      <- ref$iterative$adm_base

# Export scenario-specific normalized A_dev rasters
for (sc_name in names(scenario_results)) {
  
  A_dev_i <- scenario_results[[sc_name]]$score$A_dev
  
  terra::writeRaster(
    A_dev_i,
    file.path(dir_ras, paste0("A_dev_score_", sc_name, "_wind150m.tif")),
    overwrite = TRUE
  )
}

# Build common-scale comparable rasters across all scenarios
S_raw_ref <- scenario_results[["reference"]]$score$S_raw

S_raw_dev <- with(
  scenario_results[["development"]]$weights,
  w_wind * z2_wind - w_hp * z1_hp - w_agg * z3_agg
)
names(S_raw_dev) <- "S_raw_development"

S_raw_con <- with(
  scenario_results[["conservation"]]$weights,
  w_wind * z2_wind - w_hp * z1_hp - w_agg * z3_agg
)
names(S_raw_con) <- "S_raw_conservation"

S_all_min <- min(
  terra::global(S_raw_ref, "min", na.rm = TRUE)[1, 1],
  terra::global(S_raw_dev, "min", na.rm = TRUE)[1, 1],
  terra::global(S_raw_con, "min", na.rm = TRUE)[1, 1],
  na.rm = TRUE
)

S_all_max <- max(
  terra::global(S_raw_ref, "max", na.rm = TRUE)[1, 1],
  terra::global(S_raw_dev, "max", na.rm = TRUE)[1, 1],
  terra::global(S_raw_con, "max", na.rm = TRUE)[1, 1],
  na.rm = TRUE
)

if (!is.finite(S_all_min) || !is.finite(S_all_max) || (S_all_max - S_all_min) <= 0) {
  stop("Invalid common score range across scenarios.")
}

A_ref_common <- clamp01((S_raw_ref - S_all_min) / (S_all_max - S_all_min))
A_dev_common <- clamp01((S_raw_dev - S_all_min) / (S_all_max - S_all_min))
A_con_common <- clamp01((S_raw_con - S_all_min) / (S_all_max - S_all_min))

names(A_ref_common) <- "A_dev_common"
names(A_dev_common) <- "A_dev_common"
names(A_con_common) <- "A_dev_common"

terra::writeRaster(
  A_ref_common,
  file.path(dir_ras, "A_dev_common_reference_wind150m.tif"),
  overwrite = TRUE
)

terra::writeRaster(
  A_dev_common,
  file.path(dir_ras, "A_dev_common_development_wind150m.tif"),
  overwrite = TRUE
)

terra::writeRaster(
  A_con_common,
  file.path(dir_ras, "A_dev_common_conservation_wind150m.tif"),
  overwrite = TRUE
)

terra::writeRaster(
  A_bin,
  file.path(dir_ras, sprintf("selection_top_%dha_wind150m.tif", round(area_target_ha))),
  overwrite = TRUE
)

# 8.9 Selected scenario controls + histogram ---------------------------------

vals_A <- terra::values(A_dev, mat = FALSE, na.rm = TRUE)
df_hist <- data.frame(A_dev = as.numeric(vals_A))

p_hist <- ggplot(df_hist, aes(x = A_dev)) +
  geom_histogram(bins = 60, colour = "grey20", fill = "grey70") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Distribution of the constraint-space score",
    x = "A_dev (0–1)",
    y = "Pixel count"
  ) +
  theme_bw()

print(p_hist)

cat("\n==================== CONSTRAINT SPACE CONTROLS ====================\n")
cat("Reference scenario        :", baseline_name, "\n")
cat("Wind layer (150 m)        :", basename(wind_path_150m), "\n")
cat("Eligibility threshold     :", W_min, "m/s\n")
cat("Wind reference (W_best)   :", round(W_best, 3), "m/s (q=", wind_top_q, ")\n", sep = "")
cat("S_min / S_max             :", S_min, "/", S_max, "\n")
cat("A_dev min/max             :",
    terra::global(A_dev, "min", na.rm = TRUE)[1, 1], "/",
    terra::global(A_dev, "max", na.rm = TRUE)[1, 1], "\n")
cat("Eligible pixels           :",
    terra::global(!is.na(eligible_wind), "sum", na.rm = TRUE)[1, 1], "/",
    terra::ncell(eligible_wind), "\n")
cat("Selected area (ha)        :", round(sel_out_ha, 2), "\n")
cat("Tau (A_dev)               :", round(tau, 6), "\n")
cat("==============================================================\n")

# 8.10 Selected scenario constraint-space visualization ------------------------

set.seed(1)

n_samp <- 500000
n_plot <- 80000

samp <- terra::spatSample(
  A_dev,
  size      = n_samp,
  method    = "random",
  cells     = TRUE,
  na.rm     = TRUE,
  as.points = FALSE
)

idx <- samp[[1]]
if (!length(idx)) {
  stop("No eligible cells were sampled from A_dev.")
}

# Build a raster stack containing the four variables used in the constraint space
xyz_stack <- c(A_dev, z1_hp, z2_wind, z3_agg)
names(xyz_stack) <- c("A", "habitat_priority", "wind", "habitat_aggregation")

# Get cell-center coordinates for the sampled cells
xy <- terra::xyFromCell(A_dev, idx)

# Extract values by indexing the full value matrix with sampled cell ids
vals_mat <- terra::values(xyz_stack, mat = TRUE)
vals <- vals_mat[idx, , drop = FALSE]

df3 <- data.frame(
  x                   = xy[, 1],
  y                   = xy[, 2],
  A                   = vals[, "A"],
  habitat_priority    = vals[, "habitat_priority"],
  wind                = vals[, "wind"],
  habitat_aggregation = vals[, "habitat_aggregation"]
)

# Keep only fully valid rows
df3 <- df3[
  is.finite(df3$x) &
    is.finite(df3$y) &
    is.finite(df3$A) &
    is.finite(df3$habitat_priority) &
    is.finite(df3$wind) &
    is.finite(df3$habitat_aggregation),
]

# Clamp all constraint-space values to [0, 1]
df3$A                   <- pmin(pmax(df3$A,                   0), 1)
df3$habitat_priority    <- pmin(pmax(df3$habitat_priority,    0), 1)
df3$wind                <- pmin(pmax(df3$wind,                0), 1)
df3$habitat_aggregation <- pmin(pmax(df3$habitat_aggregation, 0), 1)

if (!nrow(df3)) {
  stop("No valid sampled points are available for constraint-space plotting.")
}

# Subsample points used for plotting
n_take <- min(n_plot, nrow(df3))
set.seed(1)
df3_plot <- df3[sample.int(nrow(df3), n_take), ]

# 3D constraint space
p3d <- plotly::plot_ly(
  data = df3_plot,
  x = ~habitat_priority,
  y = ~wind,
  z = ~habitat_aggregation,
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = 2,
    opacity = 0.45,
    color = ~A,
    colorscale = "Viridis",
    cmin = 0,
    cmax = 1,
    showscale = TRUE,
    colorbar = list(title = "A_dev")
  )
) |>
  plotly::layout(
    scene = list(
      xaxis = list(title = "Habitat priority", range = c(0, 1)),
      yaxis = list(title = "Wind opportunity", range = c(0, 1)),
      zaxis = list(title = "Habitat aggregation", range = c(0, 1))
    )
  )

p3d

htmlwidgets::saveWidget(
  p3d,
  file.path(dir_int, "constraint_space_cube_wind150m.html"),
  selfcontained = TRUE
)

# 2D spatial map of the same sampled pixels

seigneurie_sf <- aoi_sf |>
  sf::st_make_valid() |>
  sf::st_union() |>
  sf::st_as_sf()

seigneurie_outline <- sf::st_boundary(seigneurie_sf)

p_map_sampled <- ggplot() +
  geom_point(
    data = df3_plot,
    aes(x = x, y = y, color = A),
    size = 0.20,
    alpha = 0.80
  ) +
  geom_sf(
    data = seigneurie_outline,
    inherit.aes = FALSE,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  coord_sf(expand = FALSE) +
  scale_color_viridis_c(
    option = "D",
    limits = c(0, 1),
    name = "A_dev"
  ) +
  labs(
    title = "Sampled pixels used in the constraint space"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

print(p_map_sampled)

ggsave(
  filename = file.path(dir_fig, "constraint_space_sampled_pixels_map_wind150m.png"),
  plot     = p_map_sampled,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 300,
  bg       = "white",
  limitsize = FALSE
)

# 2D projections of the constraint space (three orthogonal slices)

p_sol_hp_wind <- ggplot(
  df3_plot,
  aes(x = habitat_priority, y = wind, color = A)
) +
  geom_point(size = 0.35, alpha = 0.65) +
  scale_color_viridis_c(
    option = "D",
    limits = c(0, 1),
    name   = "A_dev"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Habitat priority × wind",
    x = "Habitat priority",
    y = "Wind opportunity"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

p_sol_hp_agg <- ggplot(
  df3_plot,
  aes(x = habitat_priority, y = habitat_aggregation, color = A)
) +
  geom_point(size = 0.35, alpha = 0.65) +
  scale_color_viridis_c(
    option = "D",
    limits = c(0, 1),
    name   = "A_dev"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Habitat priority × habitat aggregation",
    x = "Habitat priority",
    y = "Habitat aggregation"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

p_sol_wind_agg <- ggplot(
  df3_plot,
  aes(x = wind, y = habitat_aggregation, color = A)
) +
  geom_point(size = 0.35, alpha = 0.65) +
  scale_color_viridis_c(
    option = "D",
    limits = c(0, 1),
    name   = "A_dev"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Wind × habitat aggregation",
    x = "Wind opportunity",
    y = "Habitat aggregation"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

p_sol2d <- (p_sol_hp_wind | p_sol_hp_agg) / p_sol_wind_agg
p_sol2d <- p_sol2d + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(p_sol2d)

ggsave(
  filename = file.path(dir_fig, "constraint_space_projections_wind150m.png"),
  plot     = p_sol2d,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

# 8.11 Selected scenario selection polygon and maps ----------------------------

# target = 2 × wind footprint

sel_poly_v <- terra::as.polygons(A_bin, dissolve = TRUE, values = TRUE, na.rm = TRUE)
nm <- names(sel_poly_v)[1]
sel_poly_v <- sel_poly_v[sel_poly_v[[nm]] == 1, ]

sel_poly_sf <- sf::st_as_sf(sel_poly_v) |>
  sf::st_make_valid()

gpkg_path <- file.path(
  dir_vec,
  sprintf("selection_top_%dha_wind150m.gpkg", round(area_target_ha))
)

sel_layer <- sprintf("top_%dha_wind150m", round(area_target_ha))

sf::st_write(sel_poly_sf, gpkg_path, layer = sel_layer, delete_layer = TRUE, quiet = TRUE)

rl_Adev <- raster::raster(A_dev)
names(rl_Adev) <- "A_dev"

rl_sel <- raster::raster(A_bin)
rl_sel <- raster::ratify(rl_sel)

rat <- raster::levels(rl_sel)[[1]]
rat <- rat[rat[, 1] == 1, , drop = FALSE]
rat$label <- sprintf("Top %d ha (A_dev ≥ %.3f)", round(area_target_ha), tau)

levels(rl_sel) <- rat
names(rl_sel) <- "Selection"

mapviewOptions(mapview.maxpixels = max(raster::ncell(rl_Adev), raster::ncell(rl_sel)))

mv_overlay_Adev <- mapview(
  rl_Adev,
  layer.name = "A_dev",
  legend     = TRUE,
  maxpixels  = 5e5,
  na.color   = "transparent"
) +
  mapview(
    rl_sel,
    layer.name    = "Selection",
    col.regions   = "#1B9E77",
    alpha.regions = 0.6,
    legend        = TRUE,
    maxpixels     = 5e5,
    na.color      = "transparent"
  )

mv_overlay_Adev

hab_cls <- terra::ifel(
  is.na(r_hab), NA,
  terra::ifel(
    r_hab < mitigation_threshold, 1L,
    terra::ifel(r_hab < protection_threshold, 2L, 3L)
  )
)
names(hab_cls) <- "Habitat_class"

rl_hab_cls <- raster::raster(hab_cls)
rl_hab_cls <- raster::ratify(rl_hab_cls)

rat_hab <- data.frame(
  ID = c(1, 2, 3),
  label = c("< mitigation", "mitigation", "protection")
)
levels(rl_hab_cls) <- rat_hab
names(rl_hab_cls) <- "Habitat (threshold classes)"

cols_hab <- c(
  "lightgrey",
  "#e8ac31",
  "#a86466"
)

mv_overlay_hab_cls <- mapview(
  rl_hab_cls,
  layer.name  = "Habitat threshold",
  col.regions = cols_hab,
  legend      = TRUE,
  maxpixels   = 5e5,
  na.color    = "transparent"
) +
  mapview(
    rl_sel,
    layer.name    = "Selection",
    col.regions   = "#1B9E77",
    alpha.regions = 0.6,
    legend        = TRUE,
    maxpixels     = 5e5,
    na.color      = "transparent"
  )

mv_overlay_hab_cls

# 8.12 Constraint-space distribution of the selected area ----------------------
# 2 × wind footprint

set.seed(1)
n_take <- min(n_plot, nrow(df3))

df3p <- df3[sample.int(nrow(df3), n_take), , drop = FALSE]
df3p$sel <- df3p$A >= tau

col_scale <- scale_color_viridis_c(
  option = "D",
  limits = c(0, 1),
  name   = "A_dev"
)

base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

p_hp_wind_sel <- ggplot(df3p, aes(habitat_priority, wind)) +
  geom_point(color = "grey70", alpha = 0.35, size = 0.6) +
  geom_point(data = df3p[df3p$sel, ], aes(color = A), alpha = 0.65, size = 0.7) +
  col_scale + base_theme +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Habitat priority", y = "Wind opportunity")

p_hp_agg_sel <- ggplot(df3p, aes(habitat_priority, habitat_aggregation)) +
  geom_point(color = "grey70", alpha = 0.35, size = 0.6) +
  geom_point(data = df3p[df3p$sel, ], aes(color = A), alpha = 0.65, size = 0.7) +
  col_scale + base_theme +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Habitat priority", y = "Habitat aggregation")

p_wind_agg_sel <- ggplot(df3p, aes(wind, habitat_aggregation)) +
  geom_point(color = "grey70", alpha = 0.35, size = 0.6) +
  geom_point(data = df3p[df3p$sel, ], aes(color = A), alpha = 0.65, size = 0.7) +
  col_scale + base_theme +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Wind opportunity", y = "Habitat aggregation")

p_sol_sel <- (p_hp_wind_sel | p_hp_agg_sel) / p_wind_agg_sel +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

p_sol_sel <- p_sol_sel +
  plot_annotation(
    title = sprintf("Constraint space: A_dev \u2265 tau (tau = %.3f)", tau)
  )

ggsave(
  filename = file.path(dir_fig, "constraint_space_selection_tau_wind150m.png"),
  plot     = p_sol_sel,
  width    = 11,
  height   = 8,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

# 8.13 Selected scenario iterative outputs ------------------------------------

sel_final <- terra::ifel(ref$iterative$sel_cum == 1L, 1L, NA)
names(sel_final) <- "sel_final_iterative_blocks"

sel_poly_v_iter <- terra::as.polygons(sel_final, dissolve = TRUE, values = TRUE, na.rm = TRUE)
nm_iter <- names(sel_poly_v_iter)[1]
sel_poly_v_iter <- sel_poly_v_iter[sel_poly_v_iter[[nm_iter]] == 1L, ]

sel_poly_sf_iter <- sf::st_as_sf(sel_poly_v_iter) |>
  sf::st_make_valid()

gpkg_path_iter <- file.path(dir_vec, "selection_cum_iter_wind150m.gpkg")

sf::st_write(
  sel_poly_sf_iter,
  gpkg_path_iter,
  layer = "selection_cum_iter_wind150m",
  delete_layer = TRUE,
  quiet = TRUE
)

write_csv_bom(
  ref$iterative$log_df,
  file.path(dir_tab, "iteration_log_blocks_fp_area_ha.csv")
)

print(ref$iterative$log_df, n = Inf)
print(loss_table, n = Inf)

write_csv_bom(
  loss_table,
  file.path(dir_tab, "cumulative_losses_by_class_ha.csv")
)

iter_targets <- c(8L, 15L)

mapviewOptions(mapview.maxpixels = 5e5)

if ("iter_8" %in% names(sel_snapshots)) {
  r_sel_8 <- terra::ifel(sel_snapshots[["iter_8"]] == 1L, 1L, NA)
  names(r_sel_8) <- "Selection"
  
  rl_sel_8 <- raster::raster(r_sel_8)
  rl_sel_8 <- raster::ratify(rl_sel_8)
  
  rat_sel_8 <- data.frame(
    ID = 1,
    label = "Cumulative selection to iteration 8"
  )
  levels(rl_sel_8) <- rat_sel_8
  names(rl_sel_8) <- "Selection"
  
  mv_overlay_iter8 <- mapview(
    rl_hab_cls,
    layer.name  = "Habitat threshold",
    col.regions = cols_hab,
    legend      = TRUE,
    maxpixels   = 5e5,
    na.color    = "transparent"
  ) +
    mapview(
      rl_sel_8,
      layer.name    = "Selection",
      col.regions   = "#1B9E77",
      alpha.regions = 0.6,
      legend        = TRUE,
      maxpixels     = 5e5,
      na.color      = "transparent"
    )
  
  mv_overlay_iter8
}

if ("iter_15" %in% names(sel_snapshots)) {
  seigneurie_sf <- sf::st_read(
    gdb_path,
    layer = "Seigneurie_de_Beaupre",
    quiet = TRUE
  ) |>
    sf::st_make_valid() |>
    sf::st_transform(terra::crs(hab_cls))
  
  seigneurie_sf <- sf::st_union(seigneurie_sf) |>
    sf::st_as_sf()
  
  seigneurie_outline <- sf::st_boundary(seigneurie_sf)
  
  r_sel_15 <- terra::ifel(sel_snapshots[["iter_15"]] == 1L, 1L, NA)
  names(r_sel_15) <- "Selection"
  
  rl_sel_15 <- raster::raster(r_sel_15)
  rl_sel_15 <- raster::ratify(rl_sel_15)
  
  rat_sel_15 <- data.frame(
    ID = 1,
    label = "Cumulative selection to iteration 15"
  )
  levels(rl_sel_15) <- rat_sel_15
  names(rl_sel_15) <- "Selection"
  
  mv_overlay_iter15 <- mapview(
    rl_hab_cls,
    layer.name  = "Habitat threshold",
    col.regions = cols_hab,
    legend      = TRUE,
    maxpixels   = 5e5,
    na.color    = "transparent"
  ) +
    mapview(
      rl_sel_15,
      layer.name    = "Selection",
      col.regions   = "#1B9E77",
      alpha.regions = 0.6,
      legend        = TRUE,
      maxpixels     = 5e5,
      na.color      = "transparent"
    ) +
    mapview(
      seigneurie_outline,
      layer.name = "Seigneurie boundary",
      color      = "black",
      lwd        = 2,
      alpha      = 1,
      legend     = FALSE
    )
  
  mv_overlay_iter15
  
  sel_poly_v_15 <- terra::as.polygons(r_sel_15, dissolve = TRUE, na.rm = TRUE)
  sel_poly_sf_15 <- sf::st_as_sf(sel_poly_v_15) |>
    sf::st_make_valid() |>
    sf::st_transform(sf::st_crs(seigneurie_sf))
  
  seigneurie_remaining <- sf::st_difference(seigneurie_sf, sel_poly_sf_15)
  area_remaining_ha <- units::set_units(sum(sf::st_area(seigneurie_remaining)), ha)
  area_remaining_ha
  
  ha_cum_iter15 <- loss_table %>%
    dplyr::filter(iteration == 15) %>%
    dplyr::pull(ha_cum)
  
  if (length(ha_cum_iter15) != 1) {
    stop("Could not extract a unique cumulative area for iteration 15.")
  }
  
  prop_sel <- ha_cum_iter15 / 159221 * 100
  round(prop_sel, 1)
}

# 8.14 Weight sensitivity summary ----------------------------------------------

scenario_summary_tbl <- dplyr::bind_rows(
  lapply(names(scenario_results), function(sc) {
    res <- scenario_results[[sc]]
    tibble::tibble(
      scenario           = sc,
      w_wind             = res$weights$w_wind,
      w_hp               = res$weights$w_hp,
      w_agg              = res$weights$w_agg,
      score_min          = res$score$S_min,
      score_max          = res$score$S_max,
      selected_area_ha   = res$fixed$selected_area_ha,
      tau                = res$fixed$tau,
      loss_mitigation_ha = res$fixed$loss_summary$loss_mitigation_ha,
      loss_protection_ha = res$fixed$loss_summary$loss_protection_ha
    )
  })
)

write_csv_bom(
  scenario_summary_tbl,
  file.path(dir_tab, "weight_sensitivity_summary.csv")
)

print(scenario_summary_tbl, n = Inf)

# 8.15 Score-level sensitivity: Spearman correlations --------------------------

sc_names <- weight_scenarios$scenario
n_sc <- length(sc_names)

corr_mat <- matrix(NA_real_, nrow = n_sc, ncol = n_sc,
                   dimnames = list(sc_names, sc_names))

for (i in seq_along(sc_names)) {
  for (j in seq_along(sc_names)) {
    
    vi <- terra::values(scenario_results[[sc_names[i]]]$score$A_dev, mat = FALSE)
    vj <- terra::values(scenario_results[[sc_names[j]]]$score$A_dev, mat = FALSE)
    
    ok <- is.finite(vi) & !is.na(vi) & is.finite(vj) & !is.na(vj)
    
    corr_mat[i, j] <- suppressWarnings(
      stats::cor(vi[ok], vj[ok], method = "spearman")
    )
  }
}

corr_df <- as.data.frame(as.table(corr_mat), stringsAsFactors = FALSE)
colnames(corr_df) <- c("scenario_1", "scenario_2", "rho_spearman")

write_csv_bom(
  corr_df,
  file.path(dir_tab, "weight_sensitivity_spearman_matrix.csv")
)

p_corr <- ggplot(corr_df, aes(x = scenario_1, y = scenario_2, fill = rho_spearman)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho_spearman)), size = 4) +
  scale_fill_viridis_c(limits = c(0, 1), option = "C") +
  coord_equal() +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(title = "Spearman correlation among weight scenarios")

print(p_corr)

ggsave(
  filename = file.path(dir_fig, "weight_sensitivity_spearman_matrix.png"),
  plot     = p_corr,
  width    = 8,
  height   = 6.5,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

# 8.16 Decision-level sensitivity ----------------------------------------------

overlap_sel_tbl <- lapply(sc_names, function(sc) {
  tibble::tibble(
    scenario = sc,
    overlap_selection_pct = selection_overlap_pct(
      scenario_results[[baseline_name]]$fixed$selection_raster,
      scenario_results[[sc]]$fixed$selection_raster
    )
  )
}) %>%
  dplyr::bind_rows()

write_csv_bom(
  overlap_sel_tbl,
  file.path(dir_tab, "weight_sensitivity_selection_overlap_vs_baseline.csv")
)

topk_tbl <- lapply(sc_names, function(sc) {
  tmp <- topk_overlap_pct(
    A_ref = scenario_results[[baseline_name]]$score$A_dev,
    A_cmp = scenario_results[[sc]]$score$A_dev,
    probs = c(0.05, 0.10, 0.20)
  )
  tmp$scenario <- sc
  tmp
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    top_label = dplyr::case_when(
      top_fraction == 0.05 ~ "Top 5%",
      top_fraction == 0.10 ~ "Top 10%",
      top_fraction == 0.20 ~ "Top 20%",
      TRUE ~ paste0("Top ", round(100 * top_fraction), "%")
    )
  )

write_csv_bom(
  topk_tbl,
  file.path(dir_tab, "weight_sensitivity_topk_overlap_vs_baseline.csv")
)

p_topk <- ggplot(
  topk_tbl,
  aes(x = scenario, y = overlap_pct, fill = top_label)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    x = NULL,
    y = "Overlap with reference (%)",
    fill = NULL,
    title = "Stability of top-ranked pixels across weight scenarios"
  )

print(p_topk)

ggsave(
  filename = file.path(dir_fig, "weight_sensitivity_topk_overlap_vs_baseline.png"),
  plot     = p_topk,
  width    = 8.5,
  height   = 5.5,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

# 8.17 Pixel-wise robustness map -----------------------------------------------

sel_stack <- terra::rast(lapply(sc_names, function(sc) scenario_results[[sc]]$fixed$selection_raster))
names(sel_stack) <- sc_names

sel_count <- terra::app(sel_stack, fun = function(x) sum(!is.na(x) & x == 1, na.rm = TRUE))
names(sel_count) <- "selection_count"

terra::writeRaster(
  sel_count,
  file.path(dir_ras, "selection_count_across_weight_scenarios.tif"),
  overwrite = TRUE
)

df_count <- as.data.frame(sel_count, xy = TRUE, na.rm = TRUE)
colnames(df_count)[3] <- "n_sel"

p_count <- ggplot(df_count, aes(x = x, y = y, fill = n_sel)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c(
    breaks = seq_len(length(sc_names)),
    limits = c(1, length(sc_names)),
    name   = "Selected in\nn scenarios"
  ) +
  theme_void() +
  labs(title = "Robustness of spatial prioritization across weight scenarios")

print(p_count)

ggsave(
  filename = file.path(dir_fig, "weight_sensitivity_selection_count_map.png"),
  plot     = p_count,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 300,
  bg       = "white"
)

# 8.18 Iterative loss trajectories across weight scenarios ---------------------

loss_table_all <- dplyr::bind_rows(
  lapply(sc_names, function(sc) {
    scenario_results[[sc]]$iterative$loss_table %>%
      dplyr::mutate(scenario = sc, .before = 1)
  })
)

print(loss_table_all, n = Inf)

write_csv_bom(
  loss_table_all,
  file.path(dir_tab, "loss_table_all_weight_scenarios.csv")
)

loss_table_compare <- loss_table_all %>%
  dplyr::select(
    scenario,
    iteration,
    ha_cum,
    loss_mitigation_ha,
    loss_protection_ha
  ) %>%
  tidyr::pivot_wider(
    names_from  = scenario,
    values_from = c(ha_cum, loss_mitigation_ha, loss_protection_ha),
    names_sep   = "__"
  )

loss_table_compare %>%
  dplyr::select(
    iteration,
    dplyr::starts_with("loss_mitigation"),
    dplyr::starts_with("loss_protection")
  ) %>%
  print(n = Inf, width = Inf)

write_csv_bom(
  loss_table_compare,
  file.path(dir_tab, "loss_table_comparison_wide.csv")
)

loss_plot_df <- loss_table_all %>%
  dplyr::mutate(
    scenario = factor(
      scenario,
      levels = c("reference", "development", "conservation")
    )
  ) %>%
  tidyr::pivot_longer(
    cols = c(loss_mitigation_ha, loss_protection_ha),
    names_to = "var",
    values_to = "loss_ha"
  ) %>%
  dplyr::mutate(
    habitat_class = dplyr::case_when(
      var == "loss_mitigation_ha" ~ "Mitigation",
      var == "loss_protection_ha" ~ "Protection"
    ),
    habitat_class = factor(
      habitat_class,
      levels = c("Mitigation", "Protection")
    )
  )

pal_scen <- c(
  "reference"    = "#1b9e77",
  "development"  = "#d95f02",
  "conservation" = "#66a61e"
)

p_loss_traj <- ggplot(
  loss_plot_df,
  aes(x = ha_cum, y = loss_ha, color = scenario)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_vline(
    xintercept = fp_area_ha,
    linetype   = "dashed",
    linewidth  = 0.9,
    color      = "black"
  ) +
  annotate(
    "text",
    x = fp_area_ha,
    y = Inf,
    label = "Total wind farm\n footprint area",
    vjust = 1.4,
    hjust = -0.1,
    size = 4,
    color = "black"
  ) +
  facet_wrap(~ habitat_class, scales = "free_y", ncol = 1) +
  scale_color_manual(values = pal_scen, name = "Scenario") +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "Cumulative development footprint (ha)",
    y = "Cumulative habitat loss (ha)",
    title = "Cumulative habitat loss trajectories across weight scenarios"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

print(p_loss_traj)

ggsave(
  filename = file.path(dir_fig, "loss_trajectories_weight_scenarios.png"),
  plot     = p_loss_traj,
  width    = 8,
  height   = 10,
  dpi      = 300
)

# ---- 9) ITERATION-BASED SCENARIO -> CLASSIFIED RASTERS + LANDSCAPEMETRICS ----
# Goal: choose one cumulative iteration from one weight scenario, apply this
# cumulative selection as an impact footprint on the same full-landscape
# pre-development habitat raster used in Section 5, and compute the same
# class-level landscape metrics. Then harmonize the resulting forest plot with
# the main forest plot from Section 5 and compare the two.

scenario_use <- "reference"
iter_use <- 8L

if (!(scenario_use %in% names(scenario_results))) {
  stop("scenario_use does not match any available scenario in scenario_results.")
}

res_use <- scenario_results[[scenario_use]]
loss_table_use <- res_use$iterative$loss_table
sel_snapshots_use <- res_use$iterative$sel_snapshots

if (!(iter_use %in% loss_table_use$iteration)) {
  stop("iter_use does not match any available iteration in the selected scenario.")
}

iter_key <- paste0("iter_", iter_use)

if (!iter_key %in% names(sel_snapshots_use)) {
  stop(
    "The requested iteration is not available in sel_snapshots for scenario: ",
    scenario_use, " / ", iter_key
  )
}

r_sel_iter <- terra::ifel(sel_snapshots_use[[iter_key]] == 1L, 1L, NA)
names(r_sel_iter) <- sprintf("selection_%s_to_iter_%02d", scenario_use, iter_use)

# Use the exact same full-landscape baseline as in Section 5
lvl_df_iter <- data.frame(
  value = 0:2,
  class = c("non_habitat", "mitigation_class", "protection_class")
)

r_class_before_iter <- r_class_pre
r_class_before_iter <- terra::as.int(r_class_before_iter)
names(r_class_before_iter) <- "habitat_class_before_iter"
levels(r_class_before_iter) <- list(lvl_df_iter)

# Apply the iteration selection to the full-landscape baseline
r_class_after_iter <- r_class_before_iter
r_class_after_iter[!is.na(r_sel_iter)] <- 0L
r_class_after_iter <- terra::as.int(r_class_after_iter)
names(r_class_after_iter) <- sprintf("habitat_class_after_iter_%02d", iter_use)
levels(r_class_after_iter) <- list(lvl_df_iter)

if (!terra::compareGeom(r_class_before_iter, r_class_after_iter, stopOnError = FALSE)) {
  stop("Geometry mismatch between r_class_before_iter and r_class_after_iter.")
}

out_before_iter <- file.path(
  dir_ras,
  sprintf("habitat_class_before_%s_iter_%02d.tif", scenario_use, iter_use)
)
out_after_iter <- file.path(
  dir_ras,
  sprintf("habitat_class_after_%s_iter_%02d.tif", scenario_use, iter_use)
)
out_sel_iter <- file.path(
  dir_ras,
  sprintf("selection_cumulative_%s_iter_%02d.tif", scenario_use, iter_use)
)

terra::writeRaster(r_class_before_iter, out_before_iter, overwrite = TRUE)
terra::writeRaster(r_class_after_iter,  out_after_iter,  overwrite = TRUE)
terra::writeRaster(r_sel_iter,          out_sel_iter,    overwrite = TRUE)

metrics_before_iter <- calc_metrics_class(r_class_before_iter)
metrics_after_iter  <- calc_metrics_class(r_class_after_iter)

before_df_iter <- do.call(rbind, lapply(metrics_before_iter, norm_lsm, state_label = "before"))
after_df_iter  <- do.call(rbind, lapply(metrics_after_iter,  norm_lsm, state_label = "after"))
all_df_iter    <- rbind(before_df_iter, after_df_iter)

all_df_iter <- merge(all_df_iter, metric_dict,  by = "metric", all.x = TRUE)
all_df_iter <- merge(all_df_iter, metric_units, by = "metric", all.x = TRUE)

if (any(is.na(all_df_iter$category))) {
  stop("Some metrics could not be matched to a category.")
}
if (any(is.na(all_df_iter$unit))) {
  stop("Some metrics could not be matched to a unit.")
}

all_df_iter <- all_df_iter %>%
  dplyr::filter(is.na(class) | class != 0)

wide_iter <- reshape(
  all_df_iter[, c("category", "metric", "unit", "level", "class", "state", "value")],
  timevar   = "state",
  idvar     = c("category", "metric", "unit", "level", "class"),
  direction = "wide"
)

wide_iter <- wide_iter[, c("category", "metric", "unit", "level", "class", "value.before", "value.after")]
wide_iter$delta <- wide_iter$value.after - wide_iter$value.before

cat_order <- c("Amount", "Fragmentation", "Connectivity", "Heterogeneity")
wide_iter$category <- factor(wide_iter$category, levels = cat_order, ordered = TRUE)
wide_iter <- wide_iter[order(wide_iter$category, wide_iter$metric, wide_iter$class), ]

wide_iter_tbl <- tibble::as_tibble(wide_iter)

print(wide_iter_tbl, n = Inf)

write_csv_bom(
  wide_iter_tbl,
  file.path(dir_tab, sprintf("table_all_metrics_%s_iter_%02d.csv", scenario_use, iter_use))
)

forest_df_iter <- wide_iter_tbl %>%
  dplyr::filter(level == "class", class %in% c(1, 2)) %>%
  dplyr::left_join(metric_labels_df, by = "metric") %>%
  dplyr::mutate(
    class = factor(
      as.character(class),
      levels = c("1", "2"),
      labels = c("Mitigation", "Protection")
    ),
    delta_pct = dplyr::if_else(
      is.finite(value.before) & !is.na(value.before) & value.before > 0,
      (value.after - value.before) / value.before * 100,
      NA_real_
    ),
    out_5 = !is.na(delta_pct) & abs(delta_pct) > 5,
    force = abs(delta_pct)
  )

metric_order_iter <- forest_df_iter %>%
  dplyr::filter(!is.na(delta_pct)) %>%
  dplyr::group_by(category, metric_label) %>%
  dplyr::summarise(score_ord = max(abs(delta_pct), na.rm = TRUE), .groups = "drop")

forest_df_iter <- forest_df_iter %>%
  dplyr::left_join(metric_order_iter, by = c("category", "metric_label")) %>%
  dplyr::mutate(
    score_ord = dplyr::coalesce(score_ord, -Inf),
    metric_lab = tidytext::reorder_within(metric_label, -score_ord, category)
  )

p_forest_iter <- ggplot(
  forest_df_iter,
  aes(x = delta_pct, y = metric_lab, color = class)
) +
  annotate(
    "rect",
    xmin = -5, xmax = 5,
    ymin = -Inf, ymax = Inf,
    alpha = 0.30,
    fill  = "grey60"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_point(
    data = dplyr::filter(forest_df_iter, !is.na(delta_pct)),
    aes(alpha = out_5, size = force)
  ) +
  scale_alpha_manual(values = c(`FALSE` = 0.65, `TRUE` = 1), guide = "none") +
  scale_size_continuous(
    range = c(2.5, 10),
    name  = "Change magnitude (%)",
    guide = guide_legend(
      override.aes = list(color = "grey30", alpha = 1, size = 6)
    )
  ) +
  facet_wrap(~ category, scales = "free_y", ncol = 1) +
  tidytext::scale_y_reordered() +
  scale_color_manual(
    values = pal_class,
    name   = "Habitat class",
    drop   = FALSE,
    guide  = guide_legend(
      override.aes = list(size = 6, alpha = 1)
    )
  ) +
  labs(
    x = "Δ (%)",
    y = NULL,
    title = sprintf("%s scenario — iteration %02d vs baseline", scenario_use, iter_use)
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4)
  )

print(p_forest_iter)

ggsave(
  filename = file.path(
    dir_fig,
    sprintf("forest_plot_metrics_%s_iter_%02d.png", scenario_use, iter_use)
  ),
  plot     = p_forest_iter,
  width    = 44,
  height   = 32,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

write_csv_bom(
  forest_df_iter,
  file.path(
    dir_tab,
    sprintf("table_forest_plot_metrics_%s_iter_%02d.csv", scenario_use, iter_use)
  )
)

impact_summary_iter <- loss_table_use %>%
  dplyr::filter(iteration == iter_use)

if (nrow(impact_summary_iter) != 1) {
  stop("Could not extract a unique iteration summary from the selected scenario loss table.")
}

print(impact_summary_iter)

cat("\nScenario used            :", scenario_use, "\n")
cat("Iteration used           :", iter_use, "\n")
cat("Cumulative area (ha)     :", impact_summary_iter$ha_cum, "\n")
cat("Mitigation loss (ha)     :", impact_summary_iter$loss_mitigation_ha, "\n")
cat("Protection loss (ha)     :", impact_summary_iter$loss_protection_ha, "\n")
cat("Outputs written to       :", dir_tab, "\n")

# Harmonize forest plots: axes, point sizes, and metric order

forest_df <- wide_tbl %>%
  dplyr::filter(level == "class", class %in% c(1, 2)) %>%
  dplyr::left_join(metric_labels_df, by = "metric") %>%
  dplyr::mutate(
    class = factor(
      as.character(class),
      levels = c("1", "2"),
      labels = c("Mitigation", "Protection")
    ),
    delta_pct = dplyr::if_else(
      is.finite(value.before) & !is.na(value.before) & value.before > 0,
      (value.after - value.before) / value.before * 100,
      NA_real_
    ),
    out_5 = !is.na(delta_pct) & abs(delta_pct) > 5,
    force = abs(delta_pct)
  )

metric_order <- forest_df %>%
  dplyr::filter(!is.na(delta_pct)) %>%
  dplyr::group_by(category, metric_label) %>%
  dplyr::summarise(score_ord = max(abs(delta_pct), na.rm = TRUE), .groups = "drop")

forest_df <- forest_df %>%
  dplyr::left_join(metric_order, by = c("category", "metric_label")) %>%
  dplyr::mutate(
    score_ord = dplyr::coalesce(score_ord, -Inf),
    metric_lab = tidytext::reorder_within(metric_label, -score_ord, category)
  )

x_vals_common <- c(forest_df$delta_pct, forest_df_iter$delta_pct)
x_vals_common <- x_vals_common[is.finite(x_vals_common)]

x_lim_common <- range(x_vals_common, na.rm = TRUE)
x_lim_common <- c(
  floor(x_lim_common[1] / 5) * 5,
  ceiling(x_lim_common[2] / 5) * 5
)

if (diff(x_lim_common) == 0) {
  x_lim_common <- x_lim_common + c(-5, 5)
}

x_breaks_common <- seq(x_lim_common[1], x_lim_common[2], by = 10)
if (length(x_breaks_common) < 2) {
  x_breaks_common <- c(x_lim_common[1], x_lim_common[2])
}

force_vals_common <- c(forest_df$force, forest_df_iter$force)
force_vals_common <- force_vals_common[is.finite(force_vals_common)]

force_lim_common <- c(0, ceiling(max(force_vals_common, na.rm = TRUE) / 10) * 10)
if (!is.finite(force_lim_common[2]) || force_lim_common[2] == 0) {
  force_lim_common <- c(0, 10)
}

force_breaks_common <- seq(20, force_lim_common[2], by = 20)
if (length(force_breaks_common) == 0) {
  force_breaks_common <- c(0, force_lim_common[2])
}

metric_order_fixed <- gsub("___.*$", "", levels(forest_df$metric_lab))

forest_df_iter$metric_lab <- factor(
  gsub("___.*$", "", as.character(forest_df_iter$metric_lab)),
  levels = metric_order_fixed
)

p_forest_sync <- make_forest_plot(
  df_plot = forest_df,
  plot_title = NULL,
  reordered = TRUE,
  x_lim_common = x_lim_common,
  x_breaks_common = x_breaks_common,
  force_lim_common = force_lim_common,
  force_breaks_common = force_breaks_common,
  pal_class = pal_class
)

p_forest_iter_sync <- make_forest_plot(
  df_plot = forest_df_iter,
  plot_title = sprintf("%s scenario — iteration %02d vs baseline", scenario_use, iter_use),
  reordered = FALSE,
  x_lim_common = x_lim_common,
  x_breaks_common = x_breaks_common,
  force_lim_common = force_lim_common,
  force_breaks_common = force_breaks_common,
  pal_class = pal_class
)

print(p_forest_sync)
print(p_forest_iter_sync)

ggsave(
  filename = file.path(dir_fig, "forest_plot_metrics_synced.png"),
  plot     = p_forest_sync,
  width    = 44,
  height   = 32,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

ggsave(
  filename = file.path(
    dir_fig,
    sprintf("forest_plot_metrics_%s_iter_%02d_synced.png", scenario_use, iter_use)
  ),
  plot     = p_forest_iter_sync,
  width    = 44,
  height   = 32,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

# Comparison of impact magnitude

compare_figs_df <- forest_df %>%
  dplyr::select(category, metric, metric_label, class, delta_pct) %>%
  dplyr::rename(delta_pct_fig1 = delta_pct) %>%
  dplyr::inner_join(
    forest_df_iter %>%
      dplyr::select(category, metric, metric_label, class, delta_pct) %>%
      dplyr::rename(delta_pct_fig2 = delta_pct),
    by = c("category", "metric", "metric_label", "class")
  ) %>%
  dplyr::mutate(
    abs_fig1 = abs(delta_pct_fig1),
    abs_fig2 = abs(delta_pct_fig2),
    gap_abs  = abs_fig1 - abs_fig2,
    interpretation = dplyr::case_when(
      gap_abs > 0  ~ "Figure 2 less impactful",
      gap_abs < 0  ~ "Figure 2 more impactful",
      TRUE         ~ "Same magnitude"
    )
  ) %>%
  dplyr::arrange(category, metric, class)

print(compare_figs_df, n = Inf, width = Inf)

write_csv_bom(
  compare_figs_df,
  file.path(
    dir_fig,
    sprintf("comparison_fig1_vs_fig2_%s_iter_%02d.csv", scenario_use, iter_use)
  )
)

compare_figs_plot_df <- compare_figs_df %>%
  dplyr::mutate(
    metric_lab = factor(metric_label, levels = rev(unique(metric_label)))
  )

p_gap_abs <- ggplot(
  compare_figs_plot_df,
  aes(x = gap_abs, y = metric_lab, color = class, size = abs_fig2)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_point(alpha = 0.9) +
  facet_wrap(~ category, scales = "free_y", ncol = 1) +
  scale_color_manual(
    values = pal_class,
    name   = "Habitat class"
  ) +
  scale_size_continuous(
    name = "Magnitude in Fig. 2 (%)",
    range = c(3, 10)
  ) +
  labs(
    x = "Difference in absolute impact magnitude (%) = |Fig. 1| - |Fig. 2|",
    y = NULL,
    title = sprintf(
      "Impact reduction in Figure 2 relative to Figure 1 (%s, iteration %02d)",
      scenario_use, iter_use
    ),
    subtitle = "Positive values indicate that Figure 2 is less impactful"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

print(p_gap_abs)

ggsave(
  filename = file.path(
    dir_fig,
    sprintf("gap_abs_fig1_vs_fig2_%s_iter_%02d.png", scenario_use, iter_use)
  ),
  plot     = p_gap_abs,
  width    = 34,
  height   = 26,
  units    = "cm",
  dpi      = 300,
  bg       = "white"
)

compare_plot_df <- compare_figs_df %>%
  dplyr::mutate(
    change_dir = dplyr::case_when(
      gap_abs > 0  ~ "Reduced in Figure 2",
      gap_abs < 0  ~ "Increased in Figure 2",
      TRUE         ~ "No change"
    ),
    metric_label = factor(metric_label, levels = rev(unique(metric_label)))
  )

pal_change <- c(
  "Reduced in Figure 2"   = "#2E8B57",
  "Increased in Figure 2" = "#B22222",
  "No change"             = "grey70"
)

p_compare <- ggplot(
  compare_plot_df,
  aes(x = metric_label, y = gap_abs, fill = change_dir)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "black") +
  geom_col(width = 0.72) +
  coord_flip() +
  facet_grid(. ~ class, scales = "fixed", space = "fixed") +
  scale_fill_manual(values = pal_change) +
  labs(
    title = "Differences in impact between actual and simulated",
    x = NULL,
    y = expression(Delta*" in landscape metrics")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(colour = "#4d4d4d"),
    panel.spacing.x = grid::unit(1.1, "lines")
  )

print(p_compare)

ggsave(
  filename = file.path(dir_fig, "fig_compare_scenarios.png"),
  plot     = p_compare,
  bg       = "white",
  width    = 10,
  height   = 6,
  dpi      = 600
)

# ---- 10) RESTART AFTER CRASH -------------------------------------------------
# This section assumes that Sections 0 and 1 have already been executed in the
# current R session. Those sections load all required packages, define project
# paths, and initialize utility functions used throughout the workflow.
#
# After a crash or a fresh R session, the correct restart order is:
#
#   0) LIBRARIES
#   1) PATHS AND GLOBAL UTILITIES
#   10) RESTART AFTER CRASH
#
# Executing this section without first running Sections 0–1 will result in
# missing objects (paths, functions, or packages) and the restart procedure
# will fail.
#
# This section reloads the main raster and vector objects written to disk after
# Section 4, allowing the workflow to resume from Section 5 onward without
# rerunning the computationally expensive steps (interpolation and SDM
# reconstruction).

# 10.1 Paths -------------------------------------------------------------------

# Recreate parameters needed to rebuild file paths
footprint_buffer_m <- 10

# Reporting thresholds from Tremblay et al. (2025)
mitigation_threshold <- 0.103
protection_threshold <- 0.194

# Automatically find the injected pre-development SDM written in Section 3
predev_full_files <- list.files(
  dir_ras,
  pattern = "^SDM_predev_injected_North_IDW_beta_.*_n_.*\\.tif$",
  full.names = TRUE
)

if (length(predev_full_files) == 0) {
  stop("No injected pre-development SDM was found in: ", dir_ras)
}
if (length(predev_full_files) > 1) {
  stop(
    "Multiple injected pre-development SDM files were found in: ", dir_ras,
    "\nPlease keep only one file or manually specify which one to reload."
  )
}

predev_full_path <- predev_full_files[1]

class_pre_path  <- file.path(dir_ras, "habitat_class_2thresholds_predev.tif")
class_post_path <- file.path(dir_ras, "habitat_class_2thresholds_postdev_impacted.tif")
footprint_gpkg  <- file.path(
  dir_vec,
  sprintf("wind_footprint_total_reconstructed_buffer%sm.gpkg", footprint_buffer_m)
)
footprint_raster_path <- file.path(
  dir_ras,
  sprintf("wind_footprint_total_reconstructed_buffer%sm_raster.tif", footprint_buffer_m)
)

predev_full_masked_path <- file.path(dir_ras, "SDM_predev_injected_North_masked_water.tif")
pred_masked_path        <- file.path(dir_ras, "prediction_2024_beaupre_masked_water.tif")
water_mask_path         <- file.path(dir_ras, "water_mask_SeigneurieBeaupre_COTQ2022.tif")

# 10.2 Checks ------------------------------------------------------------------

if (!file.exists(predev_full_path)) {
  stop("Missing injected pre-development SDM: ", predev_full_path)
}
if (!file.exists(predev_full_masked_path)) {
  stop("Missing masked pre-development SDM: ", predev_full_masked_path)
}
if (!file.exists(pred_masked_path)) {
  stop("Missing masked post-development SDM: ", pred_masked_path)
}
if (!file.exists(class_pre_path)) {
  stop("Missing pre-development classified raster: ", class_pre_path)
}
if (!file.exists(class_post_path)) {
  stop("Missing post-development classified raster: ", class_post_path)
}
if (!file.exists(footprint_gpkg)) {
  stop("Missing total wind-footprint geopackage: ", footprint_gpkg)
}
if (!file.exists(footprint_raster_path)) {
  stop("Missing total wind-footprint raster: ", footprint_raster_path)
}
if (!file.exists(water_mask_path)) {
  stop("Missing water mask raster: ", water_mask_path)
}

# 10.3 Reload core rasters -----------------------------------------------------

r_pred           <- rast(pred_masked_path)
r_predev_full    <- rast(predev_full_masked_path)
r_water          <- rast(water_mask_path)
r_class_pre      <- rast(class_pre_path)
r_class_post_imp <- rast(class_post_path)
r_footprint_all  <- rast(footprint_raster_path)

names(r_pred)           <- "prediction_2024_beaupre"
names(r_predev_full)    <- "predev_SDM"
names(r_water)          <- "water_mask"
names(r_class_pre)      <- "habitat_class_predev"
names(r_class_post_imp) <- "habitat_class_postdev_impacted"
names(r_footprint_all)  <- "wind_footprint_all"

# 10.4 Restore raster levels ---------------------------------------------------

lvl_df_raster <- data.frame(
  value = 0:2,
  class = c("non_habitat", "mitigation_class", "protection_class")
)

levels(r_class_pre)      <- list(lvl_df_raster)
levels(r_class_post_imp) <- list(lvl_df_raster)

# 10.5 Reload footprint --------------------------------------------------------

footprint_all_sf <- st_read(footprint_gpkg, quiet = TRUE) |>
  st_make_valid()

footprint_all_v <- vect(footprint_all_sf)

# 10.6 Rebuild objects used by later sections ----------------------------------

r_before <- r_class_pre
names(r_before) <- "habitat_predev"

r_after <- r_class_post_imp
names(r_after) <- "habitat_postdev"

metric_dict <- data.frame(
  metric = c("ca","np","area_mn","lpi","ed","ai","tca","cohesion","mesh","enn_mn","area_cv","shape_cv"),
  category = c(
    "Amount",
    "Fragmentation",
    "Fragmentation",
    "Fragmentation",
    "Fragmentation",
    "Fragmentation",
    "Fragmentation",
    "Connectivity",
    "Connectivity",
    "Connectivity",
    "Heterogeneity",
    "Heterogeneity"
  ),
  stringsAsFactors = FALSE
)

metric_units <- data.frame(
  metric = c("ca","np","area_mn","lpi","ed","ai","tca","cohesion","mesh","enn_mn","area_cv","shape_cv"),
  unit = c("ha","count","ha","%","m/ha","%","ha","%","ha","m","-","-"),
  stringsAsFactors = FALSE
)

metric_labels <- c(
  ca        = "CA (ha)",
  np        = "NP (count)",
  area_mn   = "AREA_MN (ha)",
  lpi       = "LPI (%)",
  ed        = "ED (m/ha)",
  ai        = "AI (%)",
  tca       = "TCA (ha)",
  cohesion  = "COHESION (%)",
  mesh      = "MESH (ha)",
  enn_mn    = "ENN_MN (m)",
  area_cv   = "AREA_CV (-)",
  shape_cv  = "SHAPE_CV (-)"
)

metric_labels_df <- tibble(
  metric = names(metric_labels),
  metric_label = unname(metric_labels)
)

pal_class <- c(
  "Mitigation" = "#e8ac31",
  "Protection" = "#a86466"
)

# 10.7 Informative message -----------------------------------------------------

cat("\nRestart objects successfully reloaded.\n")
cat("Injected SDM reloaded from:\n", predev_full_masked_path, "\n", sep = "")
cat("Masked post-development SDM reloaded from:\n", pred_masked_path, "\n", sep = "")
cat("Water mask reloaded from:\n", water_mask_path, "\n", sep = "")
cat("Footprint raster reloaded from:\n", footprint_raster_path, "\n", sep = "")
cat("You can now resume from Section 5.\n")
# ---- END ----