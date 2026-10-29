#!/usr/bin/env Rscript
# =============================================================================
# 34_marine_distances.R - Least-cost marine distances for the Kivalliq
#
# Produces:
#   1. 02_info/connectivity_covariates.csv  -> input for 35_variance_partitioning.R
#   2. 02_info/hydro_dist_pops.csv          -> fills the Hydro_km column that has
#                                              been NA in Table_IBD_pairwise.csv
#
# Does NOT use marmap / gdistance / lc.dist. Reads the ETOPO GeoTIFF with terra,
# builds an 8-neighbour graph over marine cells with great-circle edge weights,
# and runs Dijkstra with igraph. No PROJ reprojection is involved.
#
# Two things learned the hard way, both handled below:
#   - A 1 arc-minute grid does not resolve narrow inlets, so the marine mask
#     breaks into ~140 disconnected pools. Snapping a site to the NEAREST
#     marine cell can drop it into an isolated pool and every distance becomes
#     Inf. Sites are therefore snapped to the nearest cell ON THE LARGEST
#     CONNECTED COMPONENT.
#   - Anchor placement matters more than it looks. An anchor a few km off the
#     open channel lands on land or in a pool. Snap distances are reported so
#     the placement can be checked.
#
# HOR (Baker Lake) is freshwater and cannot be routed over marine cells. It is
# handled as two explicit segments: a marine leg from the mouth of Chesterfield
# Inlet, plus a fixed freshwater leg. Both are reported separately so the
# manuscript can distinguish isolation by distance from isolation by the
# freshwater transition.
#
# Input: ETOPO 2022 bathymetry at 1 arc-minute, cropped to the study area and
#   saved as 02_info/bathy/kivalliq_etopo.tif (other location: BATHY_TIF=<path>)
#
# Run from the pipeline root (BASE_DIR):
#   $RSCRIPT 01_scripts/R/34_marine_distances.R
# =============================================================================

if (nzchar(Sys.getenv("RLIB"))) .libPaths(c(Sys.getenv("RLIB"), .libPaths()))
suppressPackageStartupMessages({
  library(terra); library(igraph); library(data.table)
})

BASE <- Sys.getenv("BASE_DIR", getwd())
TIF  <- Sys.getenv("BATHY_TIF", file.path(BASE, "02_info/bathy/kivalliq_etopo.tif"))
INFO <- file.path(BASE, "02_info")

# NA = every marine cell is passable. A depth cap (e.g. -300) is biologically
# defensible for a coastal migrant, but on this grid it severs passages and
# fragments the graph further. Left off; revisit only with a finer raster.
MAX_DEPTH <- NA

# HOR: marine leg from the mouth of Chesterfield Inlet, then freshwater upstream.
HOR_MOUTH_LON <- -90.70
HOR_MOUTH_LAT <-  63.35
HOR_RIVER_KM  <- 300

MAX_SNAP_WARN <- 30   # km; flag sites that land far from their true position

# ---------------------------------------------------------------------------
sites <- data.table(
  pop = c("AKL","AUL","CRB","DIA","MEL",
          "ITI","KGJ","NOP","PAM","SUP","TIN","WHI","HOR"),
  region = c(rep("Rankin",5), rep("Naujaat",7), "Baker"),
  lat = c(62.837, 62.812, 62.467, 62.833, 62.867,
          66.312, 66.459, 66.633, 67.083, 66.383, 66.383, 65.967, 64.723),
  lon = c(-91.307, -91.370, -92.333, -92.444, -92.117,
          -84.465, -85.254, -86.883, -87.056, -86.717, -84.440, -85.043, -97.934))

# Anchors: the two plausible entry routes for the Atlantic lineage into the
# region. Placed in open channel, not against a shoreline.
anchors <- data.table(
  name = c("HudsonStrait", "FoxeBasin"),
  lat  = c(61.50, 68.00),
  lon  = c(-78.50, -80.50))

# ---------------------------------------------------------------------------
cat("Reading bathymetry...\n")
r <- rast(TIF)
print(ext(r))
m  <- as.matrix(r, wide = TRUE)
nr <- nrow(m); nc <- ncol(m)

water <- m < 0
if (!is.na(MAX_DEPTH)) water <- water & (m >= MAX_DEPTH)
n_water <- sum(water, na.rm = TRUE)
cat(sprintf("Marine cells: %d of %d (%.1f%%)\n",
            n_water, nr * nc, 100 * n_water / (nr * nc)))
if (n_water < 1000) stop("Too few marine cells - check MAX_DEPTH or the raster")

id <- matrix(NA_integer_, nr, nc)
id[which(water)] <- seq_len(n_water)

lat_row <- yFromRow(r, 1:nr)
lon_col <- xFromCol(r, 1:nc)
dlat_km <- abs(yres(r)) * 110.574
dlon_km <- abs(xres(r)) * 111.320 * cos(lat_row * pi / 180)   # varies by row

# --- 8-neighbour edge list -------------------------------------------------
cat("Building graph...\n")
add_offset <- function(dr, dc, wfun) {
  r1 <- max(1, 1 - dr):min(nr, nr - dr)
  c1 <- max(1, 1 - dc):min(nc, nc - dc)
  a <- id[r1, c1, drop = FALSE]
  b <- id[r1 + dr, c1 + dc, drop = FALSE]
  ok <- !is.na(a) & !is.na(b)
  if (!any(ok)) return(NULL)
  rows <- row(a)[ok]
  data.table(from = a[ok], to = b[ok], w = wfun(r1[rows]))
}
edges <- list(
  add_offset(0,  1, function(i) dlon_km[i]),                       # E
  add_offset(1,  0, function(i) rep(dlat_km, length(i))),          # S
  add_offset(1,  1, function(i) sqrt(dlon_km[i]^2 + dlat_km^2)),   # SE
  add_offset(1, -1, function(i) sqrt(dlon_km[i]^2 + dlat_km^2)))   # SW
E <- rbindlist(Filter(Negate(is.null), edges))
cat(sprintf("Edges: %d\n", nrow(E)))

g <- graph_from_data_frame(E[, .(from, to)], directed = FALSE,
                           vertices = data.frame(name = seq_len(n_water)))
E(g)$weight <- E$w

comp <- components(g)
big_id <- which.max(comp$csize)
cat(sprintf("Connected components: %d (largest holds %d cells, %.1f%%)\n",
            comp$no, comp$csize[big_id], 100 * comp$csize[big_id] / n_water))

# --- Snap to the nearest cell ON THE LARGEST COMPONENT ---------------------
wnodes <- id[!is.na(id)]
wr <- row(id)[!is.na(id)]
wc <- col(id)[!is.na(id)]
wlat <- lat_row[wr]
wlon <- lon_col[wc]

on_big <- comp$membership[wnodes] == big_id
bnodes <- wnodes[on_big]; blat <- wlat[on_big]; blon <- wlon[on_big]
cat(sprintf("Cells available for snapping: %d\n", length(bnodes)))

snap <- function(lon, lat) {
  d <- sqrt(((blon - lon) * cos(lat * pi / 180) * 111.32)^2 +
            ((blat - lat) * 110.574)^2)
  i <- which.min(d)
  list(node = bnodes[i], km = d[i], lon = blon[i], lat = blat[i])
}

pts <- rbind(
  sites[pop != "HOR", .(name = pop, lon, lat)],
  data.table(name = "HOR_mouth", lon = HOR_MOUTH_LON, lat = HOR_MOUTH_LAT),
  anchors[, .(name, lon, lat)])

sn <- rbindlist(lapply(seq_len(nrow(pts)), function(i) {
  s <- snap(pts$lon[i], pts$lat[i])
  data.table(name = pts$name[i], node = s$node, snap_km = s$km,
             snap_lon = round(s$lon, 3), snap_lat = round(s$lat, 3))
}))
cat("\nSnap distance to the connected marine network (km):\n")
print(sn[order(-snap_km)])

far <- sn[snap_km > MAX_SNAP_WARN]
if (nrow(far) > 0) {
  cat(sprintf("\nNOTE: %d point(s) snapped more than %d km. A 1 arc-minute grid\n",
              nrow(far), MAX_SNAP_WARN))
  cat("does not resolve narrow inlets, so upstream sites are displaced to the\n")
  cat("open coast. Disclose this in Methods; it inflates short distances more\n")
  cat("than long ones, which is conservative for a gradient test.\n")
}

if (any(duplicated(sn$node)))
  cat("\nNOTE: some sites snapped to the SAME cell - their distances are identical:\n  ",
      paste(sn$name[sn$node %in% sn$node[duplicated(sn$node)]], collapse = ", "), "\n")

# --- Dijkstra --------------------------------------------------------------
cat("\nComputing least-cost distances...\n")
D <- distances(g, v = as.character(sn$node), to = as.character(sn$node),
               algorithm = "dijkstra")
dimnames(D) <- list(sn$name, sn$name)

if (any(!is.finite(D)))
  stop("Infinite distances remain despite snapping to the largest component - ",
       "inspect the anchor placement.")

out <- sites[, .(pop, region, lat, lon)]
out[, `:=`(d_HS = NA_real_, d_FH = NA_real_, freshwater_km = 0)]

for (i in seq_len(nrow(out))) {
  p <- out$pop[i]
  src <- if (p == "HOR") "HOR_mouth" else p
  extra <- if (p == "HOR") HOR_RIVER_KM else 0
  out$d_HS[i] <- D[src, "HudsonStrait"] + extra
  out$d_FH[i] <- D[src, "FoxeBasin"]    + extra
  out$freshwater_km[i] <- extra
}
out[, d_min  := pmin(d_HS, d_FH)]
out[, d_diff := d_HS - d_FH]

# Column names expected by 33_connectivity_gradient.R
out[, marine_corridor_km := d_HS]
out[, coast_km := freshwater_km]

fwrite(out, file.path(INFO, "connectivity_covariates.csv"))
cat("\n=== connectivity_covariates.csv ===\n")
print(out[, .(pop, region, lat, d_HS = round(d_HS), d_FH = round(d_FH),
              d_min = round(d_min), d_diff = round(d_diff))])

# --- Collinearity: the whole argument depends on this ----------------------
safe_cor <- function(dt) {
  x <- as.matrix(dt)
  keep <- apply(x, 2, function(v) is.finite(sd(v)) && sd(v) > 0)
  if (sum(keep) < 2) return(NULL)
  round(cor(x[, keep, drop = FALSE]), 3)
}
cat("\nCorrelations among candidate predictors (all 13 populations):\n")
print(safe_cor(out[, .(lat, d_HS, d_FH, d_min, d_diff)]))
cat("\nMarine populations only (HOR excluded):\n")
print(safe_cor(out[pop != "HOR", .(lat, d_HS, d_FH, d_min, d_diff)]))

r_lat <- suppressWarnings(cor(out$lat, out$d_HS, use = "complete.obs"))
if (is.finite(r_lat)) {
  cat(sprintf("\ncor(latitude, corridor distance) = %.3f\n", r_lat))
  if (abs(r_lat) > 0.8)
    cat("WARNING: latitude and corridor distance are collinear. With 13\n",
        "populations the two hypotheses cannot be separated - report this\n",
        "rather than picking whichever fits better.\n", sep = "")
  else
    cat("The two predictors are separable - the gradient test is informative.\n")
}

cat("\nWhich anchor is closer, by region:\n")
print(out[, .(mean_d_HS = round(mean(d_HS)), mean_d_FH = round(mean(d_FH)),
              mean_diff = round(mean(d_diff))), by = region])

# --- Pairwise marine distances, for the IBD Hydro_km column ---------------
marine_pops <- setdiff(sites$pop, "HOR")
pop_names <- c(marine_pops, "HOR_mouth")
M <- D[pop_names, pop_names]
rownames(M) <- colnames(M) <- c(marine_pops, "HOR")
M["HOR", ] <- M["HOR", ] + HOR_RIVER_KM
M[, "HOR"] <- M[, "HOR"] + HOR_RIVER_KM
diag(M) <- 0
write.csv(round(M, 2), file.path(INFO, "hydro_dist_pops.csv"))

cat("\nPairwise marine distances -> hydro_dist_pops.csv\n")
cat(sprintf("Range: %.0f - %.0f km\n", min(M[M > 0]), max(M)))
cat("\nNote: the isolation-by-distance analysis of the manuscript uses the\n")
cat("distances of compute_hydro_distances.R (09_ibd.sh), not this matrix.\n")