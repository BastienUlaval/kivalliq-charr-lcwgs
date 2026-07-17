#!/usr/bin/env Rscript
# =============================================================================
# compute_hydro_distances.R
# In-water (least-cost) connectivity distances between sampling sites, marmap.
#
# All populations treated as ANADROMOUS (per local / Inuit knowledge),
# including the lake populations HOR, PAM, NOP. For a lake fish the marmap
# marine distance alone is not enough: it misses the channel/inlet the fish
# swims between the lake and open sea. Each population is therefore modelled by
#   - a MARINE ENTRY POINT (river/inlet mouth on open sea), and
#   - an ACCESS distance (km) along the channel/inlet from the site to that
#     entry point (0 for coastal / near-shore river sites).
# Total distance(i,j) = marine_lc(entry_i, entry_j) + access_i + access_j.
#
# Files (all in 14_ibd/):
#   bathy_kivalliq.csv   <- bathymetry grid (cached after first online run)
#   snap_QC.csv          <- entry point, snap displacement, access km per pop
#   hydro_distances.csv  <- Pop1, Pop2, Hydro_km  (input for ibd_analysis.R)
#
# RUN (env with marmap):
#   conda activate marmap
#   Rscript 01_scripts/R/compute_hydro_distances.R <base_dir>
#     base_dir = pipeline root (script writes/reads under <base_dir>/14_ibd/)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript compute_hydro_distances.R <base_dir>")
}
BASE_DIR <- args[1]

suppressPackageStartupMessages({
  library(marmap)
  library(dplyr)
})

# --- Tunables ---------------------------------------------------------------
RES_ARCMIN   <- 0.5
SNAP_WARN_KM <- 10

# --- Paths (derived from base_dir argument) ----------------------------------
IBD_DIR    <- file.path(BASE_DIR, "14_ibd")
BATHY_FILE <- file.path(IBD_DIR, "bathy_kivalliq.csv")
QC_FILE    <- file.path(IBD_DIR, "snap_QC.csv")
OUT_FILE   <- file.path(IBD_DIR, "hydro_distances.csv")
dir.create(IBD_DIR, showWarnings = FALSE, recursive = TRUE)

# --- 1. Sampling coordinates (all 13 pops; all anadromous) ------------------
coords <- data.frame(
  Pop = c("AKL","AUL","CRB","DIA","MEL","HOR","ITI","KGJ","NOP","PAM","TIN","SUP","WHI"),
  Lat = c(62.8372049, 62.905473, 62.466667, 62.833333, 62.866667, 64.723358,
          66.311737, 66.459357, 66.633333, 67.083333, 66.3826439, 66.383333, 65.96663),
  Lon = c(-91.3071148, -91.15114147, -92.333333, -92.44444, -92.116667, -97.933631,
          -84.4652392, -85.254424, -86.883333, -87.055555, -84.4398478, -86.716667, -85.0426395),
  stringsAsFactors = FALSE
)

# --- 2. Sites needing a MANUAL marine entry + access channel length ---------
# Fill once in QGIS / Google Earth, measuring ALONG the water course:
#   Lon_mouth, Lat_mouth = where the drainage reaches OPEN SEA (NOT deep inside
#                          a narrow inlet, which marmap cannot route through)
#   access_km            = channel/inlet length from the site to that open-sea
#                          entry point, following the water
# HOR: Baker Lake -> Chesterfield Inlet -> Hudson Bay. Put the mouth near the
#      ENTRANCE of Chesterfield Inlet into Hudson Bay; access_km then covers the
#      entire Baker Lake -> inlet-entrance run (long: inlet is ~180 km).
# You can also add any RIVER site here if its QC snap looks wrong.
manual_entry <- data.frame(
  Pop       = c("HOR", "PAM", "NOP"),
  Lon_mouth = c(-90.5847847, -86.4051877, -86.7424741),  # HOR, PAM, NOP
  Lat_mouth = c( 63.3731790,  66.5252368,  66.5333919),  # HOR, PAM, NOP
  access_km = c(402, 64.24, 13.58),              # HOR, PAM, NOP (measured in Google Earth)
  stringsAsFactors = FALSE
)

# --- helper: Haversine (km) -------------------------------------------------
hav_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  dlat <- (lat2 - lat1) * pi/180; dlon <- (lon2 - lon1) * pi/180
  a <- sin(dlat/2)^2 + cos(lat1*pi/180)*cos(lat2*pi/180)*sin(dlon/2)^2
  2 * R * asin(sqrt(a))
}

# --- 3. Build marine entry point + access km for every population -----------
coords$entry_lon <- coords$Lon
coords$entry_lat <- coords$Lat
coords$access_km <- 0
for (lp in manual_entry$Pop) {
  me  <- manual_entry[manual_entry$Pop == lp, ]
  idx <- which(coords$Pop == lp)
  if (any(is.na(c(me$Lon_mouth, me$Lat_mouth, me$access_km)))) {
    stop("Manual entry not filled for ", lp,
         " -- set Lon_mouth/Lat_mouth/access_km in 'manual_entry'.")
  }
  coords$entry_lon[idx] <- me$Lon_mouth
  coords$entry_lat[idx] <- me$Lat_mouth
  coords$access_km[idx] <- me$access_km
}

# --- 4. Bathymetry: from disk if present, else download ---------------------
if (file.exists(BATHY_FILE)) {
  cat("Reading cached bathymetry:", BATHY_FILE, "\n")
  bathy <- read.bathy(BATHY_FILE, header = TRUE, sep = ",")
} else {
  cat("No cached bathymetry; downloading from NOAA (res =", RES_ARCMIN, ") ...\n")
  bathy <- getNOAA.bathy(lon1 = -99, lon2 = -82, lat1 = 61, lat2 = 68,
                         resolution = RES_ARCMIN)
  write.csv(as.xyz(bathy), BATHY_FILE, row.names = FALSE)
  cat("Saved bathymetry to:", BATHY_FILE, "\n")
}

# --- 5. Snap every entry point to the nearest SEA cell (depth < 0) ----------
xyz <- as.xyz(bathy); names(xyz) <- c("lon", "lat", "depth")
sea <- xyz[is.finite(xyz$depth) & xyz$depth < 0, ]
snap_one <- function(lon, lat) {
  d2 <- ((sea$lon - lon) * cos(lat * pi/180))^2 + (sea$lat - lat)^2
  k  <- which.min(d2)
  c(lon = sea$lon[k], lat = sea$lat[k])
}
snapped <- t(mapply(snap_one, coords$entry_lon, coords$entry_lat))
coords$snap_lon <- snapped[, "lon"]
coords$snap_lat <- snapped[, "lat"]
coords$snap_km  <- round(mapply(hav_km, coords$entry_lat, coords$entry_lon,
                                coords$snap_lat, coords$snap_lon), 2)

qc <- coords[, c("Pop","Lat","Lon","entry_lat","entry_lon",
                 "snap_lat","snap_lon","snap_km","access_km")]
write.csv(qc, QC_FILE, row.names = FALSE)
cat("\n--- QC: entry snap displacement + access channel (km) ---\n")
print(qc[order(-qc$snap_km), c("Pop","snap_km","access_km")], row.names = FALSE)
flagged <- qc$Pop[qc$snap_km > SNAP_WARN_KM]
if (length(flagged) > 0)
  cat("\nWARNING: entry snapped >", SNAP_WARN_KM, "km -> verify on a map:\n  ",
      paste(flagged, collapse = ", "), "\n")

# --- 6. Marine least-cost distances between snapped entry points ------------
trans  <- trans.mat(bathy, min.depth = 0, max.depth = NULL)
loc    <- coords[, c("snap_lon", "snap_lat")]
lc     <- lc.dist(trans, loc, res = "dist")
lc_mat <- as.matrix(lc)
rownames(lc_mat) <- colnames(lc_mat) <- coords$Pop
if (any(!is.finite(lc_mat)))
  cat("\nWARNING: some marine distances non-finite; inspect snap_QC.csv\n")

# --- 7. Total distance = marine + access_i + access_j -----------------------
acc  <- setNames(coords$access_km, coords$Pop)
out  <- data.frame()
pops <- coords$Pop
for (i in seq_along(pops)) {
  for (j in seq_along(pops)) {
    if (i < j) {
      total <- lc_mat[pops[i], pops[j]] + acc[pops[i]] + acc[pops[j]]
      out <- rbind(out, data.frame(
        Pop1 = pops[i], Pop2 = pops[j],
        Hydro_km = round(as.numeric(total), 2),
        stringsAsFactors = FALSE
      ))
    }
  }
}
write.csv(out, OUT_FILE, row.names = FALSE)
cat("\nWrote", OUT_FILE, "with", nrow(out), "pairs.\n")
cat("Distance range:", round(min(out$Hydro_km), 1), "-",
    round(max(out$Hydro_km), 1), "km\n")