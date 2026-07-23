################################################################################
# BayesianTrapData.R
#
# Builds the input arrays for the Bayesian occupancy model from the combined,
# all-trap-types dataset produced by summarize_data_all_traptypes.R.
#
# FILTERS
#   1. Records from 2018 onward (year >= 2018).
#   2. Sites whose coordinate falls within the study extent
#      (Zones/study_extent.shp), via point-in-polygon.
#
# OUTPUTS (exactly three files, all written to Data/Processed/)
#   * pa_array.rds        - 3-D binary presence/absence array,
#                           dims [site, trap_type, year].
#                           1 = EGC caught in that site/trap/year; 0 otherwise.
#                           Note: a 0 means either "sampled, none caught" OR
#                           "not sampled" - use traptype_array.rds to tell them
#                           apart.
#   * cpue_fukui.rds      - 2-D CPUE matrix, dims [site, year], FUKUI TRAPS ONLY
#                           (total Fukui EGC caught / total Fukui traps set).
#                           Site-years with no Fukui effort are NA.
#   * traptype_array.rds  - 3-D binary trap-type covariate array,
#                           dims [site, trap_type, year].
#                           1 = that trap type was deployed at that site in that
#                           year (effort > 0); 0 = not deployed.
#
# All three arrays share the SAME site, trap_type and year dimension labels and
# ordering, so pa[i, t, y] and traptype[i, t, y] refer to the same cell, and
# cpue_fukui[i, y] uses the same site/year indices.
#
# INPUT: Data/data_all_sources_all_traptypes.csv
#
# Requires the `sf` package for the spatial filter.
################################################################################

library(tidyverse)
library(sf)

# ---------------------------------------------------------------------------
# Locate the repo root so paths work regardless of the working directory.
# This script lives in Code/, so the repo root is its parent folder.
# ---------------------------------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p)))
  }
  NULL
}
script_dir <- tryCatch(get_script_dir(), error = function(e) NULL)
if (!is.null(script_dir)) setwd(normalizePath(file.path(script_dir, "..")))

if (!file.exists("Data/data_all_sources_all_traptypes.csv")) {
  stop("Could not find Data/data_all_sources_all_traptypes.csv relative to '",
       getwd(), "'. Run summarize_data_all_traptypes.R first and/or set the ",
       "working directory to the repo root.")
}

outdir <- "Data/Processed"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load combined data
# ---------------------------------------------------------------------------
dat <- read.csv("Data/data_all_sources_all_traptypes.csv",
                stringsAsFactors = FALSE)
# expected columns:
# site_name, year, trap_type, data_sources, catch, effort, cpue, latitude, longitude

# ---------------------------------------------------------------------------
# 1) Filter to 2018 onward
# ---------------------------------------------------------------------------
dat <- dat %>% filter(!is.na(year), year >= 2018)

# ---------------------------------------------------------------------------
# 2) Filter to sites within the study extent (point-in-polygon)
# ---------------------------------------------------------------------------
extent <- st_read("Zones/study_extent.shp", quiet = TRUE)
extent <- st_make_valid(extent)

# one coordinate per site (already unique per site in the input, but guard
# against NA coordinates, which cannot be placed and are therefore dropped)
site_pts <- dat %>%
  distinct(site_name, latitude, longitude) %>%
  filter(!is.na(latitude), !is.na(longitude))

site_sf <- st_as_sf(site_pts, coords = c("longitude", "latitude"),
                    crs = 4326)
site_sf <- st_transform(site_sf, st_crs(extent))

inside <- lengths(st_intersects(site_sf, extent)) > 0
sites_in_extent <- site_pts$site_name[inside]

dropped <- setdiff(unique(dat$site_name), sites_in_extent)
message(length(sites_in_extent), " sites within study extent; ",
        length(dropped), " sites dropped (outside extent or missing coords).")

dat <- dat %>% filter(site_name %in% sites_in_extent)

# ---------------------------------------------------------------------------
# Shared dimension labels. Every array below is indexed by these, in this
# order, so the three outputs line up cell for cell.
# ---------------------------------------------------------------------------
sites <- sort(unique(dat$site_name))
traps <- sort(unique(dat$trap_type))
years <- sort(unique(dat$year))

nsite <- length(sites)
ntrap <- length(traps)
nyear <- length(years)

dnames <- list(site      = sites,
               trap_type = traps,
               year      = as.character(years))

# row/col/slice index of every record, reused by both 3-D arrays
i_site <- match(dat$site_name, sites)
i_trap <- match(dat$trap_type, traps)
i_year <- match(dat$year,      years)
idx    <- cbind(i_site, i_trap, i_year)

# ---------------------------------------------------------------------------
# 3-D binary presence/absence array [site, trap_type, year]
# 1 = at least one EGC caught. pmax() guards against duplicate rows for the
# same site/trap/year (a 1 must never be overwritten by a 0).
# ---------------------------------------------------------------------------
present <- as.integer(!is.na(dat$catch) & dat$catch > 0)

pa <- array(0L, dim = c(nsite, ntrap, nyear), dimnames = dnames)
pa[idx] <- pmax(pa[idx], present)

saveRDS(pa, file.path(outdir, "pa_array.rds"))

# ---------------------------------------------------------------------------
# 3-D binary trap-type covariate array [site, trap_type, year]
# 1 = that trap type was actually deployed at that site in that year.
# ---------------------------------------------------------------------------
deployed <- as.integer(!is.na(dat$effort) & dat$effort > 0)

traptype <- array(0L, dim = c(nsite, ntrap, nyear), dimnames = dnames)
traptype[idx] <- pmax(traptype[idx], deployed)

saveRDS(traptype, file.path(outdir, "traptype_array.rds"))

# ---------------------------------------------------------------------------
# 2-D CPUE matrix [site, year] - FUKUI TRAPS ONLY.
# Rows/cols keep the full site and year vectors above, so this matrix is
# index-compatible with the two 3-D arrays. Site-years with no Fukui effort
# stay NA.
# ---------------------------------------------------------------------------
fukui <- dat %>%
  filter(trap_type == "Fukui") %>%
  group_by(site_name, year) %>%
  summarise(catch  = sum(catch,  na.rm = TRUE),
            effort = sum(effort, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(cpue = ifelse(effort > 0, catch / effort, NA_real_))

cpue_fukui <- matrix(NA_real_, nrow = nsite, ncol = nyear,
                     dimnames = list(site = sites, year = as.character(years)))
cpue_fukui[cbind(match(fukui$site_name, sites),
                 match(fukui$year,      years))] <- fukui$cpue

saveRDS(cpue_fukui, file.path(outdir, "cpue_fukui.rds"))

# ---------------------------------------------------------------------------
message("Done. ", nsite, " sites x ", ntrap, " trap types x ", nyear,
        " years (", min(years), "-", max(years), ").")
message("  pa_array.rds       [site, trap_type, year] - ", sum(pa), " detections")
message("  traptype_array.rds [site, trap_type, year] - ", sum(traptype),
        " sampled site/trap/years")
message("  cpue_fukui.rds     [site, year]            - ",
        sum(!is.na(cpue_fukui)), " site-years with Fukui effort")
