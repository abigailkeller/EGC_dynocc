################################################################################
# BayesianDataSitesAdded.R
#
# Extends the replicate-level model products from 02_BayesianTrapData.R with
# sites that 02 deliberately excludes because they only have season-pooled
# catch/effort totals, never individual trap-level results: WSG (Washington
# Sea Grant), Padilla Bay, and Drayton Harbor (tagged "NWStraits" in
# data/data_all_sources_all_traptypes.csv). Same filters as 02: year >= 2018,
# summer only (already baked into data_all_sources_all_traptypes.csv upstream
# in 01_summarize_data_all_traptypes.R), canonical traps only (Fukui/Shrimp/
# Minnow), and sites within data/SpatialData/study_extent.shp. This script
# only ADDS sites to 02's data -- it does not change 02 or anything
# downstream of it (the Bayesian model code itself is untouched).
#
# We don't have site-specific (per-trap) results for these sites, only a
# season total catch and a season total trap-effort count per
# (site, year, trap type). To keep the same [site, year, replicate] array
# shape as the real data, each site-year-trap-type's total effort becomes
# that many replicate slots -- e.g. an effort of 24 becomes 24 replicate
# entries -- and every one of those slots is set to the SAME value: whether
# the pooled catch for that site-year-trap-type was greater than zero. This
# is an explicit imputation (we cannot know which specific deployment(s)
# caught anything), not a real replicate-level detection history. Treat
# occupancy/detection estimates that depend heavily on these added sites
# with that caveat in mind.
#
# Outputs in data/model_data (same names as 02_BayesianTrapData.R, with
# "_xtraSites" appended) -- same shape and NA conventions as the originals,
# just with the added sites' rows filled in:
#   PresenceArrayBinary_xtraSites.rds
#       integer [site, year, replicate]: 1 detected, 0 not detected, NA
#       unused.
#   TrapTypeArraysNew_xtraSites.rds
#       named list of Fukui/Shrimp/Minnow one-hot integer arrays, same
#       dimensions and NA pattern as PresenceArrayBinary_xtraSites.rds.
#   cpue_fukui_xtraSites.rds
#       Fukui catch per trap deployment, [site, year] -- for the added
#       sites this is the season-pooled Fukui catch / season-pooled Fukui
#       effort (imputed), i.e. exactly the same quantity 02 computes for the
#       real sites, just not decomposable to the trap level.
#   cpue_xtraSites.rds
#       Catch per trap deployment pooled across all three trap types, same
#       [site, year] shape and coverage as PresenceArrayBinary_xtraSites.rds
#       -- the counterpart to 02's cpue.rds.
################################################################################

project_library <- file.path(getwd(), ".r-library")
if (dir.exists(project_library)) {
  .libPaths(c(normalizePath(project_library), .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
})

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it with install.packages('here').")
}
here::i_am("code/data_prep/BayesianDataSitesAdded.R")

FIRST_YEAR <- 2018L  # same cutoff 02_BayesianTrapData.R uses
CANONICAL_SOURCES <- c("WDFW", "Makah", "DNWR", "DFO")
CANONICAL_TRAPS <- c("Fukui", "Shrimp", "Minnow")

outdir <- here::here("data", "model_data")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Rebuild the real, trap-level baseline exactly as 02_BayesianTrapData.R
#    does, without touching its on-disk outputs, so this script always
#    starts from the current real data rather than a possibly-stale .rds.
# -----------------------------------------------------------------------------
prep_file <- here::here("code", "data_prep", "02_BayesianTrapData.R")
prep_code <- readLines(prep_file, warn = FALSE)
prep_code <- sub(
  'outdir <- here::here\\("data", "model_data"\\)',
  'outdir <- file.path(tempdir(), "model_data_baseline")',
  prep_code
)
baseline_env <- new.env(parent = globalenv())
eval(parse(text = prep_code), envir = baseline_env)

presence_base  <- baseline_env$presence
traps_base     <- baseline_env$trap_indicators
cpue_fukui_base <- baseline_env$cpue_fukui
cpue_base      <- baseline_env$cpue
base_sites     <- dimnames(presence_base)$site
base_years     <- dimnames(presence_base)$year
base_nrep      <- dim(presence_base)[3]

# -----------------------------------------------------------------------------
# 2. Identify the aggregate-only, in-study-extent sites and pull their
#    season-pooled (site, year, trap_type) catch/effort.
# -----------------------------------------------------------------------------
all_sources <- read.csv(
  here::here("data", "data_all_sources_all_traptypes.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    data_sources = strsplit(as.character(data_sources), ",\\s*")
  )

row_has_canonical_source <- vapply(all_sources$data_sources, function(s) {
  any(trimws(s) %in% CANONICAL_SOURCES)
}, logical(1))

# Exclude a site entirely if ANY of its rows (any year, any trap type) ever
# cite a canonical trap-level source -- not just the rows being kept here.
# A site with even one WDFW/Makah/DNWR/DFO record already belongs to the
# real trap-level product (used, or filtered out there for its own reasons),
# and shouldn't also show up as an "aggregate-only" addition.
sites_with_canonical_source <- unique(all_sources$site_name[row_has_canonical_source])

aggregate_only <- all_sources %>%
  filter(!site_name %in% sites_with_canonical_source,
         trap_type %in% CANONICAL_TRAPS,
         is.finite(catch), is.finite(effort), effort > 0,
         is.finite(year), year >= FIRST_YEAR)

site_coords <- aggregate_only %>%
  group_by(site_name) %>%
  summarise(
    latitude = mean(latitude, na.rm = TRUE),
    longitude = mean(longitude, na.rm = TRUE),
    source_tag = paste(sort(unique(trimws(unlist(data_sources)))), collapse = "+"),
    .groups = "drop"
  ) %>%
  filter(is.finite(latitude), is.finite(longitude))

extent <- st_make_valid(st_read(here::here("data", "SpatialData", "study_extent.shp"), quiet = TRUE))
site_sf <- st_as_sf(site_coords, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(st_crs(extent))
in_extent_names <- site_sf$site_name[lengths(st_intersects(site_sf, extent)) > 0L]

xtra <- aggregate_only %>%
  filter(site_name %in% in_extent_names) %>%
  left_join(select(site_coords, site_name, source_tag), by = "site_name") %>%
  mutate(
    site = paste0(source_tag, "::", site_name),
    year = as.character(year)
  ) %>%
  group_by(site, site_name, source_tag, year, trap_type) %>%
  summarise(catch = sum(catch), effort = sum(effort), .groups = "drop") %>%
  mutate(
    effort = as.integer(round(effort)),
    presence = as.integer(catch > 0)
  )

message("Aggregate-only, in-extent sites added: ", length(unique(xtra$site)),
        " (", nrow(xtra), " site-year-trap_type season totals, ",
        sum(xtra$effort), " synthetic replicate slots)")

# -----------------------------------------------------------------------------
# 3. Expand the site/year/replicate dimensions to fit both the real baseline
#    and the added sites, then fill in the added sites' synthetic replicates.
# -----------------------------------------------------------------------------
xtra_sites <- sort(unique(xtra$site))
xtra_years <- sort(unique(xtra$year))

all_sites <- sort(unique(c(base_sites, xtra_sites)))
all_years <- sort(unique(c(base_years, xtra_years)))

xtra_max_rep <- xtra %>%
  group_by(site, year) %>%
  summarise(total_effort = sum(effort), .groups = "drop") %>%
  pull(total_effort) %>%
  max(0L)
nrep <- max(base_nrep, xtra_max_rep)

dims <- c(length(all_sites), length(all_years), nrep)
dnames <- list(site = all_sites, year = all_years, replicate = as.character(seq_len(nrep)))

presence <- array(NA_integer_, dim = dims, dimnames = dnames)
trap_indicators <- setNames(
  lapply(CANONICAL_TRAPS, function(tt) array(NA_integer_, dim = dims, dimnames = dnames)),
  CANONICAL_TRAPS
)

# Copy the real baseline in first.
presence[base_sites, base_years, seq_len(base_nrep)] <- presence_base
for (tt in CANONICAL_TRAPS) {
  trap_indicators[[tt]][base_sites, base_years, seq_len(base_nrep)] <- traps_base[[tt]]
}

# Fill in the added sites: each (site, year) gets its trap-type blocks laid
# out back-to-back in `effort` slots apiece; every slot in a block shares
# that trap-type's presence/absence call, since the source data has no
# finer resolution than the season total.
xtra_by_site_year <- xtra %>% arrange(site, year, trap_type)
for (i in seq_len(nrow(xtra_by_site_year))) {
  row <- xtra_by_site_year[i, ]
  s <- row$site; y <- row$year; tt <- row$trap_type; n <- row$effort
  filled_so_far <- sum(!is.na(presence[s, y, ]))
  slots <- seq.int(filled_so_far + 1L, filled_so_far + n)
  presence[s, y, slots] <- row$presence
  for (other in CANONICAL_TRAPS) {
    trap_indicators[[other]][s, y, slots] <- as.integer(other == tt)
  }
}

# -----------------------------------------------------------------------------
# 4. cpue_fukui: same [site, year] shape, Fukui catch / Fukui effort, pooled
#    at the season level for the added sites exactly like 02 already does
#    for the real ones.
# -----------------------------------------------------------------------------
cpue_fukui <- matrix(NA_real_, nrow = length(all_sites), ncol = length(all_years),
                     dimnames = list(site = all_sites, year = all_years))
cpue_fukui[base_sites, base_years] <- cpue_fukui_base

xtra_fukui <- xtra %>%
  filter(trap_type == "Fukui") %>%
  group_by(site, year) %>%
  summarise(cpue = sum(catch) / sum(effort), .groups = "drop")
cpue_fukui[cbind(match(xtra_fukui$site, all_sites), match(xtra_fukui$year, all_years))] <- xtra_fukui$cpue

# -----------------------------------------------------------------------------
# 5. cpue: same [site, year] shape, pooled across all three trap types --
#    the counterpart to 02's cpue.rds. Coverage matches presence/absence
#    exactly, same as cpue_fukui matches trap_indicators$Fukui.
# -----------------------------------------------------------------------------
cpue <- matrix(NA_real_, nrow = length(all_sites), ncol = length(all_years),
              dimnames = list(site = all_sites, year = all_years))
cpue[base_sites, base_years] <- cpue_base

xtra_pooled <- xtra %>%
  group_by(site, year) %>%
  summarise(cpue = sum(catch) / sum(effort), .groups = "drop")
cpue[cbind(match(xtra_pooled$site, all_sites), match(xtra_pooled$year, all_years))] <- xtra_pooled$cpue

# -----------------------------------------------------------------------------
# Validation -- same shape/contract as 02_BayesianTrapData.R's checks, plus a
# check that the real baseline sites/years pass through unchanged.
# -----------------------------------------------------------------------------
used <- !is.na(presence)
site_year_used <- apply(used, c(1, 2), any)
stopifnot(
  identical(names(dimnames(presence)), c("site", "year", "replicate")),
  all(presence[used] %in% c(0L, 1L)),
  all(apply(used, c(1, 2), function(z) all(diff(as.integer(z)) <= 0))),
  all(vapply(trap_indicators, function(a) identical(!is.na(a), used), logical(1))),
  all(Reduce(`+`, trap_indicators)[used] == 1L),
  identical(dimnames(presence)[c("site", "year")], dimnames(cpue_fukui)),
  all(unname(!is.na(cpue_fukui)) ==
        unname(apply(trap_indicators$Fukui == 1L, c(1, 2), any, na.rm = TRUE))),
  identical(dimnames(presence)[c("site", "year")], dimnames(cpue)),
  identical(unname(!is.na(cpue)), unname(site_year_used)),
  identical(presence_base, presence[base_sites, base_years, seq_len(base_nrep)]),
  identical(cpue_fukui_base, cpue_fukui[base_sites, base_years]),
  identical(cpue_base, cpue[base_sites, base_years])
)

outputs <- list(
  PresenceArrayBinary_xtraSites = presence,
  TrapTypeArraysNew_xtraSites = trap_indicators,
  cpue_fukui_xtraSites = cpue_fukui,
  cpue_xtraSites = cpue
)
for (nm in names(outputs)) saveRDS(outputs[[nm]], file.path(outdir, paste0(nm, ".rds")))

message("Saved [site, year, replicate] = ", paste(dims, collapse = " x "),
        " (", min(all_years), "-", max(all_years), ").")
message("  ", length(base_sites), " real sites + ", length(xtra_sites),
        " added (aggregate-only, in-extent) sites = ", length(all_sites), " total.")

# -----------------------------------------------------------------------------
# 6. site_zone_map_xtraSites.csv -- same [site, zone_id, zone] shape as
#    site_zone_map.csv, extended to the added sites via the identical
#    spatial-join logic 03_map_model_sites_to_zones.R uses (nearest-zone
#    fallback for points that land just outside a hand-drawn polygon).
#    Needed because main_bm.R indexes connectivity by zone, and the added
#    sites need a zone just like the real ones do.
# -----------------------------------------------------------------------------
base_zone_map <- read.csv(here::here("data", "model_data", "site_zone_map.csv"),
                          stringsAsFactors = FALSE)
stopifnot(identical(sort(base_zone_map$site), sort(base_sites)))

xtra_site_points <- site_coords %>%
  mutate(site = paste0(source_tag, "::", site_name)) %>%
  filter(site %in% xtra_sites) %>%
  select(site, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

zones_shp <- st_read(here::here("data", "SpatialData", "zones.shp"), quiet = TRUE) %>%
  st_make_valid()
xtra_site_points <- st_transform(xtra_site_points, st_crs(zones_shp))
zone_hits <- st_intersects(xtra_site_points, zones_shp)

if (any(lengths(zone_hits) > 1L)) {
  multi <- which(lengths(zone_hits) > 1L)
  for (i in multi) {
    candidates <- zone_hits[[i]]
    boundary_distance <- as.numeric(st_distance(
      xtra_site_points[i, ], st_boundary(zones_shp[candidates, ])
    ))
    zone_hits[[i]] <- candidates[which.max(boundary_distance)]
  }
}
outside <- which(lengths(zone_hits) == 0L)
if (length(outside)) {
  nearest <- st_nearest_feature(xtra_site_points[outside, ], zones_shp)
  distances_m <- as.numeric(st_distance(
    xtra_site_points[outside, ], zones_shp[nearest, ], by_element = TRUE
  ))
  near_boundary <- distances_m <= 1000
  zone_hits[outside[near_boundary]] <- lapply(nearest[near_boundary], function(i) i)
  if (any(!near_boundary)) {
    message("site_zone_map_xtraSites: outside all zones and >1km from the nearest one: ",
            paste(xtra_site_points$site[outside[!near_boundary]], collapse = "; "))
  }
}
zone_index <- vapply(zone_hits, function(x) if (length(x) == 1L) x else NA_integer_, integer(1))

xtra_zone_map <- data.frame(
  site = xtra_site_points$site,
  zone_id = zones_shp$zone_id[zone_index],
  zone = zones_shp$name[zone_index]
)

zone_map_xtraSites <- bind_rows(base_zone_map, xtra_zone_map) %>%
  arrange(match(site, all_sites))

stopifnot(
  identical(zone_map_xtraSites$site, all_sites),
  !anyDuplicated(zone_map_xtraSites$site),
  all(!is.na(zone_map_xtraSites$zone_id) == !is.na(zone_map_xtraSites$zone))
)

write.csv(zone_map_xtraSites, file.path(outdir, "site_zone_map_xtraSites.csv"), row.names = FALSE)
message("Saved ", nrow(zone_map_xtraSites), " site-to-zone mappings (",
        sum(is.na(zone_map_xtraSites$zone_id)), " unmapped) to site_zone_map_xtraSites.csv")
