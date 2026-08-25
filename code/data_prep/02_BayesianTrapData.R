################################################################################
# 02_BayesianTrapData.R
#
# Build replicate-level inputs for the dynamic occupancy model directly from
# raw trap records. The array layout is [site, year, replicate]. A replicate is
# one actual trap deployment/pull, not a trap-type summary.
#
# Filters are applied before array construction:
#   * year >= 2018
#   * summer (June-September), matching the existing workflow
#   * Fukui, Shrimp, and Minnow traps only (including named variants)
#   * sites within data/SpatialData/study_extent.shp
#
# WDFW, Makah, DNWR, and DFO have trap-level observations and are included.
# WSG and Padilla are supplied only as aggregated effort/catch totals, and the
# Drayton objects do not identify which traps generated aggregated catches;
# these sources cannot be converted honestly to replicate-level binary data and
# are therefore excluded from this product.
#
# Outputs in data/model_data (original filenames retained):
#   PresenceArrayBinary.rds
#       integer [site, year, replicate]: 1 detected, 0 not detected, NA unused
#   TrapTypeArraysNew.rds
#       named list of Fukui/Shrimp/Minnow one-hot integer arrays with the same
#       dimensions and NA pattern as PresenceArrayBinary.rds
#   cpue_fukui.rds
#       Fukui catch per trap deployment, [site, year]
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
here::i_am("code/data_prep/02_BayesianTrapData.R")

FIRST_YEAR <- 2018L
SUMMER_MONTHS <- 6:9
CANONICAL_TRAPS <- c("Fukui", "Shrimp", "Minnow")

raw_dir <- here::here("data", "raw")
extent_file <- here::here("data", "SpatialData", "study_extent.shp")
outdir <- here::here("data", "model_data")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

canonical_trap_type <- function(x) {
  cleaned <- tolower(trimws(as.character(x)))
  cleaned <- gsub("[[:punct:]]+", " ", cleaned)
  cleaned <- gsub("[[:space:]]+", " ", cleaned)
  result <- rep(NA_character_, length(cleaned))
  result[grepl("\\bfukui\\b", cleaned)] <- "Fukui"
  result[is.na(result) & grepl("\\bshrimp\\b", cleaned)] <- "Shrimp"
  result[is.na(result) & grepl("\\bminnow\\b", cleaned)] <- "Minnow"
  result
}

number <- function(x) suppressWarnings(as.numeric(as.character(x)))

standardize <- function(source, site, year, month, trap_type, catch,
                        latitude, longitude, order_key) {
  data.frame(
    source = source,
    site_name = trimws(as.character(site)),
    year = suppressWarnings(as.integer(as.character(year))),
    month = suppressWarnings(as.integer(as.character(month))),
    trap_type = canonical_trap_type(trap_type),
    catch = number(catch),
    latitude = number(latitude),
    longitude = number(longitude),
    order_key = as.character(order_key),
    stringsAsFactors = FALSE
  ) %>%
    filter(
      nzchar(site_name), !is.na(year), year >= FIRST_YEAR,
      month %in% SUMMER_MONTHS, !is.na(trap_type),
      is.finite(catch), catch >= 0
    ) %>%
    mutate(
      site = paste0(source, "::", site_name),
      presence = as.integer(catch > 0)
    )
}

# -----------------------------------------------------------------------------
# Makah: one row is one retrieved trap.
# -----------------------------------------------------------------------------
makah <- read.csv(file.path(raw_dir, "Makah_catchdata.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
makah_reps <- standardize(
  "Makah", makah[["Location"]], makah[["Year"]],
  match(makah[["Month"]], month.name), makah[["Trap Type"]],
  makah[["Count of EGC"]], makah[["Latitude"]], makah[["Longitude"]],
  paste(makah[["Date Deployed"]], makah[["Trap Number"]], seq_len(nrow(makah)))
)

# -----------------------------------------------------------------------------
# WDFW: both files contain one row per trap deployment.
# -----------------------------------------------------------------------------
clean_wdfw_site <- function(x) {
  x <- trimws(as.character(x))
  replacements <- c(
    "Duckland & 105" = "Duckland", "South" = "South",
    "Henry Island Preserve" = "Henry Island Preserve",
    "Ocean Shores - Airport" = "Ocean Shores",
    "Ocean Shores Airport" = "Ocean Shores", "Post Point Lagoon" = "Post Point Lagoon",
    "West Samish" = "Samish River", "North" = "Crandall Spit",
    "Tokeland Hotel" = "Tokeland"
  )
  hit <- match(x, names(replacements))
  x[!is.na(hit)] <- unname(replacements[hit[!is.na(hit)]])
  x
}

wdfw_old <- read.csv(
  file.path(raw_dir, "WDFW", "WDFW.EGC_2018-2022_Effort_Final.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
wdfw_old_date <- as.Date(wdfw_old[["Date_Deployed"]], "%m/%d/%Y")
wdfw_old_keep <- trimws(wdfw_old[["Waterbody"]]) != "Drayton Harbor"
wdfw_old_reps <- standardize(
  "WDFW", clean_wdfw_site(wdfw_old[["Site_Name"]]),
  format(wdfw_old_date, "%Y"), format(wdfw_old_date, "%m"),
  wdfw_old[["Trap_Type"]], wdfw_old[["CAMA_Total"]],
  wdfw_old[["Latitude"]], wdfw_old[["Longitude"]],
  paste(wdfw_old[["Date_Deployed"]], wdfw_old[["Trap_Number"]],
        seq_len(nrow(wdfw_old)))
) %>% filter(wdfw_old_keep[match(order_key, paste(
  wdfw_old[["Date_Deployed"]], wdfw_old[["Trap_Number"]], seq_len(nrow(wdfw_old))
))])

wdfw_new <- read.csv(
  file.path(raw_dir, "WDFW", "WDFW_EGC_Data_Collection_PUBLIC_09092025_pg1.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
wdfw_new_date <- as.Date(wdfw_new[["Set Date and Time"]], "%m/%d/%Y")
wdfw_new_reps <- standardize(
  "WDFW", clean_wdfw_site(wdfw_new[["Site Name"]]),
  format(wdfw_new_date, "%Y"), format(wdfw_new_date, "%m"),
  wdfw_new[["Trap Type"]], wdfw_new[["Total Catch"]],
  wdfw_new[["y"]], wdfw_new[["x"]],
  paste(wdfw_new[["Set Date and Time"]], wdfw_new[["Trap Number"]],
        seq_len(nrow(wdfw_new)))
)

# -----------------------------------------------------------------------------
# DNWR: effort rows define replicates; catches are counted by TrapIDJoin.
# -----------------------------------------------------------------------------
read_dnwr_year <- function(year) {
  effort <- read.csv(file.path(raw_dir, "DNWR", paste0("DNWR_", year, "_effort.csv")),
                     stringsAsFactors = FALSE, check.names = FALSE)
  catches <- read.csv(file.path(raw_dir, "DNWR", paste0("DNWR_", year, "_catch.csv")),
                      stringsAsFactors = FALSE, check.names = FALSE)
  # The 2018 catch workbook has dozens of trailing unnamed spreadsheet columns;
  # repair only those names so dplyr can operate while preserving useful labels.
  blank_names <- is.na(names(catches)) | !nzchar(names(catches))
  names(catches)[blank_names] <- paste0("unused_", seq_len(sum(blank_names)))
  names(catches) <- make.unique(names(catches))
  crab_col <- "Crab #"
  catches <- catches[!is.na(catches[[crab_col]]) & nzchar(trimws(catches[[crab_col]])), ]
  catch_counts <- catches %>%
    group_by(TrapIDJoin) %>%
    summarise(catch = n(), .groups = "drop")
  effort <- left_join(effort, catch_counts, by = "TrapIDJoin")
  effort$catch[is.na(effort$catch)] <- 0
  pull_date <- as.Date(as.character(effort[["Trap Pull Date"]]), "%Y%m%d")
  lon <- number(effort[["Longitude"]])
  site <- case_when(
    lon < -123.18 ~ "Dungeness Spit Base",
    lon > -123.16 & lon < -123.13 ~ "Graveyard Spit West (Channel)",
    lon > -123.13 ~ "Graveyard Spit East (Lagoon)",
    TRUE ~ NA_character_
  )
  standardize(
    "DNWR", site, year, format(pull_date, "%m"), effort[["Trap Type"]],
    effort$catch, effort[["Latitude"]], effort[["Longitude"]],
    paste(effort[["TrapIDJoin"]], seq_len(nrow(effort)))
  )
}
dnwr_reps <- bind_rows(lapply(2018:2020, read_dnwr_year))

# -----------------------------------------------------------------------------
# DFO: raw files can contain multiple species rows for the same physical trap.
# Collapse only those duplicate species records to one trap, then flag XMB.
# GearCode 51 is the documented Fukui code; other codes cannot be assigned to
# the three requested trap types and are excluded.
# -----------------------------------------------------------------------------
dfo_early <- read.csv(file.path(raw_dir, "DFO", "allegcdata_2007_2021.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
dfo_late <- read.csv(file.path(raw_dir, "DFO", "allegcdata_2022_2024.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
# Base rbind intentionally performs the harmless coercion needed for GearCode,
# which is character in the early file and numeric in the later file.
dfo <- rbind(dfo_early, dfo_late) %>%
  filter(TrapUsability %in% c(0, 15), as.character(GearCode) == "51") %>%
  mutate(
    trap_id = paste(Source, SetNum, Year, Month, Day, GeogLoc, TrapNum, sep = "::"),
    is_egc = as.character(Species) == "XMB"
  ) %>%
  group_by(trap_id) %>%
  summarise(
    site_name = first(GeogLoc), year = first(Year), month = first(Month),
    trap_type = "Fukui", catch = sum(is_egc, na.rm = TRUE),
    latitude = first(startLAT), longitude = first(startLONG),
    order_key = first(trap_id), .groups = "drop"
  )
dfo_reps <- standardize(
  "DFO", dfo$site_name, dfo$year, dfo$month, dfo$trap_type, dfo$catch,
  dfo$latitude, dfo$longitude, dfo$order_key
)

# -----------------------------------------------------------------------------
# Combine, apply spatial filter, and assign consecutive replicate slots.
# -----------------------------------------------------------------------------
reps <- bind_rows(wdfw_old_reps, wdfw_new_reps, makah_reps, dnwr_reps, dfo_reps)
if (nrow(reps) == 0L) stop("No replicate records remain after source filters.")

site_points <- reps %>%
  group_by(site) %>%
  summarise(
    latitude = if (all(!is.finite(latitude))) NA_real_ else mean(latitude[is.finite(latitude)]),
    longitude = if (all(!is.finite(longitude))) NA_real_ else mean(longitude[is.finite(longitude)]),
    .groups = "drop"
  ) %>%
  filter(is.finite(latitude), is.finite(longitude))

extent <- st_make_valid(st_read(extent_file, quiet = TRUE))
site_sf <- st_as_sf(site_points, coords = c("longitude", "latitude"),
                    crs = 4326, remove = FALSE) %>%
  st_transform(st_crs(extent))
sites_in_extent <- site_sf$site[lengths(st_intersects(site_sf, extent)) > 0L]

reps <- reps %>%
  filter(site %in% sites_in_extent) %>%
  arrange(site, year, order_key, trap_type) %>%
  group_by(site, year) %>%
  mutate(replicate = row_number()) %>%
  ungroup()
if (nrow(reps) == 0L) stop("No replicate records remain after spatial filtering.")

sites <- sort(unique(reps$site))
years <- sort(unique(reps$year))
nrep <- max(reps$replicate)
dims <- c(length(sites), length(years), nrep)
dnames <- list(site = sites, year = as.character(years),
               replicate = as.character(seq_len(nrep)))
idx <- cbind(match(reps$site, sites), match(reps$year, years), reps$replicate)

fill_integer <- function(values) {
  out <- array(NA_integer_, dim = dims, dimnames = dnames)
  out[idx] <- as.integer(values)
  out
}

presence <- fill_integer(reps$presence)
trap_indicators <- lapply(
  setNames(CANONICAL_TRAPS, CANONICAL_TRAPS),
  function(tt) fill_integer(reps$trap_type == tt)
)

make_cpue <- function(data) {
  pooled <- data %>%
    group_by(site, year) %>%
    summarise(cpue = sum(catch) / n(), .groups = "drop")
  out <- matrix(NA_real_, nrow = length(sites), ncol = length(years),
                dimnames = list(site = sites, year = as.character(years)))
  out[cbind(match(pooled$site, sites), match(pooled$year, years))] <- pooled$cpue
  out
}
cpue_fukui <- make_cpue(filter(reps, trap_type == "Fukui"))

# -----------------------------------------------------------------------------
# Validation and serialization.
# -----------------------------------------------------------------------------
used <- !is.na(presence)
n_used <- apply(used, c(1, 2), sum)
stopifnot(
  identical(names(dimnames(presence)), c("site", "year", "replicate")),
  all(presence[used] %in% c(0L, 1L)),
  all(apply(used, c(1, 2), function(z) all(diff(as.integer(z)) <= 0))),
  all(vapply(trap_indicators, function(a) identical(!is.na(a), used), logical(1))),
  all(Reduce(`+`, trap_indicators)[used] == 1L),
  identical(dimnames(presence)[c("site", "year")], dimnames(cpue_fukui)),
  all(unname(!is.na(cpue_fukui)) ==
        unname(apply(trap_indicators$Fukui == 1L, c(1, 2), any, na.rm = TRUE))),
  all(reps$year >= FIRST_YEAR),
  all(reps$trap_type %in% CANONICAL_TRAPS)
)

outputs <- list(
  PresenceArrayBinary = presence,
  TrapTypeArraysNew = trap_indicators,
  cpue_fukui = cpue_fukui
)
for (nm in names(outputs)) saveRDS(outputs[[nm]], file.path(outdir, paste0(nm, ".rds")))

message("Saved [site, year, replicate] = ", paste(dims, collapse = " x "),
        " (", min(years), "-", max(years), ").")
message("  ", nrow(reps), " actual trap replicates; ", sum(presence, na.rm = TRUE),
        " detection-positive replicates")
message("  trap types: ", paste(names(table(reps$trap_type)), table(reps$trap_type),
                                sep = "=", collapse = ", "))
message("  sources: ", paste(names(table(reps$source)), table(reps$source),
                             sep = "=", collapse = ", "))
