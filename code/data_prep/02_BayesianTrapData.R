################################################################################
# BayesianTrapData.R
#
# Build three aligned inputs for the Bayesian occupancy model from the combined
# monitoring dataset produced by summarize_data_all_traptypes.R.
#
# FILTERS
#   * 2018 onward (change FIRST_YEAR below if a literal post-2018 cutoff is
#     required).
#   * Shrimp, Fukui, and Minnow traps only. Spelling, case, punctuation, and
#     qualified variants (for example "Shrimp - C" and "Minnow - Unmodified")
#     are consolidated to the three canonical names.
#   * Sites inside data/SpatialData/study_extent.shp.
#
# OUTPUTS (data/model_data/)
#   * PresenceArrayBinary.rds: binary [site, trap_type, year] detections.
#   * TrapTypeArraysNew.rds:   binary [site, trap_type, year] deployment.
#   * cpue_fukui.rds:          numeric [site, year] Fukui CPUE.
#
# All products use exactly the same site and year labels in exactly the same
# order. The script validates this invariant before and after writing the files.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
})

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it with install.packages('here').")
}
here::i_am("code/data_prep/02_BayesianTrapData.R")

FIRST_YEAR <- 2018L
CANONICAL_TRAPS <- c("Shrimp", "Fukui", "Minnow")

canonical_trap_type <- function(x) {
  cleaned <- tolower(trimws(as.character(x)))
  cleaned <- gsub("[[:punct:]]+", " ", cleaned)
  cleaned <- gsub("[[:space:]]+", " ", cleaned)
  result <- rep(NA_character_, length(cleaned))
  result[grepl("\\bshrimp\\b", cleaned)] <- "Shrimp"
  result[is.na(result) & grepl("\\bfukui\\b", cleaned)] <- "Fukui"
  result[is.na(result) & grepl("\\bminnow\\b", cleaned)] <- "Minnow"
  result
}

assert_model_products_aligned <- function(presence, traptype, cpue,
                                          canonical_traps = CANONICAL_TRAPS) {
  stopifnot(
    is.array(presence), length(dim(presence)) == 3L,
    is.array(traptype), length(dim(traptype)) == 3L,
    is.matrix(cpue),
    identical(dim(presence), dim(traptype)),
    identical(dimnames(presence), dimnames(traptype)),
    identical(names(dimnames(presence)), c("site", "trap_type", "year")),
    identical(names(dimnames(cpue)), c("site", "year")),
    identical(dimnames(presence)$site, dimnames(cpue)$site),
    identical(dimnames(presence)$year, dimnames(cpue)$year),
    identical(dimnames(presence)$trap_type, canonical_traps),
    !anyDuplicated(dimnames(presence)$site),
    !anyDuplicated(dimnames(presence)$year),
    all(presence %in% c(0L, 1L)),
    all(traptype %in% c(0L, 1L)),
    all(presence <= traptype),
    all(is.na(cpue) | is.finite(cpue)),
    all(is.na(cpue) | cpue >= 0)
  )

  fukui_deployed <- traptype[, "Fukui", , drop = TRUE] == 1L
  stopifnot(
    identical(dim(fukui_deployed), dim(cpue)),
    identical(dimnames(fukui_deployed), dimnames(cpue)),
    identical(!is.na(cpue), fukui_deployed)
  )
  invisible(TRUE)
}

input_file <- here::here("data", "data_all_sources_all_traptypes.csv")
extent_file <- here::here("data", "SpatialData", "study_extent.shp")
outdir <- here::here("data", "model_data")

if (!file.exists(input_file)) {
  stop("Missing ", input_file, ". Run summarize_data_all_traptypes.R first.")
}
if (!file.exists(extent_file)) {
  stop("Missing ", extent_file, ".")
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c(
  "site_name", "year", "trap_type", "catch", "effort",
  "latitude", "longitude"
)
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0L) {
  stop("Input is missing required columns: ",
       paste(missing_columns, collapse = ", "))
}

dat <- dat %>%
  transmute(
    site_name = trimws(as.character(site_name)),
    year = suppressWarnings(as.integer(as.character(year))),
    trap_type = canonical_trap_type(trap_type),
    catch = suppressWarnings(as.numeric(catch)),
    effort = suppressWarnings(as.numeric(effort)),
    latitude = suppressWarnings(as.numeric(latitude)),
    longitude = suppressWarnings(as.numeric(longitude))
  ) %>%
  filter(
    nzchar(site_name),
    !is.na(year), year >= FIRST_YEAR,
    !is.na(trap_type)
  )

if (nrow(dat) == 0L) {
  stop("No records remain after the year and trap-type filters.")
}
if (any(dat$catch < 0, na.rm = TRUE) || any(dat$effort < 0, na.rm = TRUE)) {
  stop("Catch and effort must be non-negative.")
}

# Use one representative coordinate per normalized site. A site is retained
# when that representative point intersects the study extent.
site_points <- dat %>%
  group_by(site_name) %>%
  summarise(
    latitude = if (all(is.na(latitude))) NA_real_ else mean(latitude, na.rm = TRUE),
    longitude = if (all(is.na(longitude))) NA_real_ else mean(longitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(latitude), is.finite(longitude))

extent <- st_make_valid(st_read(extent_file, quiet = TRUE))
site_sf <- st_as_sf(
  site_points,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)
site_sf <- st_transform(site_sf, st_crs(extent))
inside <- lengths(st_intersects(site_sf, extent)) > 0L
sites_in_extent <- site_sf$site_name[inside]

message(
  length(sites_in_extent), " sites inside the study extent; ",
  length(setdiff(unique(dat$site_name), sites_in_extent)),
  " outside or missing coordinates."
)

# Consolidate all source and spelling variants to one row per canonical
# site/trap/year. Rows with effort but missing catch are treated as zero catch;
# rows without positive effort are not considered sampled.
cell_data <- dat %>%
  filter(site_name %in% sites_in_extent) %>%
  group_by(site_name, trap_type, year) %>%
  summarise(
    effort = sum(effort, na.rm = TRUE),
    catch = if (sum(effort, na.rm = TRUE) > 0) sum(catch, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

if (nrow(cell_data) == 0L) {
  stop("No records remain after the spatial filter.")
}

sites <- sort(unique(cell_data$site_name))
years <- sort(unique(cell_data$year))
traps <- CANONICAL_TRAPS

dnames <- list(
  site = sites,
  trap_type = traps,
  year = as.character(years)
)
idx <- cbind(
  match(cell_data$site_name, sites),
  match(cell_data$trap_type, traps),
  match(cell_data$year, years)
)

traptype <- array(0L, dim = c(length(sites), length(traps), length(years)),
                  dimnames = dnames)
traptype[idx] <- as.integer(cell_data$effort > 0)

presence <- array(0L, dim = dim(traptype), dimnames = dimnames(traptype))
presence[idx] <- as.integer(
  cell_data$effort > 0 & !is.na(cell_data$catch) & cell_data$catch > 0
)

cpue_fukui <- matrix(
  NA_real_,
  nrow = length(sites),
  ncol = length(years),
  dimnames = list(site = sites, year = as.character(years))
)
fukui <- cell_data %>%
  filter(trap_type == "Fukui", effort > 0) %>%
  mutate(cpue = catch / effort)
cpue_fukui[cbind(
  match(fukui$site_name, sites),
  match(fukui$year, years)
)] <- fukui$cpue

assert_model_products_aligned(presence, traptype, cpue_fukui)

output_files <- c(
  presence = file.path(outdir, "PresenceArrayBinary.rds"),
  traptype = file.path(outdir, "TrapTypeArraysNew.rds"),
  cpue = file.path(outdir, "cpue_fukui.rds")
)
saveRDS(presence, output_files[["presence"]])
saveRDS(traptype, output_files[["traptype"]])
saveRDS(cpue_fukui, output_files[["cpue"]])

# Re-read the serialized products and validate what downstream code will see.
saved_presence <- readRDS(output_files[["presence"]])
saved_traptype <- readRDS(output_files[["traptype"]])
saved_cpue <- readRDS(output_files[["cpue"]])
assert_model_products_aligned(saved_presence, saved_traptype, saved_cpue)

message(
  "Validated and saved: ", length(sites), " sites x ",
  length(traps), " trap types x ", length(years), " years (",
  min(years), "-", max(years), ")."
)
message("  PresenceArrayBinary.rds: ", sum(presence), " detections")
message("  TrapTypeArraysNew.rds: ", sum(traptype), " sampled cells")
message("  cpue_fukui.rds: ", sum(!is.na(cpue_fukui)),
        " Fukui-sampled site-years")
