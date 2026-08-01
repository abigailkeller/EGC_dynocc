################################################################################
# 02_BayesianTrapData.R
#
# Build the model inputs for the Bayesian dynamic occupancy model from the
# combined monitoring dataset produced by 01_summarize_data_all_traptypes.R.
#
# STRUCTURE (changed)
# -------------------
# Everything is now indexed [site, year, trap], where `trap` is a REPLICATE
# SLOT within a site-year, not a fixed trap-type position. Site-year (i, t)
# deployed n_traps[i, t] trap types; those occupy slots 1..n_traps[i, t] in
# canonical order (Fukui, Shrimp, Minnow). Slots beyond n_traps[i, t] -- and
# every slot of a site-year that was never sampled -- are NA.
#
# So a cell is "real" iff it is not NA, and the number of real cells equals the
# number of site/year/trap-type records that survive the filters (917).
#
# FILTERS (unchanged from the previous version)
#   * year >= 2018
#   * Shrimp, Fukui and Minnow traps only. Case, punctuation and qualified
#     variants ("Shrimp - C", "Minnow - Unmodified", ...) collapse to the three
#     canonical names.
#   * Sites whose representative coordinate falls inside
#     data/SpatialData/study_extent.shp.
#   * NO completeness filter: sites and years are kept however sparse they are.
#
# OUTPUTS (data/model_data/)
#   PresenceArray.rds  int [site, year, trap]  1 = EGC caught in that trap type
#                                              at that site-year, 0 = not
#                                              caught, NA = slot unused.
#   TrapType.rds       named list of three int [site, year, trap] arrays,
#                      $Fukui, $Shrimp, $Minnow. 1 = that slot is that trap
#                      type, 0 = it is one of the others, NA = slot unused.
#                      Across the three, every non-NA cell sums to exactly 1.
#                      Access as tt <- readRDS("TrapType.rds"); tt$Fukui
#   cpue.rds           num [site, year]        pooled catch / pooled effort
#                                              over all trap types deployed.
#
# Every product uses the same site and year labels in the same order; the
# script asserts this before and after writing.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
})

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it with install.packages('here').")
}
here::i_am("code/data_prep/02_BayesianTrapData.R")

FIRST_YEAR      <- 2018L
CANONICAL_TRAPS <- c("Fukui", "Shrimp", "Minnow")  # slot fill order

canonical_trap_type <- function(x) {
  cleaned <- tolower(trimws(as.character(x)))
  cleaned <- gsub("[[:punct:]]+", " ", cleaned)
  cleaned <- gsub("[[:space:]]+", " ", cleaned)
  result <- rep(NA_character_, length(cleaned))
  result[grepl("\\bshrimp\\b", cleaned)] <- "Shrimp"
  result[is.na(result) & grepl("\\bfukui\\b", cleaned)]  <- "Fukui"
  result[is.na(result) & grepl("\\bminnow\\b", cleaned)] <- "Minnow"
  result
}

input_file  <- here::here("data", "data_all_sources_all_traptypes.csv")
extent_file <- here::here("data", "SpatialData", "study_extent.shp")
outdir      <- here::here("data", "model_data")

if (!file.exists(input_file)) {
  stop("Missing ", input_file, ". Run 01_summarize_data_all_traptypes.R first.")
}
if (!file.exists(extent_file)) stop("Missing ", extent_file, ".")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Read and filter
# ---------------------------------------------------------------------------
dat <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

required_columns <- c("site_name", "year", "trap_type", "catch", "effort",
                      "latitude", "longitude")
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0L) {
  stop("Input is missing required columns: ",
       paste(missing_columns, collapse = ", "))
}

dat <- dat %>%
  transmute(
    site_name = trimws(as.character(site_name)),
    year      = suppressWarnings(as.integer(as.character(year))),
    trap_type = canonical_trap_type(trap_type),
    catch     = suppressWarnings(as.numeric(catch)),
    effort    = suppressWarnings(as.numeric(effort)),
    latitude  = suppressWarnings(as.numeric(latitude)),
    longitude = suppressWarnings(as.numeric(longitude))
  ) %>%
  filter(nzchar(site_name), !is.na(year), year >= FIRST_YEAR, !is.na(trap_type))

if (nrow(dat) == 0L) stop("No records remain after the year/trap-type filters.")
if (any(dat$catch < 0, na.rm = TRUE) || any(dat$effort < 0, na.rm = TRUE)) {
  stop("Catch and effort must be non-negative.")
}

# One representative coordinate per site; keep sites inside the study extent.
site_points <- dat %>%
  group_by(site_name) %>%
  summarise(
    latitude  = if (all(is.na(latitude)))  NA_real_ else mean(latitude,  na.rm = TRUE),
    longitude = if (all(is.na(longitude))) NA_real_ else mean(longitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(latitude), is.finite(longitude))

extent  <- st_make_valid(st_read(extent_file, quiet = TRUE))
site_sf <- st_as_sf(site_points, coords = c("longitude", "latitude"),
                    crs = 4326, remove = FALSE)
site_sf <- st_transform(site_sf, st_crs(extent))
sites_in_extent <- site_sf$site_name[lengths(st_intersects(site_sf, extent)) > 0L]

message(length(sites_in_extent), " sites inside the study extent; ",
        length(setdiff(unique(dat$site_name), sites_in_extent)),
        " outside or missing coordinates.")

# ---------------------------------------------------------------------------
# Collapse source/spelling variants to one row per site/year/trap type, then
# assign each surviving row a replicate slot within its site-year.
# Rows with effort but missing catch are treated as zero catch. Rows without
# positive effort were not sampled and get no slot.
# ---------------------------------------------------------------------------
cell_data <- dat %>%
  filter(site_name %in% sites_in_extent) %>%
  group_by(site_name, year, trap_type) %>%
  summarise(
    effort = sum(effort, na.rm = TRUE),
    catch  = sum(catch,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(effort > 0) %>%
  mutate(trap_type = factor(trap_type, levels = CANONICAL_TRAPS)) %>%
  arrange(site_name, year, trap_type) %>%
  group_by(site_name, year) %>%
  mutate(slot = row_number()) %>%
  ungroup()

if (nrow(cell_data) == 0L) stop("No records remain after the spatial filter.")

sites <- sort(unique(cell_data$site_name))
years <- sort(unique(cell_data$year))
ntrap <- max(cell_data$slot)          # maximum number of traps in any site-year

dims   <- c(length(sites), length(years), ntrap)
dnames <- list(site = sites,
               year = as.character(years),
               trap = as.character(seq_len(ntrap)))

idx <- cbind(match(cell_data$site_name, sites),
             match(cell_data$year,      years),
             cell_data$slot)

fill <- function(values, mode = "integer") {
  a <- array(if (mode == "integer") NA_integer_ else NA_real_,
             dim = dims, dimnames = dnames)
  a[idx] <- values
  a
}

# ---------------------------------------------------------------------------
# 1. Presence  [site, year, trap]
# ---------------------------------------------------------------------------
presence <- fill(as.integer(cell_data$catch > 0))

# ---------------------------------------------------------------------------
# 2. One indicator array per trap type, same shape, same NA pattern, bundled
#    into a single named list written to TrapType.rds
# ---------------------------------------------------------------------------
trap_indicators <- lapply(
  setNames(CANONICAL_TRAPS, CANONICAL_TRAPS),
  function(tt) fill(as.integer(as.character(cell_data$trap_type) == tt))
)

# ---------------------------------------------------------------------------
# 3. CPUE  [site, year] -- pooled catch / pooled effort over all trap types
#    deployed in that site-year. NA where the site-year was never sampled.
# ---------------------------------------------------------------------------
pooled <- cell_data %>%
  group_by(site_name, year) %>%
  summarise(cpue = sum(catch) / sum(effort), .groups = "drop")

cpue <- matrix(NA_real_, nrow = length(sites), ncol = length(years),
               dimnames = list(site = sites, year = as.character(years)))
cpue[cbind(match(pooled$site_name, sites),
           match(pooled$year, years))] <- pooled$cpue

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
assert_products_aligned <- function(presence, trap_indicators, cpue) {
  used     <- !is.na(presence)
  n_used   <- apply(used, c(1, 2), sum)   # slots used per site-year
  stopifnot(
    length(dim(presence)) == 3L,
    identical(names(dimnames(presence)), c("site", "year", "trap")),
    identical(names(dimnames(cpue)), c("site", "year")),
    identical(dimnames(presence)$site, dimnames(cpue)$site),
    identical(dimnames(presence)$year, dimnames(cpue)$year),
    !anyDuplicated(dimnames(presence)$site),
    !anyDuplicated(dimnames(presence)$year),
    all(presence[used] %in% c(0L, 1L)),
    # slots are filled from 1 upward with no gaps
    all(apply(used, c(1, 2), function(z) all(diff(as.integer(z)) <= 0))),
    # trap-type indicators share the NA pattern and are one-hot
    all(vapply(trap_indicators,
               function(a) identical(!is.na(a), used), logical(1))),
    all(Reduce(`+`, trap_indicators)[used] == 1L),
    # CPUE present exactly where the site-year was sampled
    all(unname(!is.na(cpue)) == unname(n_used > 0L)),
    all(is.na(cpue) | (is.finite(cpue) & cpue >= 0))
  )
  invisible(TRUE)
}

assert_products_aligned(presence, trap_indicators, cpue)

# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------
saveRDS(presence,        file.path(outdir, "PresenceArray.rds"))
saveRDS(trap_indicators, file.path(outdir, "TrapType.rds"))
saveRDS(cpue,            file.path(outdir, "cpue.rds"))

n_used <- apply(!is.na(presence), c(1, 2), sum)
message("Saved [site, year, trap] = ", dims[1], " x ", dims[2], " x ", dims[3],
        " (", min(years), "-", max(years), ").")
message("  ", sum(!is.na(presence)), " used trap slots (site/year/trap-type records)")
message("  ", sum(presence, na.rm = TRUE), " with EGC detected")
message("  ", sum(n_used > 0L), " sampled site-years")
message("  slots per site-year: ",
        paste(names(table(n_used[n_used > 0L])), table(n_used[n_used > 0L]),
              sep = "=", collapse = ", "))
