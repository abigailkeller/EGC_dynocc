################################################################################
# First-pass inventory of apparent occupancy transitions under alternative
# temporal and spatial scopes.
#
# "Apparent colonization" means no EGC detected at a site in year t-1 and at
# least one EGC detected at that site in year t. This is an exploratory proxy,
# not an estimate of true colonization under imperfect detection.
################################################################################

project_library <- file.path(getwd(), ".r-library")
if (dir.exists(project_library)) {
  .libPaths(c(normalizePath(project_library), .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
})

combined <- read.csv(
  "data/data_all_sources_all_traptypes.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  mutate(
    year = as.integer(year),
    canonical_trap = grepl("fukui|shrimp|minnow", trap_type, ignore.case = TRUE),
    replicate_capable_source = grepl("WDFW|Makah|DNWR|DFO", data_sources)
  )

# The pre-combined table intentionally used the older WDFW file for 2018-2022
# and the newer file for 2023-2024. The newer file contains some additional
# sites in the overlapping years (for example, Maynard Lagoon in 2022). Add the
# raw WDFW records here so this exploration does not miss those site-years.
clean_wdfw_site <- function(x) {
  x <- trimws(as.character(x))
  replacements <- c(
    "Duckland & 105" = "Duckland",
    "Ocean Shores - Airport" = "Ocean Shores",
    "Ocean Shores Airport" = "Ocean Shores",
    "West Samish" = "Samish River",
    "North" = "Crandall Spit",
    "Tokeland Hotel" = "Tokeland"
  )
  hit <- match(x, names(replacements))
  x[!is.na(hit)] <- unname(replacements[hit[!is.na(hit)]])
  x
}

wdfw_old <- read.csv(
  "data/raw/WDFW/WDFW.EGC_2018-2022_Effort_Final.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
wdfw_old_date <- as.Date(wdfw_old[["Date_Deployed"]], "%m/%d/%Y")
wdfw_old_extra <- tibble(
  site_name = clean_wdfw_site(wdfw_old[["Site_Name"]]),
  year = as.integer(format(wdfw_old_date, "%Y")),
  month = as.integer(format(wdfw_old_date, "%m")),
  trap_type = wdfw_old[["Trap_Type"]],
  data_sources = "WDFW raw",
  catch = suppressWarnings(as.numeric(wdfw_old[["CAMA_Total"]])),
  effort = 1,
  cpue = catch,
  latitude = suppressWarnings(as.numeric(wdfw_old[["Latitude"]])),
  longitude = suppressWarnings(as.numeric(wdfw_old[["Longitude"]]))
)

wdfw_new <- read.csv(
  "data/raw/WDFW/WDFW_EGC_Data_Collection_PUBLIC_09092025_pg1.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
wdfw_new_date <- as.Date(wdfw_new[["Set Date and Time"]], "%m/%d/%Y")
wdfw_new_extra <- tibble(
  site_name = clean_wdfw_site(wdfw_new[["Site Name"]]),
  year = as.integer(format(wdfw_new_date, "%Y")),
  month = as.integer(format(wdfw_new_date, "%m")),
  trap_type = wdfw_new[["Trap Type"]],
  data_sources = "WDFW raw",
  catch = suppressWarnings(as.numeric(wdfw_new[["Total Catch"]])),
  effort = 1,
  cpue = catch,
  latitude = suppressWarnings(as.numeric(wdfw_new[["y"]])),
  longitude = suppressWarnings(as.numeric(wdfw_new[["x"]]))
)

wdfw_extra <- bind_rows(wdfw_old_extra, wdfw_new_extra) %>%
  filter(
    month %in% 6:9,
    is.finite(catch), catch >= 0,
    nzchar(site_name),
    grepl("fukui|shrimp|minnow", trap_type, ignore.case = TRUE)
  ) %>%
  select(-month) %>%
  mutate(
    canonical_trap = TRUE,
    replicate_capable_source = TRUE
  )

combined <- bind_rows(combined, wdfw_extra)

# Assign the combined-data sites to the current study polygon and connectivity
# zones using the single site coordinate stored in the combined dataset.
site_points <- combined %>%
  group_by(site_name) %>%
  summarise(
    latitude = first(latitude),
    longitude = first(longitude),
    .groups = "drop"
  ) %>%
  filter(is.finite(latitude), is.finite(longitude))

extent <- st_make_valid(st_read("data/SpatialData/study_extent.shp", quiet = TRUE))
zones <- st_make_valid(st_read("data/SpatialData/zones.shp", quiet = TRUE))
site_sf <- st_as_sf(
  site_points,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(st_crs(extent))

site_points$in_study_extent <- lengths(st_intersects(site_sf, extent)) > 0L
zone_hits <- st_intersects(st_transform(site_sf, st_crs(zones)), zones)
site_points$zone_id <- vapply(
  zone_hits,
  function(idx) if (length(idx)) as.integer(zones$zone_id[idx[1]]) else NA_integer_,
  integer(1)
)

combined <- combined %>%
  left_join(
    site_points %>% select(site_name, in_study_extent, zone_id),
    by = "site_name"
  )

transition_counts <- function(status) {
  status <- status %>% arrange(site_name, year)
  pairs <- status %>%
    group_by(site_name) %>%
    mutate(
      previous_year = lag(year),
      previous_detected = lag(detected),
      year_gap = year - previous_year
    ) %>%
    ungroup() %>%
    filter(!is.na(previous_detected))

  adjacent <- pairs %>% filter(year_gap == 1L)
  successive <- pairs

  count_transition <- function(x, previous, current) {
    sum(x$previous_detected == previous & x$detected == current)
  }

  tibble(
    sites = n_distinct(status$site_name),
    site_years = nrow(status),
    sites_with_2plus_years = sum(table(status$site_name) >= 2L),
    adjacent_pairs = nrow(adjacent),
    adjacent_00 = count_transition(adjacent, FALSE, FALSE),
    adjacent_01_colonizations = count_transition(adjacent, FALSE, TRUE),
    adjacent_10_losses = count_transition(adjacent, TRUE, FALSE),
    adjacent_11_persistence = count_transition(adjacent, TRUE, TRUE),
    successive_visit_pairs = nrow(successive),
    successive_01_detections = count_transition(successive, FALSE, TRUE)
  )
}

summarise_scope <- function(data, scope_name) {
  status <- data %>%
    group_by(site_name, year) %>%
    summarise(
      detected = any(catch > 0, na.rm = TRUE),
      effort = sum(effort, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(effort > 0)

  bind_cols(tibble(scope = scope_name), transition_counts(status))
}

# Exact current input, retaining the source-prefixed site identities used by
# the fitted models.
presence <- readRDS("data/model_data/PresenceArrayBinary.rds")[, 2:7, , drop = FALSE]
site_zone <- read.csv("data/model_data/site_zone_map.csv", stringsAsFactors = FALSE)
keep <- !apply(presence, 1, function(x) all(is.na(x)))
presence <- presence[keep, , , drop = FALSE]
site_zone <- site_zone[keep, , drop = FALSE]
keep <- is.na(site_zone$zone_id) | site_zone$zone_id != 2L
presence <- presence[keep, , , drop = FALSE]

current_status_matrix <- apply(
  presence,
  c(1, 2),
  function(x) if (all(is.na(x))) NA else any(x == 1L, na.rm = TRUE)
)
current_idx <- which(!is.na(current_status_matrix), arr.ind = TRUE)
current_status <- tibble(
  site_name = rownames(current_status_matrix)[current_idx[, 1]],
  year = as.integer(colnames(current_status_matrix)[current_idx[, 2]]),
  detected = current_status_matrix[current_idx]
)

results <- bind_rows(
  bind_cols(tibble(scope = "Exact current model input: 2019-2024, zone 2 excluded"),
            transition_counts(current_status)),
  summarise_scope(
    combined %>% filter(canonical_trap, year >= 2018L,
                        in_study_extent, is.na(zone_id) | zone_id != 2L),
    "Combined data: 2018-2024, current extent, zone 2 excluded"
  ),
  summarise_scope(
    combined %>% filter(canonical_trap, in_study_extent,
                        is.na(zone_id) | zone_id != 2L),
    "Combined data: all years, current extent, zone 2 excluded"
  ),
  summarise_scope(
    combined %>% filter(canonical_trap, year >= 2014L, !is.na(zone_id)),
    "Combined data: 2014-2024, all connectivity zones"
  ),
  summarise_scope(
    combined %>% filter(canonical_trap, !is.na(zone_id)),
    "Combined data: all years, all connectivity zones"
  ),
  summarise_scope(
    combined %>% filter(canonical_trap, !in_study_extent),
    "Combined data: all years, outside current extent"
  ),
  summarise_scope(
    combined %>% filter(canonical_trap),
    "Combined data: all years and locations, three focal trap types"
  ),
  summarise_scope(
    combined,
    "Combined data: all years and locations, all trap types"
  ),
  summarise_scope(
    combined %>% filter(canonical_trap, replicate_capable_source),
    "Combined data: all years/locations, replicate-capable sources"
  )
)

write.csv(
  results,
  "data/processed/colonization_scope_summary.csv",
  row.names = FALSE
)

print(results, width = Inf)
