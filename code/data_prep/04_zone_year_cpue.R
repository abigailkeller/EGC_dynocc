################################################################################
# 04_zone_year_cpue.R
#
# Calculate one effort-weighted Fukui CPUE for every zone-year combination from
# the complete summer summary produced by 01_summarize_data_all_traptypes.R.
# This includes Fukui records from every source, including sources that provide
# only season-pooled catch and effort rather than trap-specific detection
# histories. Those records are suitable for a zone-level abundance covariate
# even though they cannot be used as occupancy-model replicates.
#
# Output:
#   data/model_data/cpue_zone_year.rds
#       numeric [zone, year] matrix; each observed cell is
#       Fukui sum(catch) / sum(effort). Unobserved zone-year cells are NA.
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
here::i_am("code/data_prep/04_zone_year_cpue.R")

input_file <- here::here("data", "data_all_sources_all_traptypes.csv")
output_file <- here::here("data", "model_data", "cpue_zone_year.rds")

all_data <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(
    year = as.integer(year),
    catch = as.numeric(catch),
    effort = as.numeric(effort),
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  ) %>%
  filter(
    trimws(trap_type) == "Fukui",
    is.finite(year), is.finite(catch), catch >= 0,
    is.finite(effort), effort > 0,
    is.finite(latitude), is.finite(longitude)
  )

# Use one stable point per named site. Coordinates in the combined summary can
# vary slightly among sources/years, so the mean prevents a site from changing
# zones because of harmless sampling-position jitter.
site_points <- all_data %>%
  group_by(site_name) %>%
  summarise(
    latitude = mean(latitude),
    longitude = mean(longitude),
    .groups = "drop"
  ) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

zones <- st_read(here::here("data", "SpatialData", "zones.shp"), quiet = TRUE) %>%
  st_make_valid() %>%
  arrange(zone_id)
site_points <- st_transform(site_points, st_crs(zones))
zone_hits <- st_intersects(site_points, zones)

# Resolve overlapping polygons by choosing the polygon whose boundary is
# farthest away, matching the site-to-zone mapping used by the model data.
for (i in which(lengths(zone_hits) > 1L)) {
  candidates <- zone_hits[[i]]
  boundary_distance <- as.numeric(st_distance(
    site_points[i, ], st_boundary(zones[candidates, ])
  ))
  zone_hits[[i]] <- candidates[which.max(boundary_distance)]
}

# Shoreline coordinates can fall just outside hand-drawn polygons. Apply the
# same conservative 1-km nearest-zone fallback as 03_map_model_sites_to_zones.R.
outside <- which(lengths(zone_hits) == 0L)
if (length(outside)) {
  nearest <- st_nearest_feature(site_points[outside, ], zones)
  distances_m <- as.numeric(st_distance(
    site_points[outside, ], zones[nearest, ], by_element = TRUE
  ))
  near_boundary <- distances_m <= 1000
  zone_hits[outside[near_boundary]] <- lapply(nearest[near_boundary], function(i) i)
}

zone_index <- vapply(
  zone_hits,
  function(x) if (length(x) == 1L) x else NA_integer_,
  integer(1)
)
site_zone <- data.frame(
  site_name = site_points$site_name,
  zone_id = zones$zone_id[zone_index],
  zone = zones$name[zone_index]
)

mapped_data <- all_data %>%
  left_join(site_zone, by = "site_name") %>%
  filter(!is.na(zone_id))

zone_year <- mapped_data %>%
  group_by(zone_id, zone, year) %>%
  summarise(
    catch = sum(catch),
    effort = sum(effort),
    cpue = catch / effort,
    .groups = "drop"
  )

all_years <- seq.int(min(all_data$year), max(all_data$year))
cpue_zone_year <- matrix(
  NA_real_,
  nrow = nrow(zones),
  ncol = length(all_years),
  dimnames = list(zone = zones$name, year = as.character(all_years))
)
cpue_zone_year[cbind(
  match(zone_year$zone_id, zones$zone_id),
  match(zone_year$year, all_years)
)] <- zone_year$cpue

stopifnot(
  !anyDuplicated(zones$zone_id),
  all(zone_year$effort > 0),
  all(is.finite(zone_year$cpue)),
  all(cpue_zone_year[!is.na(cpue_zone_year)] >= 0)
)

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(cpue_zone_year, output_file)

message("Saved zone-year CPUE [zone, year] = ",
        paste(dim(cpue_zone_year), collapse = " x "), " to ", output_file)
message("  Included ", nrow(mapped_data), " of ", nrow(all_data),
        " valid Fukui summary rows across all sources.")
message("  Excluded ", sum(is.na(site_zone$zone_id)),
        " sites outside the zone polygons and >1 km from a boundary.")
