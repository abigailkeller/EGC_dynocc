################################################################################
# Map the sites in the replicate-level model matrices to spatial zones.
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
here::i_am("code/data_prep/03_map_model_sites_to_zones.R")

# Recreate the filtered replicate table using the exact model-data preparation
# logic, while redirecting its serialized outputs to a temporary directory.
prep_file <- here::here("code", "data_prep", "02_BayesianTrapData.R")
prep_code <- readLines(prep_file, warn = FALSE)
prep_code <- sub(
  'outdir <- here::here\\("data", "model_data"\\)',
  'outdir <- file.path(tempdir(), "model_data")',
  prep_code
)
# Only the individual-trap records are needed from this script. The saved WSG
# artifacts are read separately below; do not rebuild them as a side effect.
prep_code <- prep_code[prep_code != paste0(
  'source(here::here("code", "data_prep", ',
  '"02b_WSGEventData.R"))'
)]
prep_env <- new.env(parent = globalenv())
eval(parse(text = prep_code), envir = prep_env)

zones <- st_read(here::here("data", "SpatialData", "zones.shp"), quiet = TRUE) %>%
  st_make_valid()

map_sites_to_zones <- function(site_data, model_sites, dataset_label) {
  stopifnot(!anyDuplicated(site_data$site))
  site_data <- site_data %>%
    filter(site %in% model_sites) %>%
    arrange(match(site, model_sites))
  if (!identical(site_data$site, model_sites)) {
    missing_sites <- setdiff(model_sites, site_data$site)
    stop(dataset_label, " sites missing usable coordinates: ",
         paste(missing_sites, collapse = ", "))
  }

  output_columns <- setdiff(names(site_data), c("latitude", "longitude"))
  site_points <- site_data %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(st_crs(zones))
  zone_hits <- st_intersects(site_points, zones)

  if (any(lengths(zone_hits) > 1L)) {
    multi <- which(lengths(zone_hits) > 1L)
    for (i in multi) {
      candidates <- zone_hits[[i]]
      boundary_distance <- as.numeric(st_distance(
        site_points[i, ], st_boundary(zones[candidates, ])
      ))
      selected <- candidates[which.max(boundary_distance)]
      message(dataset_label, " overlapping-zone resolution: ",
              site_points$site[i], " -> ", zones$name[selected],
              " (farthest from candidate boundaries)")
      zone_hits[[i]] <- selected
    }
  }

  # A few shoreline coordinates fall just outside the hand-drawn zone polygons.
  # Assign these to the nearest zone rather than leaving model sites unmapped.
  outside <- which(lengths(zone_hits) == 0L)
  if (length(outside)) {
    nearest <- st_nearest_feature(site_points[outside, ], zones)
    distances_m <- as.numeric(st_distance(
      site_points[outside, ], zones[nearest, ], by_element = TRUE
    ))
    near_boundary <- distances_m <= 1000
    zone_hits[outside[near_boundary]] <- lapply(
      nearest[near_boundary], function(i) i
    )
    if (any(near_boundary)) {
      message(dataset_label, " nearest-zone shoreline fallback: ", paste0(
        site_points$site[outside[near_boundary]], " -> ",
        zones$name[nearest[near_boundary]], " (",
        round(distances_m[near_boundary]), " m)", collapse = "; "
      ))
    }
    if (any(!near_boundary)) {
      message(dataset_label, " outside all zones: ", paste0(
        site_points$site[outside[!near_boundary]], " (nearest polygon ",
        round(distances_m[!near_boundary] / 1000, 1), " km away)",
        collapse = "; "
      ))
    }
  }

  zone_index <- vapply(zone_hits, function(x) {
    if (length(x) == 1L) x else NA_integer_
  }, integer(1))
  site_zone_map <- cbind(
    st_drop_geometry(site_points)[output_columns],
    zone_id = zones$zone_id[zone_index],
    zone = zones$name[zone_index]
  )
  stopifnot(
    identical(site_zone_map$site, model_sites),
    !anyDuplicated(site_zone_map$site),
    all(!is.na(site_zone_map$zone_id) == !is.na(site_zone_map$zone))
  )
  site_zone_map
}

model_sites <- dimnames(readRDS(here::here(
  "data", "model_data", "PresenceArrayBinary.rds"
)))$site
site_data <- prep_env$reps %>%
  group_by(site) %>%
  summarise(
    latitude = mean(latitude[is.finite(latitude)]),
    longitude = mean(longitude[is.finite(longitude)]),
    .groups = "drop"
  )
site_zone_map <- map_sites_to_zones(site_data, model_sites, "Model")

wsg_presence <- readRDS(here::here(
  "data", "model_data", "PresenceArrayBinary_WSG.rds"
))
wsg_sites <- dimnames(wsg_presence)$site
wsg_events <- readRDS(here::here(
  "data", "model_data", "SamplingEvents_WSG.rds"
))
stopifnot(all(
  wsg_events %>% group_by(site) %>% summarise(n = n_distinct(site_name)) %>%
    pull(n) == 1L
))
wsg_site_data <- wsg_events %>%
  group_by(site) %>%
  summarise(
    site_name = first(site_name),
    latitude = mean(latitude[is.finite(latitude)]),
    longitude = mean(longitude[is.finite(longitude)]),
    .groups = "drop"
  )
wsg_site_zone_map <- map_sites_to_zones(wsg_site_data, wsg_sites, "WSG")

outputs <- list(
  site_zone_map = site_zone_map,
  site_zone_map_WSG = wsg_site_zone_map
)
for (nm in names(outputs)) {
  output_file <- here::here("data", "model_data", paste0(nm, ".csv"))
  write.csv(outputs[[nm]], output_file, row.names = FALSE)
  message("Saved ", nrow(outputs[[nm]]), " site-to-zone mappings to ",
          output_file)
}
