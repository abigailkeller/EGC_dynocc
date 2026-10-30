# Separate WSG observations in the same [site, year, replicate] format as
# PresenceArrayBinary.rds. Here a replicate is a calendar month (June-September).
# Run after 02_BayesianTrapData.R (which also sources this script).
local({
  suppressPackageStartupMessages({
    library(dplyr)
    library(sf)
  })
  outdir <- here::here("data", "model_data")
  existing <- readRDS(file.path(outdir, "PresenceArrayBinary.rds"))
  years <- dimnames(existing)$year
  raw <- read.csv(here::here("data", "raw", "WSG_Traps.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
  events <- data.frame(
    source = "WSG", site_name = trimws(raw[["Site Name"]]),
    site = paste0("WSG::", raw$SiteID),
    year = as.integer(raw$year),
    date = as.Date(raw$EndTime, "%m/%d/%Y"),
    month = match(raw$Month, month.name),
    catch = as.numeric(raw$total.cama), effort = as.numeric(raw$trap.sets),
    latitude = as.numeric(raw$Latitude), longitude = as.numeric(raw$Longitude),
    raw_row = seq_len(nrow(raw))
  ) %>% filter(year %in% years, month %in% 6:9,
               is.finite(catch), catch >= 0,
               is.finite(effort), effort > 0)
  stopifnot(!anyNA(events$date), all(events$effort == floor(events$effort)))
  # Same expanded spatial rule as the individual-trap data: retain every site
  # with usable mean coordinates and record study-extent membership for audit.
  coords <- events %>% group_by(site) %>% summarise(
    source = first(source), site_name = first(site_name),
    latitude = mean(latitude, na.rm = TRUE),
    longitude = mean(longitude, na.rm = TRUE), .groups = "drop"
  ) %>% filter(is.finite(latitude), is.finite(longitude))
  extent <- st_make_valid(st_read(
    here::here("data", "SpatialData", "study_extent.shp"), quiet = TRUE))
  points <- st_transform(st_as_sf(coords, coords = c("longitude", "latitude"),
                                  crs = 4326), st_crs(extent))
  coords$in_study_extent <- lengths(st_intersects(points, extent)) > 0L
  wsg_site_coordinates <- coords %>%
    transmute(source, site, site_name, x = longitude, y = latitude,
              in_study_extent) %>%
    arrange(site_name)
  write.csv(wsg_site_coordinates,
            file.path(outdir, "site_coordinates_WSG.csv"), row.names = FALSE)
  individual_site_coordinates <- read.csv(
    file.path(outdir, "site_coordinates_individual_traps.csv"),
    stringsAsFactors = FALSE
  )
  expanded_site_coordinates <- bind_rows(
    individual_site_coordinates, wsg_site_coordinates
  ) %>% arrange(source, site_name)
  write.csv(expanded_site_coordinates,
            file.path(outdir, "site_coordinates_expanded.csv"),
            row.names = FALSE)
  write.csv(filter(expanded_site_coordinates, !in_study_extent),
            file.path(outdir, "site_coordinates_outside_extent.csv"),
            row.names = FALSE)
  months <- 6:9
  build_wsg_product <- function(data, retained_sites, suffix = "", label) {
    retained_events <- data %>%
      filter(site %in% retained_sites) %>%
      arrange(site, year, date, raw_row) %>%
      mutate(replicate = match(month, months), presence = as.integer(catch > 0))
    if (!nrow(retained_events)) stop("No ", label, " WSG events remain after filtering.")

    # Pool all records in a site/year/month: any positive catch is a detection.
    occasions <- retained_events %>%
      group_by(site, year, month, replicate) %>%
      summarise(presence = as.integer(any(presence == 1L)), .groups = "drop")
    sites <- sort(unique(retained_events$site))
    dn <- list(site = sites, year = years, replicate = as.character(months))
    idx <- cbind(match(occasions$site, sites), match(occasions$year, years),
                 occasions$replicate)
    presence <- array(NA_integer_, lengths(dn), dimnames = dn)
    presence[idx] <- occasions$presence

    # WSG uses a pooled monthly observation under the usual six-trap protocol
    # (three Minnow + three Fukui). It therefore has no individual-trap CPUE or
    # trap-type arrays comparable to the other sources.
    stopifnot(
      identical(dimnames(presence)$year, dimnames(existing)$year),
      sum(!is.na(presence)) == nrow(occasions),
      sum(presence, na.rm = TRUE) == sum(occasions$presence),
      all(presence[!is.na(presence)] %in% 0:1)
    )
    saveRDS(presence,
            file.path(outdir, paste0("PresenceArrayBinary_WSG", suffix, ".rds")))
    saveRDS(retained_events,
            file.path(outdir, paste0("SamplingEvents_WSG", suffix, ".rds")))
    write.csv(retained_events,
              file.path(outdir, paste0("SamplingEvents_WSG", suffix, ".csv")),
              row.names = FALSE)
    message("Saved ", label, " WSG [site, year, replicate] = ",
            paste(dim(presence), collapse = " x "), "; ", nrow(occasions),
            " monthly occasions, ", sum(occasions$presence), " detections; ",
            sum(retained_events$effort != 6), " events with effort other than six.")
    invisible(list(presence = presence, events = retained_events,
                   occasions = occasions))
  }

  original_wsg <- build_wsg_product(
    events, coords$site[coords$in_study_extent], label = "original"
  )
  expanded_wsg <- build_wsg_product(
    events, coords$site, suffix = "_expanded", label = "expanded"
  )
  original_wsg_index <- match(dimnames(original_wsg$presence)$site,
                              dimnames(expanded_wsg$presence)$site)
  stopifnot(
    !anyNA(original_wsg_index),
    identical(
      original_wsg$presence,
      expanded_wsg$presence[original_wsg_index, , , drop = FALSE]
    )
  )
})
