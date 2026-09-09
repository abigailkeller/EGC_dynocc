# Separate WSG observations in the same [site, year, replicate] format as
# PresenceArrayBinary.rds. Here a replicate is a pooled sampling event.
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
  # Same spatial rule as the individual-trap data: mean coordinates per site.
  coords <- events %>% group_by(site) %>% summarise(
    latitude = mean(latitude, na.rm = TRUE),
    longitude = mean(longitude, na.rm = TRUE), .groups = "drop"
  ) %>% filter(is.finite(latitude), is.finite(longitude))
  extent <- st_make_valid(st_read(
    here::here("data", "SpatialData", "study_extent.shp"), quiet = TRUE))
  points <- st_transform(st_as_sf(coords, coords = c("longitude", "latitude"),
                                  crs = 4326), st_crs(extent))
  keep <- points$site[lengths(st_intersects(points, extent)) > 0L]
  events <- events %>% filter(site %in% keep) %>%
    arrange(site, year, date, raw_row) %>% group_by(site, year) %>%
    mutate(replicate = row_number(), presence = as.integer(catch > 0)) %>% ungroup()
  if (!nrow(events)) stop("No WSG sampling events remain after filtering.")
  sites <- sort(unique(events$site))
  dn <- list(site = sites, year = years,
             replicate = as.character(seq_len(max(events$replicate))))
  idx <- cbind(match(events$site, sites), match(events$year, years), events$replicate)
  fill <- function(values) {
    a <- array(NA_integer_, lengths(dn), dimnames = dn)
    a[idx] <- as.integer(values)
    a
  }
  presence <- fill(events$presence)
  effort <- fill(events$effort)
  # The raw file supplies total effort only. Even six traps does not establish
  # the 3 Minnow + 3 Fukui mix. Leave counts unknown pending protocol validation.
  counts <- setNames(lapply(c("Fukui", "Shrimp", "Minnow"),
                           function(x) fill(rep(NA_integer_, nrow(events)))),
                     c("Fukui", "Shrimp", "Minnow"))
  stopifnot(identical(dimnames(presence)$year, dimnames(existing)$year),
            sum(!is.na(presence)) == nrow(events),
            sum(presence, na.rm = TRUE) == sum(events$catch > 0),
            identical(is.na(presence), is.na(effort)),
            all(presence[!is.na(presence)] %in% 0:1))
  saveRDS(presence, file.path(outdir, "WSGPresenceArrayBinary.rds"))
  saveRDS(effort, file.path(outdir, "WSGTrapEffortArray.rds"))
  saveRDS(counts, file.path(outdir, "WSGTrapCountArrays.rds"))
  saveRDS(events, file.path(outdir, "WSGSamplingEvents.rds"))
  write.csv(events, file.path(outdir, "WSGSamplingEvents.csv"), row.names = FALSE)
  message("Saved separate WSG [site, year, replicate] = ",
          paste(dim(presence), collapse = " x "), "; ", nrow(events),
          " events, ", sum(events$presence), " detections; ",
          sum(events$effort != 6), " events with effort other than six.")
})
