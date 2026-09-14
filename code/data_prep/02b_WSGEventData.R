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
  months <- 6:9
  events <- events %>% filter(site %in% keep) %>%
    arrange(site, year, date, raw_row) %>%
    mutate(replicate = match(month, months), presence = as.integer(catch > 0))
  if (!nrow(events)) stop("No WSG sampling events remain after filtering.")
  # Pool all records in a site/year/month: any positive catch is a detection.
  occasions <- events %>% group_by(site, year, month, replicate) %>%
    summarise(presence = as.integer(any(presence == 1L)), .groups = "drop")
  sites <- sort(unique(events$site))
  dn <- list(site = sites, year = years,
             replicate = as.character(months))
  idx <- cbind(match(occasions$site, sites), match(occasions$year, years),
               occasions$replicate)
  fill <- function(values) {
    a <- array(NA_integer_, lengths(dn), dimnames = dn)
    a[idx] <- as.integer(values)
    a
  }
  presence <- fill(occasions$presence)
  # WSG uses a pooled monthly observation under the usual six-trap protocol
  # (three Minnow + three Fukui). Retain recorded effort in the audit table;
  # no separate effort or trap-type arrays are needed for this approach.
  stopifnot(identical(dimnames(presence)$year, dimnames(existing)$year),
            sum(!is.na(presence)) == nrow(occasions),
            sum(presence, na.rm = TRUE) == sum(occasions$presence),
            all(presence[!is.na(presence)] %in% 0:1))
  saveRDS(presence, file.path(outdir, "PresenceArrayBinary_WSG.rds"))
  saveRDS(events, file.path(outdir, "SamplingEvents_WSG.rds"))
  write.csv(events, file.path(outdir, "SamplingEvents_WSG.csv"), row.names = FALSE)
  message("Saved separate WSG [site, year, replicate] = ",
          paste(dim(presence), collapse = " x "), "; ", nrow(occasions),
          " monthly occasions, ", sum(occasions$presence), " detections; ",
          sum(events$effort != 6), " events with effort other than six.")
})
