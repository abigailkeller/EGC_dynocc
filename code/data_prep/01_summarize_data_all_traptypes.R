################################################################################
# summarize_data_all_traptypes.R
#
# Adapted from FromAbby/EGC_spatialmodel/summarize_data.R
#
# WHAT CHANGED FROM THE ORIGINAL
# ------------------------------
# The original script kept ONLY Fukui traps in every data source (and, for
# Makah, converted shrimp CPUE to a "Fukui-equivalent" using a computed ratio).
# This version KEEPS ALL TRAP TYPES. Every source now carries a `trap_type`
# column and is summarized to catch / effort / CPUE SEPARATELY for each trap
# type within each site-year.
#
# WHAT DID NOT CHANGE (all other filtering is identical to the original)
# ---------------------------------------------------------------------
#   * Season: every source is still restricted to summer (June-September).
#   * DFO: still restricted to usable traps (TrapUsability %in% c(0, 15));
#          catch is still defined as Species == "XMB" (European green crab).
#   * WDFW: still excludes the "Drayton Harbor" waterbody (that site is handled
#           by the separate NW Straits / Drayton Harbor source); still drops
#           rows with NA catch totals.
#   * DNWR: still drops catch records with NA crab number (via the year-by-year
#           joins carried over from the original).
#   * new WDFW file still restricted to 2023-2024.
#   * Output is still summarized CPUE per site-year (now also per trap type),
#     not raw trap-level records.
#
# TRAP-TYPE LABELS (things to be aware of)
# ----------------------------------------
#   * DFO has no usable trap-type text field (its TrapType column is blank), so
#     trap type is derived from GearCode. GearCode 51 = Fukui. Other gear codes
#     present in the data (50, 76, 82, 90, 99, 29, A1) ARE retained but their
#     gear type is unknown here; they are labelled "GearCode_<code>". Map them
#     from DFO's gear-code lookup table if you need real names.
#   * WSG_Traps.csv has no trap-type column at all (it is pre-aggregated). Those
#     data were collected with Fukui traps and are labelled "Fukui".
#   * Drayton Harbor now sums the Minnow (m_index) and Shrimp (s_index) indices
#     in addition to Fukui (f_index). This assumes m_index/s_index share the
#     same [time, trap, year] structure as f_index (they are siblings from the
#     same preprocessing step).
#   * WDFW/Padilla/DNWR trap-type strings are only whitespace-trimmed; distinct
#     variants (e.g. "Shrimp - C", "Shrimp - R") are preserved as-is.
#
# PATHS: resolved from the repository root with the `here` package. Inputs live
# in data/raw; outputs are written to data/processed and data/.
#
# NOTE: this script was written without an R runtime available for testing.
# Column names were validated against the raw files, but please run once and
# sanity-check the outputs.
################################################################################

library(tidyverse)

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it with install.packages('here').")
}
here::i_am("code/data_prep/01_summarize_data_all_traptypes.R")

raw <- here::here("data", "raw")
procd <- here::here("data", "processed")
if (!dir.exists(raw)) {
  stop("Could not find the raw-data directory: ", raw)
}
dir.create(procd, showWarnings = FALSE, recursive = TRUE)

# standard output column order used by every source
# year, site_name, trap_type, catch, effort, cpue, latitude, longitude

#########
# Makah #
#########

data_Makah <- read.csv(file.path(raw, "Makah_catchdata.csv"))

# subset to June - September
data_Makah <- data_Makah[data_Makah$Month %in% c("June", "July",
                                                 "August", "September"), ]

# summarise per trap type (ALL trap types kept; the original Fukui/shrimp
# selection and shrimp -> Fukui conversion has been removed)
summary_Makah <- data_Makah %>%
  group_by(Year, Trap.Type, Location) %>%
  summarize(effort = n(),
            catch = sum(Count.of.EGC),
            longitude = mean(Longitude),
            latitude = mean(Latitude),
            .groups = "drop") %>%
  mutate(cpue = catch / effort) %>%
  transmute(year = Year, site_name = Location, trap_type = Trap.Type,
            catch, effort, cpue, latitude, longitude)

saveRDS(summary_Makah, file.path(procd, "Makah_summer.rds"))

########
# WDFW #
########

## 2018 - 2022 ----------------------------------------------------------------
data_WDFW <- read.csv(file.path(raw, "WDFW",
                                "WDFW.EGC_2018-2022_Effort_Final.csv"))

# trim whitespace so "Fukui " == "Fukui" etc. (ALL trap types kept)
data_WDFW$Trap_Type <- trimws(data_WDFW$Trap_Type)

# remove Drayton Harbor (handled by the NW Straits / Drayton Harbor source)
data_WDFW <- data_WDFW[data_WDFW$Waterbody != "Drayton Harbor", ]

# date columns
data_WDFW$Date_Deployed <- as.Date(data_WDFW$Date_Deployed, "%m/%d/%Y")
data_WDFW$year  <- format(data_WDFW$Date_Deployed, "%Y")
data_WDFW$month <- format(data_WDFW$Date_Deployed, "%m")

# subset to June - September
data_WDFW <- data_WDFW[data_WDFW$month %in% c("06", "07", "08", "09"), ]

# numeric catch; drop NA
data_WDFW$CAMA_Total <- as.numeric(data_WDFW$CAMA_Total)
data_WDFW <- data_WDFW[!is.na(data_WDFW$CAMA_Total), ]

# clean site names (unchanged from original)
data_WDFW[data_WDFW$Site_Name == " Duckland & 105", "Site_Name"] <- "Duckland"
data_WDFW[data_WDFW$Site_Name == " South", "Site_Name"] <- "South"
data_WDFW[data_WDFW$Site_Name == "Henry Island Preserve ", "Site_Name"] <- "Henry Island Preserve"
data_WDFW[data_WDFW$Site_Name == "Ocean Shores - Airport", "Site_Name"] <- "Ocean Shores Airport"
data_WDFW[data_WDFW$Site_Name == "Post Point Lagoon ", "Site_Name"] <- "Post Point Lagoon"
data_WDFW[data_WDFW$Site_Name == "West Samish ", "Site_Name"] <- "West Samish"
data_WDFW[data_WDFW$Site_Name == "North ", "Site_Name"] <- "Crandall Spit"
data_WDFW[data_WDFW$Site_Name == "Ocean Shores Airport", "Site_Name"] <- "Ocean Shores"
data_WDFW[data_WDFW$Site_Name == "Tokeland Hotel", "Site_Name"] <- "Tokeland"
data_WDFW[data_WDFW$Site_Name == "West Samish", "Site_Name"] <- "Samish River"

## 2023 - 2024 ----------------------------------------------------------------
new_WDFW <- read.csv(file.path(raw, "WDFW",
                               "WDFW_EGC_Data_Collection_PUBLIC_09092025_pg1.csv"))

new_WDFW$Date_Deployed <- as.Date(new_WDFW$Set.Date.and.Time, "%m/%d/%Y")
new_WDFW$year  <- format(new_WDFW$Date_Deployed, "%Y")
new_WDFW$month <- format(new_WDFW$Date_Deployed, "%m")

# trim whitespace; keep ALL trap types
new_WDFW$Trap.Type <- trimws(new_WDFW$Trap.Type)

# subset to June - September in 2023 - 2024
new_WDFW <- new_WDFW[new_WDFW$year %in% c("2023", "2024"), ]
new_WDFW <- new_WDFW[new_WDFW$month %in% c("06", "07", "08", "09"), ]

# clean site names (unchanged from original)
new_WDFW[new_WDFW$Site.Name == " Duckland & 105", "Site.Name"] <- "Duckland"
new_WDFW[new_WDFW$Site.Name == " South", "Site.Name"] <- "South"
new_WDFW[new_WDFW$Site.Name == "West Samish ", "Site.Name"] <- "West Samish"
new_WDFW[new_WDFW$Site.Name == "Ocean Shores Airport", "Site.Name"] <- "Ocean Shores"
new_WDFW[new_WDFW$Site.Name == "Tokeland Hotel", "Site.Name"] <- "Tokeland"
new_WDFW[new_WDFW$Site.Name == "West Samish", "Site.Name"] <- "Samish River"
new_WDFW[new_WDFW$Site.Name == "Pysht ", "Site.Name"] <- "Pysht"
new_WDFW[new_WDFW$Site.Name == "Right Smart Cove ", "Site.Name"] <- "Right Smart Cove"
new_WDFW[new_WDFW$Site.Name == "Point Whitney ", "Site.Name"] <- "Point Whitney"
new_WDFW[new_WDFW$Site.Name == "Shine Tidelands ", "Site.Name"] <- "Shine Tidelands"
new_WDFW[new_WDFW$Site.Name == "Walan Point ", "Site.Name"] <- "Walan Point"
new_WDFW[new_WDFW$Site.Name == "Linger Longer ", "Site.Name"] <- "Linger Longer"
new_WDFW[new_WDFW$Site.Name == "Skokomish Estuary ", "Site.Name"] <- "Skokomish Estuary"
new_WDFW[new_WDFW$Site.Name == "Chico Creek ", "Site.Name"] <- "Chico Creek"
new_WDFW[new_WDFW$Site.Name == "Belfair State Park ", "Site.Name"] <- "Belfair State Park"
new_WDFW[new_WDFW$Site.Name == "Fairmont ", "Site.Name"] <- "Fairmont"
new_WDFW[new_WDFW$Site.Name == "Big Beef Harbor ", "Site.Name"] <- "Big Beef Harbor"
new_WDFW[new_WDFW$Site.Name == "Duckland ", "Site.Name"] <- "Duckland"
new_WDFW[new_WDFW$Site.Name == "Thorndyke ", "Site.Name"] <- "Thorndyke"
new_WDFW[new_WDFW$Site.Name == "Conference Center Lagoon ",
         "Site.Name"] <- "Conference Center Lagoon"
new_WDFW[new_WDFW$Site.Name == "Maynard Lagoon ", "Site.Name"] <- "Maynard Lagoon"

# group by site, year AND trap type
wdfw_old <- data_WDFW %>%
  group_by(year, Site_Name, Trap_Type) %>%
  summarise(catch = sum(as.numeric(CAMA_Total)),
            effort = n(),
            Latitude = mean(Latitude), Longitude = mean(Longitude),
            .groups = "drop") %>%
  mutate(cpue = catch / effort) %>%
  transmute(year, site_name = Site_Name, trap_type = Trap_Type,
            catch, effort, cpue, latitude = Latitude, longitude = Longitude)

wdfw_new <- new_WDFW %>%
  group_by(year, Site.Name, Trap.Type) %>%
  summarise(catch = sum(as.numeric(Total.Catch)),
            effort = n(),
            Latitude = mean(y), Longitude = mean(x),
            .groups = "drop") %>%
  mutate(cpue = catch / effort) %>%
  transmute(year, site_name = Site.Name, trap_type = Trap.Type,
            catch, effort, cpue, latitude = Latitude, longitude = Longitude)

data_WDFW_summary_summer <- rbind(wdfw_old, wdfw_new)

saveRDS(data_WDFW_summary_summer, file.path(procd, "WDFW_summer.rds"))

##################
# Drayton Harbor #
##################

biweek <- c(59, 76, 91, 106, 121, 137, 152, 167, 182, 198,
            213, 229, 244, 259, 274, 290, 305, 320, 335)

DH_counts <- readRDS(file.path(raw, "DraytonHarbor", "counts.rds"))
DH_index  <- readRDS(file.path(raw, "DraytonHarbor", "index.rds"))

# ALL trap types: Fukui (f), Minnow (m), Shrimp (s)
DH_indices <- list(
  Fukui  = readRDS(file.path(raw, "DraytonHarbor", "f_index.rds")),
  Minnow = readRDS(file.path(raw, "DraytonHarbor", "m_index.rds")),
  Shrimp = readRDS(file.path(raw, "DraytonHarbor", "s_index.rds"))
)

years <- c("2020", "2021", "2022", "2023")

DH_rows <- list()
for (tt in names(DH_indices)) {
  idx <- DH_indices[[tt]]
  for (y in seq_along(years)) {
    summer_catch  <- 0
    summer_effort <- 0
    for (t in 1:sum(!is.na(DH_index[, y]))) {
      date <- as.Date(biweek[DH_index[t, y]] - 1,
                      origin = paste0(years[y], "-01-01"))
      if (format(date, "%m") %in% c("06", "07", "08", "09")) {
        summer_catch  <- summer_catch  + sum(idx[t, , y] * rowSums(DH_counts[t, , y, ]))
        summer_effort <- summer_effort + sum(idx[t, , y])
      }
    }
    DH_rows[[length(DH_rows) + 1]] <- data.frame(
      year = years[y], site_name = "Drayton Harbor", trap_type = tt,
      catch = summer_catch, effort = summer_effort,
      cpue = ifelse(summer_effort > 0, summer_catch / summer_effort, NA_real_),
      latitude = 48.98034, longitude = -122.757813,
      stringsAsFactors = FALSE)
  }
}
DH_data_summer <- do.call(rbind, DH_rows)

saveRDS(DH_data_summer, file.path(procd, "DH_summer.rds"))

#######
# WSG #
#######

data_WSG <- read.csv(file.path(raw, "WSG_Traps.csv"))
data_WSG$Date <- as.Date(data_WSG$EndTime, "%m/%d/%Y")

# subset to June - September
data_WSG <- data_WSG[data_WSG$Month %in% c("June", "July",
                                           "August", "September"), ]

# clean site names (unchanged from original)
data_WSG[data_WSG$Site.Name == "Biomonitoring", "Site.Name"] <- "West 90 North"
data_WSG[data_WSG$Site.Name == "Alice Bay", "Site.Name"] <- "Alice Bay Sentinel Site"
data_WSG[data_WSG$Site.Name == "Swinomish Casino", "Site.Name"] <- "Swinomish Casino Marsh"

WSG_sites <- data_WSG %>%
  group_by(Site.Name) %>%
  summarise(Latitude = mean(Latitude), Longitude = mean(Longitude),
            .groups = "drop")

# NOTE: WSG_Traps.csv has no trap-type column; these data were collected with
# Fukui traps and are labelled accordingly.
data_WSG_summary_summer <- data_WSG %>%
  group_by(year, Site.Name) %>%
  summarise(catch = sum(as.numeric(total.cama)),
            effort = n(), .groups = "drop") %>%
  mutate(cpue = catch / effort, trap_type = "Fukui") %>%
  left_join(WSG_sites, by = "Site.Name") %>%
  transmute(year, site_name = Site.Name, trap_type,
            catch, effort, cpue, latitude = Latitude, longitude = Longitude)

saveRDS(data_WSG_summary_summer, file.path(procd, "WSG_summer.rds"))

###############
# Padilla Bay #
###############

data_Padilla <- read.csv(file.path(raw, "PadillaBay.csv"))

data_Padilla$Date <- paste0(match(data_Padilla$Month, month.name),
                            "/", "1/", data_Padilla$Year)
data_Padilla$Date <- as.Date(data_Padilla$Date, "%m/%d/%Y")

# subset to June - September (ALL trap types kept)
data_Padilla <- data_Padilla[data_Padilla$Month %in% c("June", "July",
                                                       "August", "September"), ]

# trim whitespace on trap type
data_Padilla$Trap_Type <- trimws(data_Padilla$Trap_Type)

# clean site names (unchanged from original)
data_Padilla[data_Padilla$Location == "Sullivan Minor and Oswald's",
             "Location"] <- "Sullivan Minor"
data_Padilla[data_Padilla$Location == "Upper Big Indian/Little Indian",
             "Location"] <- "Big Indian Slough"
data_Padilla[data_Padilla$Location == "Joe Leary to West 90 South",
             "Location"] <- "Joe Leary"

Padilla_sites <- data_Padilla %>%
  group_by(Location) %>%
  summarise(Latitude = mean(Latitude), Longitude = mean(Longitude),
            .groups = "drop")

data_Padilla_summary_summer <- data_Padilla %>%
  group_by(Year, Location, Trap_Type) %>%
  summarise(catch = sum(Catch), effort = sum(Effort), .groups = "drop") %>%
  mutate(cpue = catch / effort) %>%
  left_join(Padilla_sites, by = "Location") %>%
  transmute(year = Year, site_name = Location, trap_type = Trap_Type,
            catch, effort, cpue, latitude = Latitude, longitude = Longitude)

saveRDS(data_Padilla_summary_summer, file.path(procd, "Padilla_summer.rds"))

#########
# USFWS / DNWR #
#########

# 2017
DNWR_catch_2017 <- read.csv(file.path(raw, "DNWR", "DNWR_2017_catch2.csv")) %>%
  group_by(TrapID) %>% summarise(total = n(), .groups = "drop")
DNWR_effort_2017 <- read.csv(file.path(raw, "DNWR", "DNWR_2017_effort2.csv"))
DNWR_effort_2017$Date <- as.Date(
  as.character(DNWR_effort_2017$StartLocDate), format = "%m/%d/%Y")
DNWR_effort_2017 <- DNWR_effort_2017 %>%
  mutate(site_name = case_when(Longitude < -123.18 ~ "Dungeness Spit Base",
                               Longitude > -123.16 & Longitude < -123.13 ~
                                 "Graveyard Spit West (Channel)",
                               Longitude > -123.13 ~
                                 "Graveyard Spit East (Lagoon)"))
DNWR_2017 <- left_join(DNWR_effort_2017, DNWR_catch_2017,
                       by = "TrapID") %>% mutate(Year = "2017") %>%
  select(TrapType, Date, Year, total, site_name, Longitude, Latitude) %>%
  rename(Trap.Type = TrapType)

# 2018
DNWR_catch_2018 <- read.csv(file.path(raw, "DNWR", "DNWR_2018_catch.csv")) %>%
  filter(!is.na(Crab..)) %>%
  group_by(TrapIDJoin) %>% summarise(total = n(), .groups = "drop")
DNWR_effort_2018 <- read.csv(file.path(raw, "DNWR", "DNWR_2018_effort.csv"))
DNWR_effort_2018$Date <- as.Date(
  as.character(DNWR_effort_2018$Trap.Pull.Date), format = "%Y%m%d")
DNWR_effort_2018 <- DNWR_effort_2018 %>%
  mutate(site_name = case_when(Longitude < -123.18 ~ "Dungeness Spit Base",
                               Longitude > -123.16 & Longitude < -123.13 ~
                                 "Graveyard Spit West (Channel)",
                               Longitude > -123.13 ~
                                 "Graveyard Spit East (Lagoon)"))
DNWR_2018 <- left_join(DNWR_effort_2018, DNWR_catch_2018,
                       by = "TrapIDJoin") %>% mutate(Year = "2018") %>%
  select(Trap.Type, Date, Year, total, site_name, Longitude, Latitude)

# 2019
DNWR_catch_2019 <- read.csv(file.path(raw, "DNWR", "DNWR_2019_catch.csv")) %>%
  filter(!is.na(Crab..)) %>%
  group_by(TrapIDJoin) %>% summarise(total = n(), .groups = "drop")
DNWR_effort_2019 <- read.csv(file.path(raw, "DNWR", "DNWR_2019_effort.csv"))
DNWR_effort_2019$Date <- as.Date(
  as.character(DNWR_effort_2019$Trap.Pull.Date), format = "%Y%m%d")
DNWR_effort_2019 <- DNWR_effort_2019 %>%
  mutate(site_name = case_when(Longitude < -123.18 ~ "Dungeness Spit Base",
                               Longitude > -123.16 & Longitude < -123.13 ~
                                 "Graveyard Spit West (Channel)",
                               Longitude > -123.13 ~
                                 "Graveyard Spit East (Lagoon)"))
DNWR_2019 <- left_join(DNWR_effort_2019, DNWR_catch_2019,
                       by = "TrapIDJoin") %>% mutate(Year = "2019") %>%
  select(Trap.Type, Date, Year, total, site_name, Longitude, Latitude)

# 2020
DNWR_catch_2020 <- read.csv(file.path(raw, "DNWR", "DNWR_2020_catch.csv")) %>%
  filter(!is.na(Crab..)) %>%
  group_by(TrapIDJoin) %>% summarise(total = n(), .groups = "drop")
DNWR_effort_2020 <- read.csv(file.path(raw, "DNWR", "DNWR_2020_effort.csv"))
DNWR_effort_2020$Date <- as.Date(
  as.character(DNWR_effort_2020$Trap.Pull.Date), format = "%Y%m%d")
DNWR_effort_2020 <- DNWR_effort_2020 %>%
  mutate(site_name = case_when(Longitude < -123.18 ~ "Dungeness Spit Base",
                               Longitude > -123.16 & Longitude < -123.13 ~
                                 "Graveyard Spit West (Channel)",
                               Longitude > -123.13 ~
                                 "Graveyard Spit East (Lagoon)"))
DNWR_2020 <- left_join(DNWR_effort_2020, DNWR_catch_2020,
                       by = "TrapIDJoin") %>% mutate(Year = "2020") %>%
  select(Trap.Type, Date, Year, total, site_name, Longitude, Latitude)

# combine all four years
DNWR <- rbind(DNWR_2017, DNWR_2018, DNWR_2019, DNWR_2020)

# ALL trap types kept (original Fukui-only filter removed); just trim whitespace
DNWR$Trap.Type <- trimws(DNWR$Trap.Type)

# replace NA total with 0
DNWR[is.na(DNWR$total), ]$total <- 0

# subset to just summer
DNWR$month <- format(DNWR$Date, "%m")
DNWR <- DNWR[DNWR$month %in% c("06", "07", "08", "09"), ]

# group by site, year AND trap type
data_DNWR_summary_summer <- DNWR %>%
  group_by(Year, site_name, Trap.Type) %>%
  summarise(catch = sum(total),
            effort = n(),
            latitude = mean(Latitude),
            longitude = mean(Longitude),
            .groups = "drop") %>%
  mutate(cpue = catch / effort) %>%
  transmute(year = Year, site_name, trap_type = Trap.Type,
            catch, effort, cpue, latitude, longitude)

saveRDS(data_DNWR_summary_summer, file.path(procd, "DNWR_summer.rds"))

#######
# DFO #
#######

DFO_early <- read.csv(file.path(raw, "DFO", "allegcdata_2007_2021.csv"))
DFO_late  <- read.csv(file.path(raw, "DFO", "allegcdata_2022_2024.csv"))
DFO_all   <- rbind(DFO_early, DFO_late)

# Keep ALL gear/trap types (original kept only GearCode == "51", i.e. Fukui).
# Still restrict to usable traps (unchanged from original).
DFO_all <- DFO_all %>% filter(TrapUsability %in% c(0, 15))

# trap-type label from gear code (51 = Fukui; other codes unknown here)
DFO_all <- DFO_all %>%
  mutate(trap_type = ifelse(as.character(GearCode) == "51", "Fukui",
                            paste0("GearCode_", GearCode)))

# unique trap sets (now also grouped by gear code / trap type)
unique_traps <- DFO_all %>%
  group_by(SetNum, Year, Month, Day, HoursSoak, GeogLoc, startLAT, startLONG,
           endLAT, endLONG, TrapSpacing, TrapNum, TrapUsability,
           GearCode, trap_type) %>%
  summarize(nind = n(), .groups = "drop") %>%
  unite("TrapID", GeogLoc, SetNum, Month, Day, TrapNum,
        sep = "_", remove = FALSE)

# EGC captures (catch still defined as Species == "XMB")
EGC_captures <- DFO_all %>%
  filter(Species == "XMB") %>%
  unite("TrapID", GeogLoc, SetNum, Month, Day, TrapNum,
        sep = "_", remove = FALSE) %>%
  group_by(TrapID, GeogLoc, SetNum, Year, Month, Day, TrapNum) %>%
  summarize(count = n(), .groups = "drop")

# add captures to effort
unique_traps <- left_join(unique_traps, EGC_captures,
                          by = c("TrapID", "Year", "Month", "Day",
                                 "SetNum", "GeogLoc", "TrapNum"))

# replace NA count with 0
unique_traps[is.na(unique_traps$count), "count"] <- 0

# subset to summer, summarise per site-year-traptype
DFO_final <- unique_traps[unique_traps$Month %in% c(6, 7, 8, 9), ] %>%
  group_by(Year, GeogLoc, trap_type) %>%
  summarise(catch = sum(count),
            effort = n(),
            latitude = mean(startLAT),
            longitude = mean(startLONG),
            .groups = "drop") %>%
  mutate(cpue = catch / effort) %>%
  transmute(year = Year, site_name = GeogLoc, trap_type,
            catch, effort, cpue, latitude, longitude)

saveRDS(DFO_final, file.path(procd, "DFO_summer.rds"))

############################
# combine all data sources #
############################

tag <- function(df, src) {
  # year is character in some sources and integer in others; unify to integer
  df$year <- suppressWarnings(as.integer(as.character(df$year)))
  dplyr::mutate(df, data_source = src)
}

all_sources <- dplyr::bind_rows(
  tag(data_WDFW_summary_summer,   "WDFW"),
  tag(DH_data_summer,             "NWStraits"),
  tag(data_WSG_summary_summer,    "WSG"),
  tag(data_Padilla_summary_summer,"Padilla"),
  tag(data_DNWR_summary_summer,   "DNWR"),
  tag(DFO_final,                  "DFO"),
  tag(summary_Makah,              "Makah")
)
all_sources$year <- as.integer(all_sources$year)

# ---------------------------------------------------------------------------
# Normalize trap-type labels.
#
# Sources spell the same trap several different ways. This function reconciles
# differences that are PURELY casing, whitespace, punctuation or word order --
# i.e. cases where the two strings are obviously the same physical trap. It
# does NOT merge labels that differ in substance.
#
#   casing        "FUKUI" / "SHRIMP" / "MINNOW"  -> "Fukui" / "Shrimp" / "Minnow"
#                 (Makah writes trap types in all caps)
#   dashes        en/em dashes -> plain hyphen
#                 "Shrimp - Other" (en dash)     -> "Shrimp - Other"
#   hyphen space  uniform " - " spacing
#                 "Shrimp-C"                     -> "Shrimp - C"   [merges with
#                                                   the existing "Shrimp - C"]
#   word order    '1/2" Shrimp'                  -> 'Shrimp 1/2"'  [merges with
#                                                   the existing 'Shrimp 1/2"']
#
# DELIBERATELY NOT MERGED (these are, or may be, genuinely different traps --
# collapsing them would destroy real information):
#   * "Shrimp" vs "Shrimp - C" / "- R" / "- K"   generic vs. specific variant
#   * "Shrimp - C" vs "Collapsible Shrimp"       C is unconfirmed; may be a
#                                                crew/site code, not "collapsible"
#   * "Minnow" vs "Minnow - Unmodified"          the qualifier implies a
#                                                modified counterpart exists
#   * "Pitfall" vs "Pitfall - Small" / "- Large" size classes
#   * "GearCode_*"                               DFO codes; left untouched until
#                                                mapped from DFO's lookup table
# ---------------------------------------------------------------------------
normalize_trap_type <- function(x) {
  x    <- trimws(x)
  keep <- grepl("^GearCode_", x)   # DFO gear codes pass through unchanged
  y    <- x[!keep]

  y <- str_to_title(y)                       # casing
  y <- gsub("‒|–|—|―", "-", y)  # en/em dash -> hyphen
  y <- gsub("\\s+", " ", y)                  # collapse internal whitespace
  y <- gsub("\\s*-\\s*", " - ", y)           # uniform " - " around hyphens
  y <- trimws(y)
  y <- ifelse(y == '1/2" Shrimp', 'Shrimp 1/2"', y)  # word-order variant

  x[!keep] <- y
  x
}
all_sources$trap_type <- normalize_trap_type(all_sources$trap_type)

# ---------------------------------------------------------------------------
# *** LOSSY STEP - SET THE FLAG BELOW TO FALSE TO UNDO ***
#
# Collapse the Pitfall, Minnow and Shrimp variants into three parent labels.
# Unlike normalize_trap_type() above, this DISCARDS REAL INFORMATION: the
# variant suffixes are thrown away and their catch/effort are pooled into the
# parent trap type. It is applied because the variants are too sparse to
# estimate separate trap-type effects (most have < 25 records), not because
# they are known to be the same trap.
#
#   Pitfall  <- "Pitfall - Small", "Pitfall - Large"
#   Minnow   <- "Minnow - Unmodified"
#   Shrimp   <- "Shrimp - C", "Shrimp - K", "Shrimp - R", "Shrimp - Other",
#               'Shrimp 1/2"', "Collapsible Shrimp"
#
# NOT collapsed: Fukui, Hand, Burrow, Russell, Crayfish, Seine Net,
# Holding Cage, GearCode_*.
#
# TO UNDO: set collapse_trap_variants <- FALSE and re-run this script, then
# re-run Code/BayesianTrapData.R. Nothing downstream hard-codes the collapsed
# labels, so that is the only change required. Note that the unmodified
# variant labels are still preserved in the per-source .rds files written to
# Data/Processed by this script, so nothing is lost on disk either way.
# ---------------------------------------------------------------------------
collapse_trap_variants <- TRUE

if (collapse_trap_variants) {
  collapse_trap_type <- function(x) {
    dplyr::case_when(
      grepl("^GearCode_", x)   ~ x,          # DFO codes untouched
      grepl("Shrimp",  x, ignore.case = TRUE) ~ "Shrimp",
      grepl("^Minnow", x, ignore.case = TRUE) ~ "Minnow",
      grepl("^Pitfall", x, ignore.case = TRUE) ~ "Pitfall",
      TRUE ~ x
    )
  }
  all_sources$trap_type <- collapse_trap_type(all_sources$trap_type)
}

# one representative coordinate per site: mean lat/lon across ALL of that site's
# records, so every site has a single X,Y in the output regardless of year,
# trap type, or data source.
site_coords <- all_sources %>%
  group_by(site_name) %>%
  summarise(latitude = mean(latitude, na.rm = TRUE),
            longitude = mean(longitude, na.rm = TRUE),
            .groups = "drop")

# summarize duplicate sites across data sources (now also per trap type)
all_sites <- unique(all_sources[, c("site_name", "data_source",
                                    "year", "trap_type")]) %>%
  group_by(site_name, year, trap_type) %>%
  summarise(data_sources = paste(data_source, collapse = ", "),
            .groups = "drop")

all_sources_final <- left_join(all_sources, all_sites,
                               by = c("site_name", "year", "trap_type")) %>%
  group_by(site_name, year, trap_type, data_sources) %>%
  summarise(catch = sum(catch),
            effort = sum(effort),
            cpue = catch / effort,
            .groups = "drop") %>%
  left_join(site_coords, by = "site_name")  # single X,Y per site

saveRDS(all_sources_final, file.path(procd, "all_sources_summer.rds"))
write.csv(all_sources_final,
          here::here("data", "data_all_sources_all_traptypes.csv"),
          row.names = FALSE)
