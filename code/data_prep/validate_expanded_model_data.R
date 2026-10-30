################################################################################
# Validate the expanded replicate-level occupancy products against the original
# in-extent products and report noteworthy (not necessarily erroneous) values.
################################################################################

model_dir <- file.path("data", "model_data")
read_model <- function(name) readRDS(file.path(model_dir, name))

original_presence <- read_model("PresenceArrayBinary.rds")
expanded_presence <- read_model("PresenceArrayBinary_expanded.rds")
original_traps <- read_model("TrapTypeArraysNew.rds")
expanded_traps <- read_model("TrapTypeArraysNew_expanded.rds")
original_cpue <- read_model("cpue_fukui.rds")
expanded_cpue <- read_model("cpue_fukui_expanded.rds")
original_wsg <- read_model("PresenceArrayBinary_WSG.rds")
expanded_wsg <- read_model("PresenceArrayBinary_WSG_expanded.rds")

checks <- data.frame(check = character(), passed = logical(), detail = character())
add_check <- function(name, passed, detail) {
  checks <<- rbind(checks, data.frame(
    check = name, passed = isTRUE(passed), detail = detail,
    stringsAsFactors = FALSE
  ))
}

original_sites <- dimnames(original_presence)$site
expanded_sites <- dimnames(expanded_presence)$site
site_index <- match(original_sites, expanded_sites)
original_reps <- dimnames(original_presence)$replicate

add_check(
  "All original individual-trap sites occur in expanded data",
  !anyNA(site_index),
  sprintf("%d original; %d expanded; %d added", length(original_sites),
          length(expanded_sites), length(setdiff(expanded_sites, original_sites)))
)
add_check(
  "Original presence values are unchanged in expanded data",
  identical(
    original_presence,
    expanded_presence[site_index, dimnames(original_presence)$year,
                      original_reps, drop = FALSE]
  ),
  "Exact equality including NA pattern and dimension names"
)
add_check(
  "Original Fukui CPUE values are unchanged in expanded data",
  identical(
    original_cpue,
    expanded_cpue[site_index, dimnames(original_cpue)[[2]], drop = FALSE]
  ),
  "Exact equality including NA pattern and dimension names"
)

for (trap_name in names(expanded_traps)) {
  add_check(
    paste("Original", trap_name, "indicator values are unchanged"),
    identical(
      original_traps[[trap_name]],
      expanded_traps[[trap_name]][
        site_index, dimnames(original_presence)$year, original_reps, drop = FALSE
      ]
    ),
    "Exact equality including NA pattern and dimension names"
  )
}

used <- !is.na(expanded_presence)
trap_used <- lapply(expanded_traps, function(x) !is.na(x))
trap_sum <- Reduce(`+`, expanded_traps)
add_check(
  "Expanded presence is binary",
  all(expanded_presence[used] %in% 0:1),
  sprintf("%d observed trap replicates", sum(used))
)
add_check(
  "Trap indicator NA patterns match expanded presence",
  all(vapply(trap_used, identical, logical(1), y = used)),
  "Fukui, Minnow, and Shrimp checked"
)
add_check(
  "Exactly one trap type is active per observed replicate",
  all(trap_sum[used] == 1L),
  sprintf("%d observed trap replicates checked", sum(used))
)
add_check(
  "Fukui CPUE is finite and nonnegative when present",
  all(is.finite(expanded_cpue[!is.na(expanded_cpue)])) &&
    all(expanded_cpue[!is.na(expanded_cpue)] >= 0),
  sprintf("%d sampled site-years", sum(!is.na(expanded_cpue)))
)
add_check(
  "Fukui CPUE availability matches Fukui sampling",
  identical(
    unname(!is.na(expanded_cpue)),
    unname(apply(expanded_traps$Fukui == 1L, c(1, 2), any, na.rm = TRUE))
  ),
  "Every and only site-years with Fukui deployments have CPUE"
)

wsg_index <- match(dimnames(original_wsg)$site, dimnames(expanded_wsg)$site)
add_check(
  "Original WSG sites and values are unchanged",
  !anyNA(wsg_index) && identical(
    original_wsg, expanded_wsg[wsg_index, , , drop = FALSE]
  ),
  sprintf("%d original; %d expanded; %d added", dim(original_wsg)[1],
          dim(expanded_wsg)[1], dim(expanded_wsg)[1] - dim(original_wsg)[1])
)

coordinates <- read.csv(file.path(model_dir, "site_coordinates_expanded.csv"),
                        stringsAsFactors = FALSE)
add_check(
  "Every expanded site has one valid coordinate pair",
  !anyDuplicated(coordinates$site) &&
    all(is.finite(coordinates$x)) && all(is.finite(coordinates$y)) &&
    setequal(coordinates$site, c(expanded_sites, dimnames(expanded_wsg)$site)),
  sprintf("%d source-specific sites", nrow(coordinates))
)
add_check(
  "Explicitly excluded DFO possible-GPS-error site is absent",
  !"DFO::possibleGPSerror" %in% coordinates$site &&
    !"DFO::possibleGPSerror" %in% expanded_sites,
  "Checked coordinate table and expanded individual-trap site axis"
)

zone_map <- read.csv(file.path(model_dir, "site_zone_map_expanded.csv"),
                     stringsAsFactors = FALSE)
wsg_zone_map <- read.csv(file.path(model_dir, "site_zone_map_WSG_expanded.csv"),
                         stringsAsFactors = FALSE)
add_check(
  "Expanded zone-map rows align with array site axes",
  identical(zone_map$site, expanded_sites) &&
    identical(wsg_zone_map$site, dimnames(expanded_wsg)$site),
  "Order and identifiers checked"
)

write.csv(checks, file.path(model_dir, "expanded_data_validation.csv"),
          row.names = FALSE)

replicate_counts <- apply(!is.na(expanded_presence), c(1, 2), sum)
max_reps <- max(replicate_counts)
max_reps_at <- which(replicate_counts == max_reps, arr.ind = TRUE)
max_cpue <- max(expanded_cpue, na.rm = TRUE)
max_cpue_at <- which(expanded_cpue == max_cpue, arr.ind = TRUE)
expanded_events <- read_model("SamplingEvents_WSG_expanded.rds")

flags <- data.frame(
  flag = c(
    "Very large replicate dimension",
    "High Fukui CPUE",
    "Individual-trap sites outside current zone polygons",
    "WSG sites outside current zone polygons",
    "WSG events with effort other than six"
  ),
  detail = c(
    sprintf("Maximum is %d traps at %s in %s; this sets the third array dimension for every site-year.",
            max_reps, rownames(replicate_counts)[max_reps_at[1, 1]],
            colnames(replicate_counts)[max_reps_at[1, 2]]),
    sprintf("Maximum is %.5f at %s in %s.", max_cpue,
            rownames(expanded_cpue)[max_cpue_at[1, 1]],
            colnames(expanded_cpue)[max_cpue_at[1, 2]]),
    sprintf("%d of %d sites have no zone_id.", sum(is.na(zone_map$zone_id)),
            nrow(zone_map)),
    sprintf("%d of %d sites have no zone_id.", sum(is.na(wsg_zone_map$zone_id)),
            nrow(wsg_zone_map)),
    sprintf("%d events; observed effort ranges from %g to %g traps.",
            sum(expanded_events$effort != 6), min(expanded_events$effort),
            max(expanded_events$effort))
  ),
  stringsAsFactors = FALSE
)
write.csv(flags, file.path(model_dir, "expanded_data_flags.csv"), row.names = FALSE)

if (!all(checks$passed)) {
  stop("Expanded-data validation failed: ",
       paste(checks$check[!checks$passed], collapse = "; "))
}
message("All ", nrow(checks), " expanded-data validation checks passed.")
message("Wrote expanded_data_validation.csv and expanded_data_flags.csv.")
