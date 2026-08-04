################################################################################
# Inspect the replicate-level model inputs produced by
# code/data_prep/02_BayesianTrapData.R.
################################################################################

project_library <- file.path(getwd(), ".r-library")
if (dir.exists(project_library)) {
  .libPaths(c(normalizePath(project_library), .libPaths()))
}

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required.")
}
here::i_am("sandbox.R")

model_dir <- here::here("data", "model_data")
presence <- readRDS(file.path(model_dir, "PresenceArrayBinary.rds"))
traptype <- readRDS(file.path(model_dir, "TrapTypeArraysNew.rds"))
cpue_fukui <- readRDS(file.path(model_dir, "cpue_fukui.rds"))

check <- function(description, condition) {
  passed <- isTRUE(condition)
  cat(sprintf("[%s] %s\n", if (passed) "PASS" else "FAIL", description))
  if (!passed) stop("Validation failed: ", description, call. = FALSE)
  invisible(passed)
}

used <- !is.na(presence)
n_used <- apply(used, c(1, 2), sum)

cat("Replicate-level model products\n")
cat("Dimensions: ", paste(names(dimnames(presence)), dim(presence),
                           sep = "=", collapse = " x "), "\n", sep = "")

check("Presence is [site, year, replicate]",
      identical(names(dimnames(presence)), c("site", "year", "replicate")))
check("Presence contains only binary values in used cells",
      all(presence[used] %in% c(0L, 1L)))
check("Replicate slots form a contiguous prefix in every site-year",
      all(apply(used, c(1, 2), function(z) all(diff(as.integer(z)) <= 0))))
check("Trap type is a named list for the three canonical types",
      identical(names(traptype), c("Fukui", "Shrimp", "Minnow")))
check("Trap-type arrays have the same dimensions, labels, and NA pattern",
      all(vapply(traptype, function(a) {
        identical(dimnames(a), dimnames(presence)) && identical(!is.na(a), used)
      }, logical(1))))
check("Each real replicate has exactly one trap type",
      all(Reduce(`+`, traptype)[used] == 1L))
check("Fukui CPUE aligns with the site and year axes",
      identical(dimnames(presence)[c("site", "year")], dimnames(cpue_fukui)))
check("Fukui CPUE exists exactly where Fukui replicates were sampled",
      all(unname(!is.na(cpue_fukui)) ==
            unname(apply(traptype$Fukui == 1L, c(1, 2), any, na.rm = TRUE))))
check("All years are 2018 or later",
      all(as.integer(dimnames(presence)[["year"]]) >= 2018L))

cat("\nContent summary\n")
cat("Sites: ", dim(presence)[1], "\n", sep = "")
cat("Years: ", paste(dimnames(presence)[["year"]], collapse = ", "), "\n", sep = "")
cat("Maximum replicates in a site-year: ", dim(presence)[3], "\n", sep = "")
cat("Actual trap replicates: ", sum(used), "\n", sep = "")
cat("Detection-positive replicates: ", sum(presence, na.rm = TRUE), "\n", sep = "")
cat("Replicates by trap type:\n")
print(vapply(traptype, function(a) sum(a, na.rm = TRUE), integer(1)))
cat("Replicates per sampled site-year:\n")
print(summary(n_used[n_used > 0L]))
