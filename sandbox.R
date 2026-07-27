################################################################################
# sandbox.R
#
# Inspect and cross-check the three model-data products created by
# code/data_prep/02_BayesianTrapData.R.
#
# Run from RStudio with Source, or from the repository root with:
#   Rscript sandbox.R
################################################################################

options(stringsAsFactors = FALSE)

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it with install.packages('here').")
}
here::i_am("sandbox.R")

model_dir <- here::here("data", "model_data")
files <- c(
  presence = file.path(model_dir, "PresenceArrayBinary.rds"),
  traptype = file.path(model_dir, "TrapTypeArraysNew.rds"),
  cpue = file.path(model_dir, "cpue_fukui.rds")
)

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing model product(s):\n  ",
    paste(missing_files, collapse = "\n  "),
    "\nRun code/data_prep/02_BayesianTrapData.R first."
  )
}

presence <- readRDS(files[["presence"]])
traptype <- readRDS(files[["traptype"]])
cpue_fukui <- readRDS(files[["cpue"]])

section <- function(title) {
  cat("\n", paste(rep("=", 78L), collapse = ""), "\n", sep = "")
  cat(title, "\n")
  cat(paste(rep("=", 78L), collapse = ""), "\n", sep = "")
}

check <- function(description, condition) {
  passed <- isTRUE(condition)
  cat(sprintf("[%s] %s\n", if (passed) "PASS" else "FAIL", description))
  passed
}

format_dimensions <- function(x) {
  if (is.null(dim(x))) return("<no dimensions>")
  labels <- names(dimnames(x))
  if (is.null(labels)) labels <- rep("unnamed", length(dim(x)))
  paste(paste0(labels, "=", dim(x)), collapse = " x ")
}

section("1. OBJECT STRUCTURE")
cat("Presence:  ", format_dimensions(presence), "\n", sep = "")
cat("Trap type: ", format_dimensions(traptype), "\n", sep = "")
cat("Fukui CPUE:", format_dimensions(cpue_fukui), "\n", sep = "")
cat("\nObject classes:\n")
cat("  Presence:  ", paste(class(presence), collapse = ", "), "\n", sep = "")
cat("  Trap type: ", paste(class(traptype), collapse = ", "), "\n", sep = "")
cat("  Fukui CPUE:", paste(class(cpue_fukui), collapse = ", "), "\n", sep = "")

section("2. ALIGNMENT CHECKS")

checks <- c(
  check(
    "Presence is a three-dimensional array",
    is.array(presence) && length(dim(presence)) == 3L
  ),
  check(
    "Trap-type deployment is a three-dimensional array",
    is.array(traptype) && length(dim(traptype)) == 3L
  ),
  check("Fukui CPUE is a matrix", is.matrix(cpue_fukui))
)

new_format <- all(checks)

if (new_format) {
  presence_names <- names(dimnames(presence))
  traptype_names <- names(dimnames(traptype))
  cpue_names <- names(dimnames(cpue_fukui))

  checks <- c(
    checks,
    check(
      "Presence dimensions are named site, trap_type, year",
      identical(presence_names, c("site", "trap_type", "year"))
    ),
    check(
      "Trap-type dimensions are named site, trap_type, year",
      identical(traptype_names, c("site", "trap_type", "year"))
    ),
    check(
      "CPUE dimensions are named site, year",
      identical(cpue_names, c("site", "year"))
    ),
    check(
      "Presence and trap-type dimensions and labels are identical",
      identical(dim(presence), dim(traptype)) &&
        identical(dimnames(presence), dimnames(traptype))
    ),
    check(
      "Presence and CPUE site labels and ordering are identical",
      identical(dimnames(presence)[["site"]], dimnames(cpue_fukui)[["site"]])
    ),
    check(
      "Presence and CPUE year labels and ordering are identical",
      identical(dimnames(presence)[["year"]], dimnames(cpue_fukui)[["year"]])
    ),
    check(
      "Only Shrimp, Fukui, and Minnow are present",
      identical(
        dimnames(presence)[["trap_type"]],
        c("Shrimp", "Fukui", "Minnow")
      )
    ),
    check(
      "Presence contains only zero and one",
      all(presence %in% c(0L, 1L))
    ),
    check(
      "Deployment contains only zero and one",
      all(traptype %in% c(0L, 1L))
    ),
    check(
      "Every detection has corresponding sampling effort",
      all(presence <= traptype)
    ),
    check(
      "All non-missing CPUE values are finite and non-negative",
      all(is.na(cpue_fukui) |
            (is.finite(cpue_fukui) & cpue_fukui >= 0))
    )
  )

  site_names <- dimnames(presence)[["site"]]
  year_names <- dimnames(presence)[["year"]]
  trap_names <- dimnames(presence)[["trap_type"]]

  duplicate_sites <- unique(site_names[
    duplicated(tolower(trimws(site_names))) |
      duplicated(tolower(trimws(site_names)), fromLast = TRUE)
  ])

  checks <- c(
    checks,
    check(
      "No duplicate site names after trimming whitespace and ignoring case",
      length(duplicate_sites) == 0L
    ),
    check(
      "All included years are 2018 or later",
      all(suppressWarnings(as.integer(year_names)) >= 2018L)
    )
  )

  fukui_deployed <- traptype[, "Fukui", , drop = TRUE] == 1L
  checks <- c(
    checks,
    check(
      "CPUE is present exactly when Fukui traps were deployed",
      identical(dim(fukui_deployed), dim(cpue_fukui)) &&
        identical(dimnames(fukui_deployed), dimnames(cpue_fukui)) &&
        identical(fukui_deployed, !is.na(cpue_fukui))
    )
  )

  if (length(duplicate_sites) > 0L) {
    cat("\nDuplicate-looking sites:\n")
    print(duplicate_sites)
  }

  section("3. CONTENT SUMMARY")
  cat("Sites:      ", length(site_names), "\n", sep = "")
  cat("Years:      ", paste(year_names, collapse = ", "), "\n", sep = "")
  cat("Trap types: ", paste(trap_names, collapse = ", "), "\n", sep = "")

  summary_table <- data.frame(
    trap_type = trap_names,
    sampled_cells = apply(traptype, 2L, sum),
    detection_cells = apply(presence, 2L, sum),
    detection_rate_when_sampled = vapply(
      seq_along(trap_names),
      function(i) {
        sampled <- traptype[, i, ] == 1L
        if (!any(sampled)) return(NA_real_)
        mean(presence[, i, ][sampled])
      },
      numeric(1L)
    ),
    row.names = NULL
  )
  print(summary_table)

  cat("\nFukui CPUE summary:\n")
  print(summary(as.vector(cpue_fukui)))

  cat("\nSites with no sampling by any of the three traps: ")
  sampled_by_site <- apply(traptype, 1L, sum)
  cat(sum(sampled_by_site == 0L), "\n", sep = "")

  cat("Sites with at least one detection: ")
  detected_by_site <- apply(presence, 1L, sum)
  cat(sum(detected_by_site > 0L), "\n", sep = "")

  section("4. FIRST FEW LABELS AND VALUES")
  cat("First 20 sites:\n")
  print(utils::head(site_names, 20L))

  cat("\nPresence array for the first three sites:\n")
  print(presence[seq_len(min(3L, length(site_names))), , , drop = FALSE])

  cat("\nDeployment array for the first three sites:\n")
  print(traptype[seq_len(min(3L, length(site_names))), , , drop = FALSE])

  cat("\nFukui CPUE for the first ten sites:\n")
  print(cpue_fukui[seq_len(min(10L, length(site_names))), , drop = FALSE])

  section("5. OPTIONAL DIAGNOSTIC PLOTS")
  if (interactive()) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mfrow = c(1L, 2L), mar = c(8, 4, 3, 1))

    barplot(
      summary_table$sampled_cells,
      names.arg = summary_table$trap_type,
      las = 2L,
      ylab = "Sampled site-years",
      main = "Sampling effort"
    )

    cpue_values <- cpue_fukui[!is.na(cpue_fukui)]
    if (length(cpue_values) > 0L) {
      hist(
        cpue_values,
        breaks = "FD",
        xlab = "Fukui CPUE",
        main = "Fukui CPUE distribution"
      )
    } else {
      plot.new()
      title("No non-missing Fukui CPUE")
    }
    cat("Two diagnostic plots were sent to the RStudio Plots pane.\n")
  } else {
    cat("Run with Source in RStudio to display the diagnostic plots.\n")
  }

  section("FINAL RESULT")
  if (all(checks)) {
    cat("ALL CHECKS PASSED: the three products are structurally aligned.\n")
  } else {
    cat("ONE OR MORE CHECKS FAILED. Review each [FAIL] line above.\n")
  }
} else {
  section("FINAL RESULT")
  cat(
    "These files use an older or unexpected structure, so the detailed\n",
    "alignment checks could not run. Regenerate them with:\n",
    "  code/data_prep/02_BayesianTrapData.R\n",
    sep = ""
  )
}
