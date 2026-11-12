################################################################################
# Compare modeled larval supply with apparent colonizations at different lags.
#
# For lag L, larvae released in year Y are matched to an apparent colonization
# detected in year Y + L: no detection at the site in Y + L - 1, and at least
# one detection in Y + L. Source zones count as "occupied" if EGC was detected
# there in year Y (or in any year up to Y, as a looser proxy).
#
# This uses raw detections (no correction for imperfect detection), so it is
# descriptive only.
#
# Usage:  Rscript code/explore_colonization_lag.R 1 2 3
################################################################################

lags <- as.integer(commandArgs(trailingOnly = TRUE))
if (length(lags) == 0) lags <- 1:3

#############
# larval data (2013 - 2025)
#############

conn_years <- 2013:2025
conn_paths <- sprintf("data/connectivity/%d_connectivity_counts_365.csv", conn_years)
conn_2013 <- read.csv(conn_paths[1])
site_names <- conn_2013$settlement_site_id
nsite <- length(site_names)
ncy <- length(conn_years)

larv_S <- array(NA, dim = c(nsite, nsite, ncy))
for (k in seq_along(conn_paths)) {
  larv_S[, , k] <- as.matrix(read.csv(conn_paths[k])[, 3:(nsite + 2)])
}
larv_R <- as.matrix(read.csv("data/connectivity/larvae_released.csv",
                             row.names = 1)[, 2:(ncy + 1)])
larv_E <- as.matrix(read.csv("data/connectivity/larvae_external_settlers.csv",
                             row.names = 1)[, 2:(ncy + 1)])
stopifnot(all(rownames(larv_R) == site_names), all(rownames(larv_E) == site_names))

# per-source dispersal probability, own zone excluded (unscaled)
P_int <- larv_S
for (k in 1:ncy) {
  R <- larv_R[, k]
  P_int[, , k] <- sweep(larv_S[, , k], 2, ifelse(R > 0, R, 1), "/")
  P_int[, R == 0, k] <- 0
  diag(P_int[, , k]) <- 0
}
int_all <- apply(P_int, c(1, 3), sum)        # supply if every zone were occupied
ever_src <- apply(P_int > 0, 1, any)         # zone receives internal larvae in any year

#############
# observed occupancy (2014 - 2024), main traps + WSG
#############

detections <- readRDS("data/model_data/PresenceArrayBinary_expanded.rds")
detections_WSG <- readRDS("data/model_data/PresenceArrayBinary_WSG_expanded.rds")
type <- readRDS("data/model_data/TrapTypeArraysNew_expanded.rds")
obs_years <- as.integer(dimnames(detections)[[2]])
nobs <- length(obs_years)

keep <- rownames(detections[, 1, ]) %in% site_names
detections <- detections[keep, , ]
types <- lapply(type[c("Minnow", "Fukui", "Shrimp")], function(x) x[keep, , ])
larv_ind <- match(rownames(detections[, 1, ]), site_names)

# a site-year counts as surveyed only if it has traps of a modeled type,
# matching how the model uses the data
obs <- matrix(NA, nsite, nobs)
any_trap <- Reduce(`|`, lapply(types, function(x) !is.na(detections) & !is.na(x) & x == 1))
n_tr <- rowSums(any_trap, dims = 2)
y_tr <- rowSums(any_trap & detections == 1, dims = 2)
obs[larv_ind, ][n_tr > 0] <- (y_tr > 0)[n_tr > 0]

wsg_ind <- match(rownames(detections_WSG), site_names)
n_w <- rowSums(!is.na(detections_WSG), dims = 2)
y_w <- rowSums(detections_WSG, dims = 2, na.rm = TRUE)
for (r in seq_along(wsg_ind)) for (t in 1:nobs) if (n_w[r, t] > 0) {
  obs[wsg_ind[r], t] <- max(obs[wsg_ind[r], t], y_w[r, t] > 0, na.rm = TRUE)
}

det <- matrix(obs %in% 1, nsite, nobs)                 # detected in year
ever_det <- t(apply(det, 1, cumsum)) > 0               # detected in any year up to this one

#############
# helpers
#############

rate_tab <- function(x, g) {
  out <- do.call(rbind, tapply(x, g, function(v) c(n = length(v), col = sum(v),
                                                   rate = round(mean(v), 3))))
  out[!is.na(out[, 1]), , drop = FALSE]
}
cut_nonzero <- function(x) {
  q <- quantile(x[x > 0], c(1/3, 2/3))
  if (sum(x > 0) < 3 || any(duplicated(q))) return(factor(ifelse(x > 0, ">0", "0")))
  cut(x, c(-Inf, 0, q, Inf), labels = c("0", "low", "mid", "high"))
}

#############
# analysis per lag
#############

for (L in lags) {
  rows <- list()
  for (Y in conn_years) {
    t_now <- match(Y + L, obs_years); t_prev <- match(Y + L - 1, obs_years)
    t_src <- match(Y, obs_years); k <- match(Y, conn_years)
    if (is.na(t_now) || is.na(t_prev)) next
    src_det  <- if (is.na(t_src)) rep(0, nsite) else det[, t_src]
    src_ever <- if (is.na(t_src)) rep(0, nsite) else ever_det[, t_src]
    rows[[length(rows) + 1]] <- data.frame(
      i = 1:nsite, release = Y, detect = Y + L,
      prev = obs[, t_prev], now = obs[, t_now],
      int_all = int_all[, k],
      conn_det = c(P_int[, , k] %*% src_det),
      conn_ever = c(P_int[, , k] %*% src_ever),
      E = larv_E[, k])
  }
  d <- do.call(rbind, rows)
  dc <- d[d$prev %in% 0 & !is.na(d$now), ]
  dc$col <- as.integer(dc$now)
  dc$ever_src <- ever_src[dc$i]
  dc$program <- sub("::.*", "", site_names[dc$i])

  cat("\n################################################################\n")
  cat(sprintf("LAG %d YEAR(S): larvae released %d-%d, colonizations detected %d-%d\n",
              L, min(dc$release), max(dc$release), min(dc$detect), max(dc$detect)))
  cat(sprintf("site-years undetected -> surveyed: %d, apparent colonizations: %d\n",
              nrow(dc), sum(dc$col)))

  cat("\n-- larval supply at colonized vs not\n")
  f <- function(x) c(n = nrow(x),
                     no_internal = sum(x$int_all == 0),
                     never_internal = sum(!x$ever_src),
                     no_external = sum(x$E == 0),
                     no_supply = sum(x$int_all == 0 & x$E == 0))
  print(rbind(colonized = f(dc[dc$col == 1, ]), not = f(dc[dc$col == 0, ])))
  cat("\n-- colonization rate: any modeled supply (internal or external)?\n")
  print(rate_tab(dc$col, ifelse(dc$int_all > 0 | dc$E > 0, "some supply", "no supply")))

  cat("\n-- by internal supply from sources detected in release year\n")
  print(rate_tab(dc$col, cut_nonzero(dc$conn_det)))
  cat("\n-- by internal supply from sources detected in any year up to release\n")
  print(rate_tab(dc$col, cut_nonzero(dc$conn_ever)))
  cat("\n-- by external settlers in release year\n")
  print(rate_tab(dc$col, cut(dc$E, c(-1, 0, 5, 20, Inf), labels = c("0", "1-5", "6-20", ">20"))))

  cat("\n-- naive logistic regressions (log(x + offset) covariates)\n")
  fit <- glm(col ~ log(conn_ever + 0.001) + log1p(E), binomial, dc)
  print(round(summary(fit)$coef, 3))

  cat("\n-- colonized site-years\n")
  cz <- dc[dc$col == 1, ]
  cz$site <- site_names[cz$i]
  print(cz[order(cz$site, cz$detect),
           c("site", "release", "detect", "int_all", "conn_det", "conn_ever", "E")],
        digits = 2, row.names = FALSE)
}
