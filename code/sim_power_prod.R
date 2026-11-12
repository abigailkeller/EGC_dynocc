################################################################################
# Simulation check for the prod model: can beta1 (internal connectivity from
# occupied sources) and beta2 (external settlers) be recovered given the real
# sites, connectivity (P_src), survey design (agg / agg_w) and ~10-11 years?
#
# For each scenario (true beta1, beta2) and replicate:
#   1. simulate z and detections from the prod model at the true parameter values
#      (same sites, same traps per site-year and trap type, same WSG effort)
#   2. refit the prod model to the simulated detections
#   3. record posterior median / 95% CI for each parameter
#
# Usage:  Rscript code/sim_power_prod.R 1year    (or 2year)
################################################################################

library(nimble)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
lag <- if (length(args) > 0) args[1] else "1year"
prod_file <- sprintf("code/main_rescale_%s_prod.R", lag)

################
# settings     #
################

# true values for the parameters not being varied
# (posterior medians from posterior_samples_20261111_<lag>_prod.rds)
truth_fixed <- list(
  "1year" = list(psi = 0.405, p_minnow = 0.206, p_fukui = 0.291, p_shrimp = 0.644,
                 beta0 = -2.506, beta3 = 1.860),
  "2year" = list(psi = 0.363, p_minnow = 0.206, p_fukui = 0.281, p_shrimp = 0.644,
                 beta0 = -2.433, beta3 = 1.871)
)[[lag]]

# scenarios: true effect sizes for beta1 (internal) and beta2 (external)
scenarios <- expand.grid(beta1 = c(0, 0.5, 1, 2),
                         beta2 = c(0, 0.5, 1))

n_rep    <- 20      # replicates per scenario
n_iter   <- 6000    # MCMC iterations per fit (one chain per replicate)
n_burnin <- 2000
n_cores  <- min(8, detectCores())

out_file <- sprintf("data/sim_power/sim_power_%s.rds", lag)
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)

############################################
# data prep + model code from prod script  #
############################################

# run everything in the prod script up to (not including) the MCMC section,
# which defines model_code, constants, data, inits and zobs
src <- readLines(prod_file)
end <- grep("^# run MCMC in parallel", src) - 2
eval(parse(text = src[1:end]))

params <- c("psi", "p_minnow", "p_fukui", "p_shrimp",
            "beta0", "beta1", "beta2", "beta3")

# one task per scenario x replicate
tasks <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(s)
  data.frame(scenario = s, rep = seq_len(n_rep),
             beta1 = scenarios$beta1[s], beta2 = scenarios$beta2[s])))
tasks$seed <- 1000 + seq_len(nrow(tasks))

#############################################
# worker setup: build + compile model once  #
#############################################

cl <- makeCluster(n_cores)
clusterExport(cl, c("model_code", "constants", "data", "inits", "zobs",
                    "params", "truth_fixed", "n_iter", "n_burnin", "agg", "agg_w"))

invisible(clusterEvalQ(cl, {
  library(nimble)
  set.seed(1)
  m <- nimbleModel(model_code, constants = constants, data = data, inits = inits())
  conf <- configureMCMC(m, monitors = params)
  cm <- compileNimble(m)
  cmcmc <- compileNimble(buildMCMC(conf), project = m)

  # everything downstream of the parameters: gamma, epsilon, conn_*, z, p, y
  sim_nodes <- m$getDependencies(params, self = FALSE, downstream = TRUE)

  # site-years with any survey effort (main traps or WSG)
  surveyed <- matrix(FALSE, constants$nSites, constants$nYears)
  surveyed[cbind(agg$site, agg$year)] <- TRUE
  surveyed[cbind(agg_w$site, agg_w$year)] <- TRUE
  NULL
}))

run_rep <- function(task) {
  set.seed(task$seed)

  # 1. set true values and simulate latent states + detections
  truth <- c(truth_fixed, list(beta1 = task$beta1, beta2 = task$beta2))
  for (nm in names(truth)) cm[[nm]] <- truth[[nm]]
  cm$simulate(sim_nodes, includeData = TRUE)
  cm$calculate()

  z_true <- cm$z
  nT <- ncol(z_true)
  prev <- z_true[, -nT]; now <- z_true[, -1]
  surv_now <- surveyed[, -1]
  conn <- cm$conn_occ

  # observed (naive) occupancy from simulated detections
  obs <- matrix(NA, nrow(z_true), nT)
  y_det <- tapply(cm$y_agg > 0, list(factor(agg$site, 1:nrow(z_true)),
                                     factor(agg$year, 1:nT)), any)
  y_w   <- tapply(cm$y_wsg > 0, list(factor(agg_w$site, 1:nrow(z_true)),
                                     factor(agg_w$year, 1:nT)), any)
  obs[!is.na(y_det)] <- y_det[!is.na(y_det)]
  obs[!is.na(y_w)] <- pmax(obs[!is.na(y_w)], y_w[!is.na(y_w)], na.rm = TRUE)

  stats <- data.frame(
    occ_true        = mean(z_true),
    col_true        = sum(prev == 0 & now == 1),                    # all site-years
    col_true_surv   = sum(prev == 0 & now == 1 & surv_now),         # at surveyed site-years
    col_obs         = sum(obs[, -nT] %in% 0 & obs[, -1] %in% 1),    # naive 0 -> 1 detected
    prop_conn_pos   = mean(conn[prev == 0] > 0)                     # empty zones with occupied sources
  )

  # 2. refit from fresh inits
  ini <- inits()
  for (nm in params) cm[[nm]] <- ini[[nm]]
  cm$z <- zobs
  cm$calculate()
  cmcmc$run(n_iter, nburnin = n_burnin, reset = TRUE, progressBar = FALSE)
  s <- as.matrix(cmcmc$mvSamples)

  # 3. summarise
  res <- do.call(rbind, lapply(params, function(p) data.frame(
    param = p, true = truth[[p]],
    median = median(s[, p]),
    lo = unname(quantile(s[, p], 0.025)),
    hi = unname(quantile(s[, p], 0.975)))))
  cbind(task[rep(1, nrow(res)), c("scenario", "rep")], res, stats[rep(1, nrow(res)), ],
        row.names = NULL)
}

clusterExport(cl, "run_rep")
t0 <- Sys.time()
results <- do.call(rbind, parLapplyLB(cl, split(tasks, seq_len(nrow(tasks))), run_rep))
stopCluster(cl)
cat("elapsed:", format(Sys.time() - t0), "\n")

results <- merge(results, cbind(scenario = seq_len(nrow(scenarios)),
                                setNames(scenarios, c("true_beta1", "true_beta2"))))
saveRDS(results, out_file)

#############
# summary   #
#############

# power = CI excludes 0 (in the direction of the true effect, when nonzero)
# coverage = CI contains the true value
summ <- do.call(rbind, lapply(split(results, list(results$scenario, results$param)),
  function(d) {
    tv <- d$true[1]
    data.frame(true_beta1 = d$true_beta1[1], true_beta2 = d$true_beta2[1],
               param = d$param[1], true = tv,
               median_est = median(d$median),
               bias = mean(d$median - tv),
               ci_width = mean(d$hi - d$lo),
               coverage = mean(d$lo <= tv & d$hi >= tv),
               excl_zero = mean(if (tv > 0) d$lo > 0 else if (tv < 0) d$hi < 0
                                else d$lo > 0 | d$hi < 0),
               col_true_surv = mean(d$col_true_surv),
               col_obs = mean(d$col_obs),
               prop_conn_pos = mean(d$prop_conn_pos))
  }))
summ <- summ[summ$param %in% c("beta1", "beta2"), ]
summ <- summ[order(summ$param, summ$true_beta1, summ$true_beta2), ]
rownames(summ) <- NULL
print(summ, digits = 2)
write.csv(summ, sub("\\.rds$", "_summary.csv", out_file), row.names = FALSE)
