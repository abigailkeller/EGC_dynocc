library(nimble)
library(MCMCvis)
library(parallel)
library(tidyverse)

################
# read in data #
################

# read in connectivity data - larvae settled
conn_paths <- c("data/connectivity/2013_connectivity_counts_365.csv",
                "data/connectivity/2014_connectivity_counts_365.csv",
                "data/connectivity/2015_connectivity_counts_365.csv",
                "data/connectivity/2016_connectivity_counts_365.csv",
                "data/connectivity/2017_connectivity_counts_365.csv",
                "data/connectivity/2018_connectivity_counts_365.csv",
                "data/connectivity/2019_connectivity_counts_365.csv",
                "data/connectivity/2020_connectivity_counts_365.csv",
                "data/connectivity/2021_connectivity_counts_365.csv",
                "data/connectivity/2022_connectivity_counts_365.csv")
conn_2013 <- read.csv(conn_paths[1])
site_names <- conn_2013$settlement_site_id
nsite <- length(site_names)
nyear <- length(conn_paths)
larv_S <- array(data = NA, dim = c(nsite, nsite, nyear))
larv_S[1:nsite, 1:nsite, 1] <- as.matrix(conn_2013[, 3:(nsite + 2)])
for (i in 2:length(conn_paths)) {
  larv_S[1:nsite, 1:nsite, i] <- as.matrix(
    read.csv(conn_paths[i])[, 3:(nsite + 2)])
}

# read in connectivity data - larvae released
larv_R <- as.matrix(read.csv("data/connectivity/larvae_released.csv",
                             row.names = 1)[, 2:11])

# read in connectivity data - external settlers
larv_E <- as.matrix(read.csv("data/connectivity/larvae_external_settlers.csv",
                             row.names = 1)[, 2:11])

# rescale connectivity data
# larv_S[i, j, t]: settlers into zone i (row) released from zone j (column) in year t
# (larv_R rows/cols already align with site_names and conn years 2013-2022)

# per-source dispersal probability: P[i, j, t] = share of larvae released from j that settle in i
P <- larv_S
for (t in 1:nyear) {
  R <- larv_R[, t]
  P[, , t] <- sweep(larv_S[, , t], 2, ifelse(R > 0, R, 1), "/")
  P[, R == 0, t] <- 0                 # sources that released nothing contribute nothing
}

# log-transform (small offset for zeros) and z-scale
zscale_prop <- function(m) { l <- log(m + 0.001); (l - mean(l)) / sd(l) }
zscale <- function(m) { l <- log1p(m); (l - mean(l)) / sd(l) }

# internal immigration: summed dispersal probability from all other modeled zones
imm_int <- sapply(1:nyear, function(t) rowSums(P[, , t]) - diag(P[, , t]))
z_imm_int <- zscale_prop(imm_int)

# external settlers: no release total to divide by, so keep as counts
z_imm_ext <- zscale(larv_E)

# retention: share of larvae released from zone i that settle back in i
ret_prop <- sapply(1:nyear, function(t) diag(P[, , t]))
z_ret <- zscale_prop(ret_prop)

# 2013 - 2024

#cpue <- readRDS("data/model_data/cpue_fukui_expanded.rds")[, 1:10]
detections <- readRDS("data/model_data/PresenceArrayBinary_expanded.rds")[, 2:11, ]
detections_WSG <- readRDS("data/model_data/PresenceArrayBinary_WSG_expanded.rds")[, 2:11, ]
type <- readRDS("data/model_data/TrapTypeArraysNew_expanded.rds")
#wsg_map <- read.csv("data/model_data/wsg_map.csv")

# split up trap types
type_M <- type$Minnow[, 2:11, ]
type_F <- type$Fukui[, 2:11, ]
type_S <- type$Shrimp[, 2:11, ]

# remove sites not in larval data
site_names_main <- rownames(detections[, 1, ])
keep <- which(site_names_main %in% site_names)
detections <- detections[keep, , ]
type_M <- type_M[keep, , ]; type_F <- type_F[keep, , ]; type_S <- type_S[keep, , ]
site_names_main <- site_names_main[keep]
nsite_main <- length(site_names_main)

# get ntraps for each site/year for non-WSG detections
ntraps <- matrix(NA, nrow = nsite_main, ncol = nyear)
for (i in 1:nsite_main) {
  for (t in 1:nyear) {
    ntraps[i, t] <- sum(!is.na(detections[i, t,]))
  }
}

##
# flatten observation data to vectors
##

# map main/WSG site rows to the larval-site index
larv_ind <- match(site_names_main, site_names)

# trap data: one row per site x year x trap type
has_obs <- !is.na(detections)
types   <- list(type_M, type_F, type_S)   # 1 = minnow, 2 = fukui, 3 = shrimp

agg <- do.call(rbind, lapply(1:3, function(k) {
  is_k <- has_obs & !is.na(types[[k]]) & types[[k]] == 1
  n <- rowSums(is_k, dims = 2)                      # traps of type k
  y <- rowSums(is_k & detections == 1, dims = 2)    # detections in those traps
  idx <- which(n > 0, arr.ind = TRUE)
  data.frame(site = larv_ind[idx[, 1]], year = idx[, 2], type = k,
             n = n[idx], y = y[idx])
}))

# WSG data: one row per site x year
n_w <- rowSums(!is.na(detections_WSG), dims = 2)
y_w <- rowSums(detections_WSG, dims = 2, na.rm = TRUE)
idx <- which(n_w > 0, arr.ind = TRUE)
agg_w <- data.frame(site = match(rownames(detections_WSG), site_names)[idx[, 1]],
                    year = idx[, 2], n = n_w[idx], y = y_w[idx])



##############
# model code #
##############

model_code <- nimbleCode({
  
  # --- Priors ---
  
  # Probability of detection
  p_minnow ~ dunif(0, 1)
  p_fukui ~ dunif(0, 1)
  p_shrimp ~ dunif(0, 1)
  
  # Initial occupancy probability
  psi ~ dunif(0, 1)
  
  beta0 ~ dunif(-5, 5)
  beta1 ~ dunif(-5, 5) 
  beta2 ~ dunif(-5, 5) 
  beta3 ~ dunif(-5, 5)  
  beta4 ~ dunif(-5, 5)
  
  # --- Initial occupancy (t = 1) ---
  for (i in 1:nSites) {
    z[i, 1] ~ dbern(psi)
  }
  
  # connectivity
  for (i in 1:nSites) {
    for (t in 2:nYears) {
      # colonization
      logit(gamma[i, t - 1]) <- beta0 + beta1 * z_imm_int[i, t] +
        beta4 * z_imm_ext[i, t]
      # persistence
      logit(epsilon[i, t - 1]) <- beta2 + beta3 * z_ret[i, t]
    }
  }
  
  
  
  # --- Temporal dynamics (t = 2, ..., T) ---
  for (i in 1:nSites) {
    
    for (t in 2:nYears) {
        
        # Occupancy state
        z[i, t] ~ dbern(
          # probability occupied and persisted
          z[i, t - 1] * epsilon[i, t - 1] +
            # probability unoccupied and colonized
            (1 - z[i, t - 1]) * gamma[i, t - 1]
        )
      
    }
  }
  
  # --- Observation model ---
  p[1] <- p_minnow
  p[2] <- p_fukui
  p[3] <- p_shrimp
  
  # main data
  for (k in 1:nAgg) {
    y_agg[k] ~ dbin(z[agg_site[k], agg_year[k]] * p[agg_type[k]], n_agg[k])
  }
  
  # wsg data
  p_star <- 1 - (1 - p_minnow) ^ 3 * (1 - p_fukui) ^ 3
  for (k in 1:nAgg_w) {
    y_wsg[k] ~ dbin(z[wsg_site[k], wsg_year[k]] * p_star, n_wsg[k])
  }
  
})

# Package data and constants
constants <- list(nSites = nsite, 
                  nYears = nyear,
                  nAgg = nrow(agg), 
                  agg_site = agg$site, 
                  agg_year = agg$year,
                  agg_type = agg$type, 
                  n_agg = agg$n,
                  nAgg_w = nrow(agg_w), 
                  wsg_site = agg_w$site,
                  wsg_year = agg_w$year, 
                  n_wsg = agg_w$n)

data <- list(y_agg = agg$y, 
             y_wsg = agg_w$y, 
             z_imm_int = z_imm_int, 
             z_imm_ext = z_imm_ext, 
             z_ret = z_ret)

# occupancy inits
zobs <- matrix(1, nrow = nsite, ncol = nyear)

# Initial values
inits  <- function() {
  list(psi = runif(1, 0, 1),
       beta0 = runif(1, -1, 1),
       beta1 = runif(1, 0, 1),
       beta2 = runif(1, 0, 1),
       beta3 = runif(1, -1, 1),
       beta4 = runif(1, -1, 1),
       p_minnow = runif(1, 0, 1),
       p_fukui = runif(1, 0, 1),
       p_shrimp = runif(1, 0, 1),
       z = zobs)
}

########################
# run MCMC in parallel #
########################

cl <- makeCluster(8)

set.seed(10120)

clusterExport(cl, c("model_code", "inits", "data", "constants", "zobs"))

# parallelize running MCMC
out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  
  # build model
  myModel <- nimbleModel(code = model_code,
                         data = data,
                         constants = constants,
                         inits = inits())
  
  
  # build the MCMC
  mcmcConf_myModel <- configureMCMC(
    myModel,
    monitors = c("psi", "beta0", "beta1", "beta2", "beta3", "beta4",
                 "p_minnow", "p_fukui", "p_shrimp", "z",
                 "gamma", "epsilon"),
    enableWAIC = TRUE
    )
  
  # build MCMC
  myMCMC <- buildMCMC(mcmcConf_myModel)
  
  # compile the model and MCMC
  CmyModel <- compileNimble(myModel)
  
  # compile the MCMC
  cmodel_mcmc <- compileNimble(myMCMC, project = myModel)
  
  # run MCMC
  # cmodel_mcmc$run(1000000, thin = 1000,
  #                 reset = FALSE)
  cmodel_mcmc$run(10000, thin = 10,
                  reset = FALSE)
  
  samples <- as.mcmc(as.matrix(cmodel_mcmc$mvSamples))
  
  return(samples)
})

# discard burnin
lower <- 200
upper <- dim(out[[1]])[1]
sequence <- seq(lower, upper, 1)
out_sub <- list(out[[1]][sequence, ], out[[2]][sequence, ],
                out[[3]][sequence, ], out[[4]][sequence, ],
                out[[5]][sequence, ], out[[6]][sequence, ],
                out[[7]][sequence, ], out[[8]][sequence, ])

# save samples
saveRDS(out_sub, "data/posterior_samples/posterior_samples_20261111_rescale_2year.rds")

stopCluster(cl)

