library(nimble)
library(MCMCvis)
library(parallel)
library(tidyverse)

################
# read in data #
################

# 2014 - 2024

cpue <- readRDS("data/model_data/cpue_zone_year.rds")[c(1, 3:14), 8:18]
detections <- readRDS("data/model_data/PresenceArrayBinary.rds")#[, 2:7, ]
detections_WSG <- readRDS("data/model_data/PresenceArrayBinary_WSG.rds")#[, 2:7, ]
type <- readRDS("data/model_data/TrapTypeArraysNew.rds")
zones <- read.csv("data/model_data/site_zone_map.csv")[, "zone_id"]
wsg_zones <- read.csv("data/model_data/site_zone_map_WSG.csv")[, "zone_id"]
wsg_map <- read.csv("data/model_data/wsg_map.csv")

##############
# clean data #
##############

# replace cpue NA with 0
cpue[is.na(cpue)] <- 0

##
# get sites without traps
##
##

# remove <- as.integer(which(apply(detections, 1, function(s) all(is.na(s)))))
# remove_WSG <- as.integer(which(apply(detections_WSG, 1, 
#                                      function(s) all(is.na(s)))))
# 
# # remove sites without traps
# detections <- detections[-remove, , ]
# zones <- zones[-remove]
# detections_WSG <- detections_WSG[-remove_WSG, , ]
# wsg_map <- wsg_map[-remove_WSG, , ]
# wsg_zones <- wsg_zones[-remove_WSG]

# split up trap types
type_M <- type$Minnow#[-remove, 2:7, ]
type_F <- type$Fukui#[-remove, 2:7, ]
type_S <- type$Shrimp#[-remove, 2:7, ]

# NA for Campbell Slough
zones[which(is.na(zones))] <- 1

##
# remove zone 2
##
##

zone2 <- which(zones == 2)

detections <- detections[-zone2, , ]
zones <- zones[-zone2]

# trap types
type_M <- type_M[-zone2, , ]
type_F <- type_F[-zone2, , ]
type_S <- type_S[-zone2, , ]

# replace zone 14 with zone 2
zones[which(zones == 14)] <- 2
wsg_zones[which(wsg_zones == 14)] <- 2

# get constants
nyear <- dim(detections)[2]
nsite <- dim(detections)[1]
nzone <- length(unique(zones))
nyear_wsg <- dim(detections_WSG)[2]
nsite_wsg <- dim(detections_WSG)[1]

# create index for WSG sites
site_names <- rownames(detections[, 1, ])
site_names_wsg <- rownames(detections_WSG)
ind_WSG <- rep(NA, length(site_names_wsg))
site_counter <- nsite + 1
for (i in 1:length(ind_WSG)) {
  name <- wsg_map[which(wsg_map[, "WSG"] == site_names_wsg[i]), "Other"]
  if (length(name) == 1) {
    ind_WSG[i] <- which(site_names == name)
  } else {
    ind_WSG[i] <- site_counter
    site_counter <- site_counter + 1
  }
}
nsite_total <- site_counter - 1

# get ntraps for each site/year
ntraps <- matrix(NA, nrow = nsite, ncol = nyear)
for (i in 1:nsite) {
  for (t in 1:nyear) {
    ntraps[i, t] <- sum(!is.na(detections[i, t,]))
  }
}

# flatten observation data to vectors
nObs <- sum(ntraps)
obs_site <- rep(0, nObs)
obs_year <- rep(0, nObs)
obs_yM <- rep(0, nObs)
obs_yF <- rep(0, nObs)
obs_yS <- rep(0, nObs)
y_long <- rep(0, nObs)

ind <- 1
for (i in 1:nsite) {
  for (t in 1:nyear) {
    if (ntraps[i, t] > 0) {
      obs_site[ind:(ind + ntraps[i, t] - 1)] <- i
      obs_year[ind:(ind + ntraps[i, t] - 1)] <- t
      obs_yM[ind:(ind + ntraps[i, t] - 1)] <- type_M[i, t, 1:ntraps[i, t]]
      obs_yF[ind:(ind + ntraps[i, t] - 1)] <- type_F[i, t, 1:ntraps[i, t]]
      obs_yS[ind:(ind + ntraps[i, t] - 1)] <- type_S[i, t, 1:ntraps[i, t]]
      y_long[ind:(ind + ntraps[i, t] - 1)] <- detections[i, t, 1:ntraps[i, t]]
      
      ind <- ind + ntraps[i, t]
    }
  }
}

# replace NA with 0
obs_yM[is.na(obs_yM)] <- 0
obs_yF[is.na(obs_yF)] <- 0
obs_yS[is.na(obs_yS)] <- 0

# get nvisits for each site/year in WSG data
nvisits_wsg <- matrix(NA, nrow = nsite_wsg, ncol = nyear_wsg)
for (i in 1:nsite_wsg) {
  for (t in 1:nyear_wsg) {
    nvisits_wsg[i, t] <- sum(!is.na(detections_WSG[i, t,]))
  }
}

# flatten WSG observation data to vectors
nObs_wsg <- sum(nvisits_wsg)
obs_site_wsg <- rep(0, nObs_wsg)
obs_year_wsg <- rep(0, nObs_wsg)
y_long_WSG <- rep(0, nObs_wsg)

ind <- 1
for (i in 1:nsite_wsg) {
  for (t in 1:nyear_wsg) {
    keep <- which(!is.na(detections_WSG[i, t, ]))
    if (length(keep) > 0) {
      n <- length(keep)
      obs_site_wsg[ind:(ind + n - 1)] <- ind_WSG[i]
      obs_year_wsg[ind:(ind + n - 1)] <- t
      y_long_WSG[ind:(ind + n - 1)] <- detections_WSG[i, t, keep]
      ind <- ind + n
    }
  }
}

# get all zones
zones_full <- rep(NA, nsite_total)
zones_full[1:nsite] <- as.integer(zones)
for (i in seq_along(ind_WSG)) {
  if (ind_WSG[i] > nsite) zones_full[ind_WSG[i]] <- as.integer(wsg_zones[i])
}

# read in connectivity data - larvae settled
larv_S <- array(data = NA, dim = c(nzone, nzone, nyear))
conn_paths <- c("data/connectivity/_zones_yearly_connectivity_matrix_counts_2013.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2014.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2015.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2016.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2017.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2018.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2019.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2020.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2021.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2022.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2023.csv")
for (i in 1:nyear) {
  larv_S[1:nzone, 1:nzone, i] <- as.matrix(
    read.csv(conn_paths[i])[, 2:15])[c(1, 3:14), c(1, 3:14)]
}

# read in connectivity data - larvae settled
larv_R <- as.matrix(read.csv("data/SpatialData/yearly_larvae_released.csv",
                             row.names = 1)[c(1, 3:14), 1:11])


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
  
  # Colonization intercept
  beta0 ~ dunif(-10, 10)
  
  # Colonization coefficients
  beta1 ~ dunif(0, 1000)
  
  # Persistence intercept
  beta2 ~ dunif(-10, 10)
  
  # Persistence coefficient
  beta3 ~ dunif(0, 1000)
  
  # --- Initial occupancy (t = 1) ---
  for (i in 1:nSites) {
    z[i, 1] ~ dbern(psi)
  }
  
  # --- Connectivity (random variable)
  for (i in 1:nZones) {
    for (t in 1:nYears) {
      for (k in 1:nZones) {
        
        # probability of connectivity
        C[i, k, t] ~ dbeta(1, 1)
        larv_S[i, k, t] ~ dbinom(prob = C[i, k, t], size = larv_R[i, t])
        
      }
    }
    
    for (t in 2:nYears) {
      
      # Probability of colonization
      logit(gamma[i, t - 1]) <- colonization(CPUE[1:nZones, (t-1):t], 
                                             beta0, beta1,
                                             C[i, 1:nZones, (t-1):t], i)
      
      # Probability of persistence
      logit(epsilon[i, t - 1]) <- persistence(CPUE[i, t], 
                                              beta2, beta3, 
                                              C[i, i, t])
    }
  }
  
  
  
  # --- Temporal dynamics (t = 2, ..., T) ---
  for (i in 1:nSites) {
    
    for (t in 2:nYears) {
      
      # Occupancy state
      z[i, t] ~ dbern(
        # probability occupied and persisted
        z[i, t - 1] * epsilon[zones[i], t - 1] +
          # probability unoccupied and colonized
          (1 - z[i, t - 1]) * gamma[zones[i], t - 1]
      )
      
    }
  }
  
  # --- Observation model ---
  # trap-level data
  for (o in 1:nObs) {
    p_long[o] <- p_minnow * obs_yM[o] + p_fukui * obs_yF[o] + 
      p_shrimp * obs_yS[o]
    y_long[o] ~ dbern(z[obs_site[o], obs_year[o]] * p_long[o])
  }
  # aggregated WSG data
  p_star <- 1 - (1 - p_minnow) ^ 3 * (1 - p_fukui) ^ 3
  for (o in 1:nObs_wsg) {
    y_long_WSG[o] ~ dbern(z[obs_site_WSG[o], obs_year_WSG[o]] * p_star)
  }
  
})


# connectivity inits
C_hat <- array(NA, dim = c(nzone, nzone, nyear))
for (i in 1:nzone) {
  for (k in 1:nzone) {
    for (t in 1:nyear) {
      C_hat[i, k, t] <- larv_S[i, k, t] / larv_R[i, t]
    }
  }
}

# occupancy inits
# zobs <- apply(detections, c(1, 2), function(x) {
#   if (all(is.na(x))) 0 else max(x, na.rm = TRUE)
# })
# dimnames(zobs) <- NULL
zobs <- matrix(1, nrow = nsite_total, ncol = nyear)

# Package data and constants
constants <- list(
  nSites = nsite_total,
  nYears = nyear,
  nZones = nzone,
  nObs = nObs,
  zones = as.integer(zones_full),
  obs_yM = obs_yM,
  obs_yF = obs_yF,
  obs_yS = obs_yS,
  obs_site = obs_site,
  obs_year = obs_year,
  nObs_wsg = nObs_wsg,
  obs_site_WSG = obs_site_wsg,
  obs_year_WSG = obs_year_wsg
)

data <- list(y_long = y_long, 
             y_long_WSG = y_long_WSG,
             CPUE = cpue, # dimensions [zones, years]
             larv_R = larv_R, # no. released larvae [zones, years]
             larv_S = larv_S # no. settled larvae [zones_set, zones_rel, years]
) 

# Initial values
inits  <- function() {
  list(psi = runif(1, 0, 1),
       beta0 = runif(1, -1, 1),
       beta1 = runif(1, 0, 1),
       beta2 = runif(1, -1, 1),
       beta3 = runif(1, 0, 1),
       p_minnow = runif(1, 0, 1),
       p_fukui = runif(1, 0, 1),
       p_shrimp = runif(1, 0, 1),
       z = zobs,
       C = C_hat)
}

########################
# run MCMC in parallel #
########################

cl <- makeCluster(8)

set.seed(10120)

clusterExport(cl, c("model_code", "inits", "data", "constants", 
                    "C_hat", "zobs"))

# parallelize running MCMC
out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  # Define nimbleFunctions directly on each worker
  colonization <- nimbleFunction(
    run = function(CPUE = double(2), 
                   beta0 = double(0), beta1 = double(0), 
                   C = double(2), index = double(0))
    {
      returnType(double(0))
      idx <- as.integer(index)
      n <- dim(C)[1]
      gamma <- 0
      for (j in 1:n) {
        if (j != idx) {
          gamma <- gamma + beta1 * (CPUE[j, 2] * C[j, 2] +
                                      CPUE[j, 1] * C[j, 1])
        }
      }
      gamma <- gamma + beta0
      return(gamma)
    }
  )
  assign("colonization", colonization, envir = .GlobalEnv)
  
  persistence <- nimbleFunction(
    run = function(CPUE = double(0), 
                   beta3 = double(0), beta4 = double(0), 
                   C = double(0))
    {
      returnType(double(0))
      epsilon <- beta3 + CPUE * C * beta4
      return(epsilon)
    }
  )
  assign("persistence", persistence, envir = .GlobalEnv)
  
  # build model
  myModel <- nimbleModel(code = model_code,
                         data = data,
                         constants = constants,
                         inits = inits())
  
  
  # build the MCMC
  mcmcConf_myModel <- configureMCMC(
    myModel,
    monitors = c("psi", "beta0", "beta1", "beta2", "beta3",
                 "p_minnow", "p_fukui", "p_shrimp", "C", "z",
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
saveRDS(out_sub, "data/posterior_samples/model_selection/model2.rds")

stopCluster(cl)


# calculate WAIC
samples_mat <- rbind(out[[1]][sequence, ], out[[2]][sequence, ],
                     out[[3]][sequence, ], out[[4]][sequence, ],
                     out[[5]][sequence, ], out[[6]][sequence, ],
                     out[[7]][sequence, ], out[[8]][sequence, ])
calculateWAIC(samples_mat, CmyModel)
# WAIC: 16532.86
# lppd: -7381.831
# pWAIC: 884.5978
