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
                "data/connectivity/2022_connectivity_counts_365.csv",
                "data/connectivity/2023_connectivity_counts_365.csv")
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

# read in connectivity data - larvae settled
larv_R <- as.matrix(read.csv("data/connectivity/larvae_released.csv",
                             row.names = 1)[, 2:12])

# read in connectivity data - external settlers
larv_E <- as.matrix(read.csv("data/connectivity/larvae_external_settlers.csv",
                             row.names = 1)[, 2:12])

# 2014 - 2024

#cpue <- readRDS("data/model_data/cpue_fukui_expanded.rds")[, 1:10]
detections <- readRDS("data/model_data/PresenceArrayBinary_expanded.rds")
detections_WSG <- readRDS("data/model_data/PresenceArrayBinary_WSG_expanded.rds")
type <- readRDS("data/model_data/TrapTypeArraysNew_expanded.rds")
#wsg_map <- read.csv("data/model_data/wsg_map.csv")

# split up trap types
type_M <- type$Minnow
type_F <- type$Fukui
type_S <- type$Shrimp

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

# flatten observation data to vectors
nObs <- sum(ntraps)
obs_site <- rep(0, nObs)
obs_year <- rep(0, nObs)
obs_yM <- rep(0, nObs)
obs_yF <- rep(0, nObs)
obs_yS <- rep(0, nObs)
y_long <- rep(0, nObs)

ind <- 1
for (i in 1:nsite_main) {
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
nsite_wsg <- dim(detections_WSG)[1]
nvisits_wsg <- matrix(NA, nrow = nsite_wsg, ncol = nyear)
for (i in 1:nsite_wsg) {
  for (t in 1:nyear) {
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
  for (t in 1:nyear) {
    keep <- which(!is.na(detections_WSG[i, t, ]))
    if (length(keep) > 0) {
      n <- length(keep)
      obs_site_wsg[ind:(ind + n - 1)] <- i
      obs_year_wsg[ind:(ind + n - 1)] <- t
      y_long_WSG[ind:(ind + n - 1)] <- detections_WSG[i, t, keep]
      ind <- ind + n
    }
  }
}



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
  beta1 ~ dunif(0, 10)
  beta2 ~ dunif(0, 10)
  
  # Persistence intercept
  beta3 ~ dunif(-10, 10)
  
  # Persistence coefficient
  beta4 ~ dunif(0, 10)
  
  # --- Initial occupancy (t = 1) ---
  for (i in 1:nSites) {
    z[i, 1] ~ dbern(psi)
  }
  
  # --- Connectivity (random variable)
  for (i in 1:nSites) {
    for (t in 1:nYears) {
      for (k in 1:nSites) {
        
        # probability of connectivity
        C[i, k, t] ~ dbeta(1, 1)
        larv_S[i, k, t] ~ dbinom(prob = C[i, k, t], size = larv_R[i, t])
        
      }
    }
    
    for (t in 2:nYears) {
      
      # Probability of colonization
      logit(gamma[i, t - 1]) <- colonization(
        beta0, beta1, beta2,
        C[i, 1:nZones, (t-1):t], i, Kscale)
      
      # Probability of persistence
      logit(epsilon[i, t - 1]) <- persistence(
        beta3, beta4, 
        C[i, i, t], Kscale)
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

# create scale for colonization calculation
X <- matrix(NA, nzone, nyear)
for (i in 1:nzone) for (t in 1:nyear) {
  j <- setdiff(1:nzone, i)
  X[i, t] <- sum(cpue[j, t] * C_hat[i, j, t])
}
Kscale <- sd(as.vector(X), na.rm = TRUE)

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
  Kscale = Kscale,
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
       beta2 = runif(1, 0, 1),
       beta3 = runif(1, -1, 1),
       beta4 = runif(1, 0, 1),
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
    run = function(#CPUE = double(2), 
                   beta0 = double(0), beta1 = double(0), 
                   beta2 = double(0), C = double(2), index = double(0),
                   Kscale = double(0))
    {
      returnType(double(0))
      idx <- as.integer(index)
      n <- dim(C)[1]
      gamma <- 0
      for (j in 1:n) {
        if (j != idx) {
          gamma <- gamma + CPUE[j, 2] * beta1 * C[j, 2] / Kscale +
            CPUE[j, 1] * beta2 * C[j, 1] / Kscale
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
                   C = double(0), Kscale = double(0))
    {
      returnType(double(0))
      epsilon <- beta3 + CPUE * C * beta4 / Kscale
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
    monitors = c("psi", "beta0", "beta1", "beta2", "beta3", "beta4", 
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
saveRDS(out_sub, "data/posterior_samples/posterior_samples_20260917_short.rds")

stopCluster(cl)

