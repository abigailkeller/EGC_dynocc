library(nimble)
library(MCMCvis)
library(parallel)
library(tidyverse)

################
# read in data #
################

cpue <- readRDS("data/model_data/cpue_zone_year.rds")[c(1, 3:14), 13:18]
detections <- readRDS("data/model_data/PresenceArrayBinary.rds")[, 2:7, ]
type <- readRDS("data/model_data/TrapTypeArraysNew.rds")
zones <- read.csv("data/model_data/site_zone_map.csv")[, "zone_id"]

##############
# clean data #
##############

##
# get sites without traps
##
##

remove <- as.integer(which(apply(detections, 1, function(s) all(is.na(s)))))

# remove sites without traps
detections <- detections[-remove, , ]
zones <- zones[-remove]

# split up trap types
type_M <- type$Minnow[-remove, 2:7, ]
type_F <- type$Fukui[-remove, 2:7, ]
type_S <- type$Shrimp[-remove, 2:7, ]

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


# get constants
nyear <- dim(detections)[2]
nsite <- dim(detections)[1]
nzone <- length(unique(zones))

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

# read in connectivity data - larvae settled
larv_S <- array(data = NA, dim = c(nzone, nzone, nyear))
conn_paths <- c("data/connectivity/_zones_yearly_connectivity_matrix_counts_2018.csv",
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
                             row.names = 1)[c(1, 3:14), 6:11])


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
        beta0, beta1, beta2,
        C[i, 1:nZones, (t-1):t], i, Kscale)
      
      # Probability of persistence
      logit(epsilon[i, t - 1]) <- persistence(CPUE[i, t], 
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
  for (o in 1:nObs) {
    p_long[o] <- p_minnow * obs_yM[o] + p_fukui * obs_yF[o] + 
      p_shrimp * obs_yS[o]
    y_long[o] ~ dbern(z[obs_site[o], obs_year[o]] * p_long[o])
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
zobs <- apply(detections, c(1, 2), function(x) {
  if (all(is.na(x))) 0 else max(x, na.rm = TRUE)
})
dimnames(zobs) <- NULL

# Package data and constants
constants <- list(
  nSites = nsite,
  nYears = nyear,
  nZones = nzone,
  nObs = nObs,
  zones = as.integer(zones),
  obs_yM = obs_yM,
  obs_yF = obs_yF,
  obs_yS = obs_yS,
  obs_site = obs_site,
  obs_year = obs_year,
  Kscale = Kscale
)

data <- list(y_long = y_long, # dimensions [sites, years, traps]
             CPUE = cpue, # dimensions [zones, years]
             larv_R = larv_R, # no. released larvae [zones, years]
             larv_S = larv_S # no. settled larvae [zones_set, zones_rel, years]
) 

# Initial values
inits  <- function() {
  list(psi = runif(1, 0, 1),
       beta0 = runif(1, -1, 1),
       beta1 = runif(1, -1, 1),
       beta2 = runif(1, -1, 1),
       beta3 = runif(1, -1, 1),
       beta4 = runif(1, -1, 1),
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
  cmodel_mcmc$run(100000, thin = 100,
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
saveRDS(out_sub, "data/posterior_samples/posterior_samples_20260902.rds")

stopCluster(cl)

