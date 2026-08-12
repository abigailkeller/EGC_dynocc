library(nimble)
library(MCMCvis)
library(parallel)

################
# read in data #
################

cpue <- readRDS("data/model_data/cpue_fukui.rds")[, 1:6]
detections <- readRDS("data/model_data/PresenceArrayBinary.rds")[, 2:7, ]
type <- readRDS("data/model_data/TrapTypeArraysNew.rds")
zones <- read.csv("data/model_data/site_zone_map.csv")[, "zone_id"]

# NA for Campbell Slough
zones[29] <- 1

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



# split up trap types
type_M <- type$Minnow
type_F <- type$Fukui
type_S <- type$Shrimp

# replace NA with 0
detections[is.na(detections)] <- 0
type_M[is.na(type_M)] <- 0
type_F[is.na(type_F)] <- 0
type_S[is.na(type_S)] <- 0

# read in connectivity data
larv_S <- array(data = NA, dim = c(nzone, nzone, nyear))
conn_paths <- c("data/connectivity/_zones_yearly_connectivity_matrix_counts_2018.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2019.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2020.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2021.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2022.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2023.csv",
                "data/connectivity/_zones_yearly_connectivity_matrix_counts_2024.csv")
for (i in 1:nyear) {
  larv_S[1:nzone, 1:nzone, i] <- as.matrix(read.csv(conn_paths[i])[, 2:15])
}

# make arbitrary number of released
larv_R <- matrix(400, nrow = nzone, ncol = nyear)


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
  beta1 ~ dunif(-10, 10)
  beta2 ~ dunif(-10, 10)
  
  # Persistence intercept
  beta3 ~ dunif(-10, 10)
  
  # Persistence coefficient
  beta4 ~ dunif(-10, 10)
  
  # --- Initial occupancy (t = 1) ---
  for (i in 1:nSites) {
    z[i, 1] ~ dbern(psi)
  }
  
  # --- Connectivity (random variable)
  for (i in 1:nZones) {
    for (t in 1:nYears) {
      for (k in 1:nZones) {
        larv_S[i, k, t] ~ dbinom(size = larv_R[i, t], prob = C[i, k, t])
      }
    }
  }
  
  # --- expand latent zone-level C to site level ---
  for (i in 1:nSites) {
    for (t in 1:nYears) {
      C_self[i, t] <- C[zones[i], zones[i], t] # within-zone (persistence)
      for (k in 1:nSites) {
        C_site[i, k, t] <- C[zones[i], zones[k], t] # site-pair (colonization)
      }
    }
  }
  
  # --- Temporal dynamics (t = 2, ..., T) ---
  for (i in 1:nSites) {
    
    for (t in 2:nYears) {
      
      # Probability of colonization
      logit(gamma[i, t - 1]) <- colonization(CPUE[1:nSites, (t-1):t], 
                                             beta0, beta1, beta2,
                                             C_site[i, 1:nSites, (t-1):t], i)
      
      # Probability of persistence
      logit(epsilon[i, t - 1]) <- persistence(CPUE[i, t], beta3, beta4, 
                                              C_self[i, t])
        
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
  for (i in 1:nSites) {
    for (t in 1:nYears) {
      if (nTraps[i, t] > 0) {
        for (j in 1:nTraps[i, t]) {
          p[i, t, j] <- p_minnow * type_M[i, t, j] + p_fukui * type_F[i, t, j] +
            p_shrimp * type_S[i, t, j]
          y[i, t, j] ~ dbern(z[i, t] * p[i, t, j])
        }
      }
    }
  }
  
})

# Package data and constants
constants <- list(
  nSites = nsite,
  nYears = nyear,
  nTraps = ntraps,
  nZones = nzone,
  zones = as.integer(zones),
  type_M = type_M,
  type_F = type_F,
  type_S = type_S
)

data <- list(y = detections, # dimensions [sites, years, traps]
             CPUE = cpue, # dimensions [sites, years + 1 year burnin]
             larv_R = larv_R, # no. released larvae [zones, years + 1 year burnin]
             larv_S = larv_S # no. settled larvae [zones_set, zones_rel, years + 1 year burnin]
             ) 

# connectivity inits
C <- array(NA, dim = c(nzone, nzone, nyear))
for (i in 1:nzone) {
  for (k in 1:nzone) {
    for (t in 1:nyear) {
      C[i, k, t] <- larv_S[i, k, t] / larv_R[i, t]
    }
  }
}

# Initial values
inits  <- function(detections) {
  list(psi = runif(1, 0, 1),
       beta0 = runif(1, -1, 1),
       beta1 = runif(1, -1, 1),
       beta2 = runif(1, -1, 1),
       beta3 = runif(1, -1, 1),
       beta4 = runif(1, -1, 1),
       p_minnow = runif(1, 0, 1),
       p_fukui = runif(1, 0, 1),
       p_shrimp = runif(1, 0, 1),
       z = apply(detections, c(1, 2), max, na.rm = TRUE),
       C = C)
}

########################
# run MCMC in parallel #
########################

cl <- makeCluster(4)

set.seed(10120)

clusterExport(cl, c("model_code", "inits", "data", "constants", 
                    "detections", "C"))

# parallelize running MCMC
out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  # Define nimbleFunctions directly on each worker
  colonization <- nimbleFunction(
    run = function(CPUE = double(2), beta0 = double(0), beta1 = double(0), 
                   beta2 = double(0), C = double(2), index = double(0))
    {
      returnType(double(0))
      idx <- as.integer(index)
      n <- dim(CPUE)[1]
      gamma <- 0
      for (j in 1:n) {
        if (j != idx) {
          gamma <- gamma + CPUE[j, 2] * beta1 * C[j, 2] +
            CPUE[j, 1] * beta2 * C[j, 1]
        }
      }
      gamma <- gamma + beta0
      return(gamma)
    }
  )
  assign("colonization", colonization, envir = .GlobalEnv)
  
  persistence <- nimbleFunction(
    run = function(CPUE = double(0), beta3 = double(0), beta4 = double(0), 
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
                         inits = inits(detections))
  
  
  # build the MCMC
  mcmcConf_myModel <- configureMCMC(
    myModel,
    monitors = c("psi", "beta0", "beta1", "beta2", "beta3", "beta4", "p"),
    enableWAIC = TRUE
    )
  
  # build MCMC
  myMCMC <- buildMCMC(mcmcConf_myModel)
  
  # compile the model and MCMC
  CmyModel <- compileNimble(myModel)
  
  # compile the MCMC
  cmodel_mcmc <- compileNimble(myMCMC, project = myModel)
  
  # run MCMC
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
                out[[3]][sequence, ], out[[4]][sequence, ])

# save samples
saveRDS(out_sub, "posterior_samples_sim.rds")

stopCluster(cl)


#####################
# summarize samples #
#####################

summmary <- MCMCsummary(out_sub)
