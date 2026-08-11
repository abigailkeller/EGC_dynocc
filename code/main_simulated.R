library(nimble)
library(MCMCvis)
library(parallel)

#################
# simulate data #
#################

source("simulate_data.R")


##############
# model code #
##############

model_code <- nimbleCode({
  
  # --- Priors ---
  
  # Probability of detection
  p ~ dunif(0, 1) 
  
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
  for (i in 1:nSites) {
    for (t in 1:nYears) {
      for (k in 1:nSites) {
        larv_S[i, k, t] ~ dbinom(size = larv_R[i, t], prob = C[i, k, t])
      }
    }
  }
  
  # --- Temporal dynamics (t = 2, ..., T) ---
  for (i in 1:nSites) {
    
    for (t in 2:nYears) {
      
      # Probability of colonization
      logit(gamma[i, t - 1]) <- colonization_lp(CPUE[1:nSites, (t-1):t], 
                                                beta0, beta1, beta2, 
                                                C[i, 1:nSites, (t-1):t], 
                                                i)
      
      # Probability of persistence
      logit(epsilon[i, t - 1]) <- persistence_lp(CPUE[i, t], beta3, beta4,
                                                 C[i, i, t])
        
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
      for (j in 1:nTraps) {
        y[i, t, j] ~ dbern(z[i, t] * p)
      }
    }
  }
  
})

# Package data and constants
constants <- list(
  nSites = nSites,
  nYears = nYears,
  nTraps = nTraps
)

data <- list(y = y, # dimensions [sites, years, traps]
             CPUE = CPUE, # dimensions [sites, years + 1 year burnin]
             larv_R = larv_R, # no. released larvae [sites, years + 1 year burnin]
             larv_S = larv_S # no. settled larvae [sites_rel, sites_set, years + 1 year burnin]
             ) 

# connectivity inits
C <- array(NA, dim = c(nSites, nSites, nYears))
for (i in 1:nSites) {
  for (k in 1:nSites) {
    for (t in 1:nYears) {
      C[i, k, t] <- larv_S[i, k, t] / larv_R[i, t]
    }
  }
}

# Initial values
inits  <- function(y) {
  list(psi = runif(1, 0, 1),
       beta0 = runif(1, -1, 1),
       beta1 = runif(1, -1, 1),
       beta2 = runif(1, -1, 1),
       beta3 = runif(1, -1, 1),
       beta4 = runif(1, -1, 1),
       p = runif(1, 0, 1),
       z = apply(y, c(1, 2), max, na.rm = TRUE),
       C = C)
}

########################
# run MCMC in parallel #
########################

cl <- makeCluster(4)

set.seed(10120)

clusterExport(cl, c("model_code", "inits", "data", "constants", "y", "C"))

# parallelize running MCMC
out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  # Define nimbleFunctions directly on each worker
  colonization_lp <- nimbleFunction(
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
  assign("colonization_lp", colonization_lp, envir = .GlobalEnv)
  
  persistence_lp <- nimbleFunction(
    run = function(CPUE = double(0), beta3 = double(0), beta4 = double(0), 
                   C = double(0))
    {
      returnType(double(0))
      epsilon <- beta3 + CPUE * C * beta4
      return(epsilon)
    }
  )
  assign("persistence_lp", persistence_lp, envir = .GlobalEnv)
  
  # build model
  myModel <- nimbleModel(code = model_code,
                         data = data,
                         constants = constants,
                         inits = inits(y))
  
  
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
