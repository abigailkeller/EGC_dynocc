set.seed(42)

#################################
# functions for simulating data #
#################################

# calculate probability of colonization
calc_gamma <- nimbleFunction (
  # input and output types
  run = function(CPUE = double(2), beta0 = double(0), beta1 = double(0), 
                 beta2 = double(0), C = double(2), index = double(0))
  {
    returnType(double(0))
    
    idx <- as.integer(index)
    n <- dim(CPUE)[1]
    gamma <- 0
    
    for (j in 1:n) {
      if (j != idx) {
        gamma <- gamma + CPUE[j, 2] * beta1 * C[j, 2] +  # t - 1
          CPUE[j, 1] * beta2 * C[j, 1]    # t - 2
      }
    }
    
    gamma <- gamma + beta0 
    
    return(gamma)
  }
)

# calculate probability of persistence
calc_epsilon <- nimbleFunction (
  # input and output types
  run = function(CPUE = double(0), beta3 = double(0), beta4 = double(0), 
                 C = double(0))
  {
    returnType(double(0))
    
    epsilon <- beta3 + CPUE * C * beta4
    
    return(epsilon)
  }
)


# constants
nSites  <- 20
nYears  <- 5
nTraps <- 3

# CPUE
CPUE <- matrix(runif(nYears * nSites, 0, 10), 
               nrow = nSites, ncol = nYears)
# larval data
# released larvae
larv_R <- matrix(100, nrow = nSites, ncol = nYears)
# settled larvae
larv_S <- array(NA, dim = c(nSites, nSites, nYears ))
for (i in 1:nSites) {
  for (t in 1:nYears) {
    larv_S[i, , t] <- rmultinom(1, size = 100, prob = rep(1 / nSites, nSites))
  }
}

# params
psi <- 0.5
epsilon <- 0.8
p <- 0.7
beta0 <- 0
beta1 <- 0.1
beta2 <- 0.02
beta3 <- 0
beta4 <- 0.1

# calculate connectivity
C <- array(NA, dim = c(nSites, nSites, nYears))
for (i in 1:nSites) {
  for (k in 1:nSites) {
    for (t in 1:nYears) {
      C[i, k, t] <- larv_S[i, k, t] / larv_R[i, t]
    }
  }
}

# calculate gamma
gamma <- matrix(NA, nrow = nSites, ncol = nYears - 1)
for (i in 1:nSites) {
  for (t in 2:nYears) {
    gamma[i, t - 1] <- calc_gamma(CPUE[1:nSites, (t-1):t], 
                                       beta0, beta1, beta2,
                                       C[i, 1:nSites, (t-1):t], 
                                       i)
  }
}
# 
gamma_trans <- plogis(gamma)

# calculate probability of persistence
epsilon <- matrix(NA, nrow = nSites, ncol = nYears - 1)
for (i in 1:nSites) {
  for (t in 2:nYears) {
    epsilon[i, t - 1] <- calc_epsilon(CPUE[i, t], 
                                        beta3, beta4,
                                        C[i, i, t])
  }
}
# 
epsilon_trans <- plogis(epsilon)

# Simulate occupancy + detections
z <- matrix(NA, nSites, nYears)
y <- array(NA, dim = c(nSites, nYears, nTraps))
z[, 1] <- rbinom(nSites, 1, psi)

for (t in 2:nYears) {
  for (i in 1:nSites) {
    z[, t] <- rbinom(nSites, 1, z[, t - 1] * epsilon_trans[, t - 1] + 
                       (1 - z[, t - 1]) * gamma_trans[, t - 1])
  }
}
  
for (j in 1:nTraps) {
  y[, , j] <- rbinom(nSites * nYears, 1, z * p) 
}
