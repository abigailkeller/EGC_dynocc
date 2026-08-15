library(MCMCvis)
library(tidyverse)

# read in samples
samples <- readRDS("data/posterior_samples/posterior_samples.rds")

# summarize
summary <- MCMCsummary(samples)

# trace plot
param <- "z[73, 6]"

ggplot() +
  geom_line(aes(x = 1:nrow(out_sub[[1]]), y = out_sub[[1]][, param]),
            color = "blue") +
  geom_line(aes(x = 1:nrow(out_sub[[2]]), y = out_sub[[2]][, param]),
            color = "red") +
  geom_line(aes(x = 1:nrow(out_sub[[3]]), y = out_sub[[3]][, param]),
            color = "purple") +
  geom_line(aes(x = 1:nrow(out_sub[[4]]), y = out_sub[[4]][, param]),
            color = "green")


## Note: some Rhat > 1.1 (will discuss)
not_converged <- summary[which(summary$Rhat > 1.1), ]
