library(MCMCvis)
library(tidyverse)
library(patchwork)

# read in samples
samples <- readRDS("data/posterior_samples/posterior_samples_20260827.rds")

# summarize
summary <- MCMCsummary(samples)

# trace plot
param <- "beta2"

ggplot() +
  geom_line(aes(x = 1:nrow(samples[[1]]), y = samples[[1]][, param]),
            color = "blue") +
  geom_line(aes(x = 1:nrow(samples[[2]]), y = samples[[2]][, param]),
            color = "red") +
  geom_line(aes(x = 1:nrow(samples[[3]]), y = samples[[3]][, param]),
            color = "purple") +
  geom_line(aes(x = 1:nrow(samples[[4]]), y = samples[[4]][, param]),
            color = "green") + 
  geom_line(aes(x = 1:nrow(samples[[5]]), y = samples[[5]][, param]),
            color = "black") +
  geom_line(aes(x = 1:nrow(samples[[6]]), y = samples[[6]][, param]),
            color = "yellow") +
  geom_line(aes(x = 1:nrow(samples[[7]]), y = samples[[7]][, param]),
            color = "brown") +
  geom_line(aes(x = 1:nrow(samples[[8]]), y = samples[[8]][, param]),
            color = "orange") +
  labs(x = "iteration", y = "value") +
  ggtitle(param)

param1 <- "beta3"
param2 <- "beta4"

ggplot() + 
  geom_point(aes(x = c(samples[[1]][, param1]),
                 y = c(samples[[1]][, param2])))


# 


## Note: some Rhat > 1.1 (will discuss)
not_converged <- summary[which(summary$Rhat > 1.1), ]

# summarize prob detection
summary[c("p_minnow", "p_shrimp", "p_fukui"), ]

# summarize beta
summary[c("beta0", "beta1", "beta2", "beta3", "beta4"), ]

# read in detection data
detections <- readRDS("data/model_data/PresenceArrayBinary.rds")[, 2:7, ]
remove <- as.integer(which(apply(detections, 1, function(s) all(is.na(s)))))
detections <- detections[-remove, , ]

# example posterior of a z
ggplot() +
  geom_histogram(aes(x = c(samples[[1]][, "z[74, 1]"], 
                           samples[[2]][, "z[74, 1]"],
                           samples[[3]][, "z[74, 1]"],
                           samples[[4]][, "z[74, 1]"]))) +
  labs(x = "z, site 74, year 1", y = "count") +
  theme_minimal()

# look at specific sites
nicks_lagoon <- rowSums(detections[74, , ], na.rm = TRUE)
seabeck_conf <- rowSums(detections[103, , ], na.rm = TRUE)
seabeck_marina <- rowSums(detections[104, , ], na.rm = TRUE)
ind_island <- rowSums(detections[58, , ], na.rm = TRUE)
jimmy <- rowSums(detections[60, , ], na.rm = TRUE)

site_ind <- 74
summary[c(paste0("z[", site_ind, ", 1]"), paste0("z[", site_ind, ", 2]"),
          paste0("z[", site_ind, ", 3]"), paste0("z[", site_ind, ", 4]"), 
          paste0("z[", site_ind, ", 5]"), paste0("z[", site_ind, ", 6]")), ]

# function to plot observations and posterior means
plot_z <- function(index, name) {
  
  det_summary <- rowSums(detections[index, , ], na.rm = TRUE)
  
  detections <- ggplot() +
    geom_point(aes(x = 1:6, y = det_summary)) +
    labs(x = "time", y = "# detections") +
    ggtitle(name) +
    theme_minimal()
  
  post_summary <- summary[c(paste0("z[", index, ", 1]"), 
                            paste0("z[", index, ", 2]"),
                            paste0("z[", index, ", 3]"), 
                            paste0("z[", index, ", 4]"), 
                            paste0("z[", index, ", 5]"), 
                            paste0("z[", index, ", 6]")), ]
  
  post <- ggplot() +
    geom_point(aes(x = 1:6, y = post_summary[, "mean"])) +
    geom_errorbar(aes(x = 1:6, ymin = post_summary[, "2.5%"],
                      ymax = post_summary[, "97.5%"])) +
    labs(x = "time", y = "z (occupancy state)") +
    theme_minimal()
  
  return(detections + post + plot_layout(ncol = 1))
    
}

plot_z(74, "Nick's Laggon")
plot_z(103, "Seabeck Conference")
plot_z(104, "Seabeck Marina")
plot_z(58, "Indian Island")
plot_z(60, "Jimmycomelately")


