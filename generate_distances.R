# Calculate the environmental and composition distances between the training
# and test data of each forecast window set in config.yaml, without fitting any
# models. Saves env_distances.rds and composition_distances.rds.
# sliding_window_fits.R also does this at the end of a run when
# `calc_dist: yes`.
source("R/split_train_test.R")
source("R/windows.R")
source("R/distances.R")

analysis_config <- config::get(file = "config.yaml")
data_all <- readRDS("data_heteromyid.rds")

calculate_distances(data_all, analysis_config)
