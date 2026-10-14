library(dplyr)
library(tidyr)
library(purrr)
library(dad)
library(fitdistrplus)
source("R/split_train_test.R")

analysis_config <- config::get(file = "config.yaml")

data_all <- readRDS("data_heteromyid.rds")
newmoon_min <- min(data_all$newmoonnumber)
newmoon_max <- max(data_all$newmoonnumber)
train_win_width <- analysis_config$train_win_width

get_composition_distance <- function(
    comp_data_train,
    comp_data_test,
    test_start
) {
  get_probability_masses <- function(comp_data, split) {
    cleaned <- comp_data |> drop_na()
    tryCatch(
      {
        nb_fits <- cleaned |>
          group_by(species) |>
          summarize(fit = list(fitdist(abundance, "nbinom"))) |>
          mutate(
            size = map_dbl(fit, ~ .x$estimate["size"]),
            mu = map_dbl(fit, ~ .x$estimate["mu"])
          ) |>
          rowwise() |>
          mutate(prob_mass = list(dnbinom(0:200, mu = mu, size = size)))
        unlist(nb_fits$prob_mass) / nrow(nb_fits)
      },
      error = function(e) {
        species_abund <- cleaned |>
          group_by(species) |>
          summarize(vals = paste(abundance, collapse = ","), .groups = "drop")
        abund_str <- paste(
          paste0(species_abund$species, ": ", species_abund$vals),
          collapse = "; "
        )
        cat(
          format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"),
          sprintf(
            "test_start=%s, split=%s, abundances: %s\n",
            test_start,
            split,
            abund_str
          ),
          file = "comp_dist_poisson_fallbacks.log",
          append = TRUE
        )
        pois_fits <- cleaned |>
          group_by(species) |>
          summarize(fit = list(fitdist(abundance, "pois", method = "mme"))) |>
          mutate(lambda = map_dbl(fit, ~ .x$estimate["lambda"])) |>
          rowwise() |>
          mutate(prob_mass = list(dpois(0:200, lambda = lambda)))
        unlist(pois_fits$prob_mass) / nrow(pois_fits)
      }
    )
  }
  
  train_prob_vector <- get_probability_masses(comp_data_train, "train")
  test_prob_vector <- get_probability_masses(comp_data_test, "test")
  
  return(ddhellingerpar(train_prob_vector, test_prob_vector))
}


# Run all possible windows unless limited by config.yaml
if (is.null(analysis_config$test_starts)) {
  train_starts <- newmoon_min:(newmoon_max - train_win_width - 13 + 1)
} else {
  test_starts <- seq(
    from = analysis_config$test_starts$from,
    to = analysis_config$test_starts$to,
    by = analysis_config$test_starts$by %||% 1
  )
  train_starts <- test_starts - train_win_width
}

distances <- function(train_start) {

train_end <- train_start + train_win_width - 1
test_start <- train_end + 1
test_end <- test_start + 13 - 1

data_split <- split_train_test(
  data_all,
  gap = 0,
  train_start = train_start,
  train_end = train_end,
  test_start = test_start,
  test_end = test_end
)
  
  comp_data_train <- bind_cols(
    species = data_split$train$series,
    abundance = data_split$train$y
  )
  comp_data_test <- bind_cols(
    species = data_split$test$series,
    abundance = data_split$test$y
  )
  composition_distance <- data.frame(
    comp_dist = get_composition_distance(
      comp_data_train,
      comp_data_test,
      test_start
    )
  )
  composition_distance$test_start_newmoonnumber <- test_start
  composition_distance$species_list <- paste(data_split$species_list, collapse = "_")
  
  env_train <- data.frame(
    ndvi = data_split$train$ndvi,
    mintemp = data_split$train$meantemp_lag_1
  )
  env_test <- data.frame(
    ndvi = data_split$test$ndvi,
    mintemp = data_split$test$meantemp_lag_1
  )
  env_distance <- data.frame(enviro_dist = hellinger(env_train, env_test))
  env_distance$test_start_newmoonnumber <- test_start
  env_distance$species_list <- paste(data_split$species_list, collapse = "_")

  list(composition = composition_distance, env = env_distance)
}

quiet_distances <- quietly(distances)
all_distances <- map(train_starts, quiet_distances)

composition_distances <- map(all_distances, \(x) x$result$composition) |> list_rbind()
env_distances <- map(all_distances,  \(x) x$result$env) |> list_rbind()

saveRDS(env_distances, "env_distances.rds")
saveRDS(composition_distances, "composition_distances.rds")
