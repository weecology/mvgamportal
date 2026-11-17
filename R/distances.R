library(dad)
library(dplyr)
library(fitdistrplus)
library(purrr)

get_probability_masses <- function(abundances, test_start, split) {
  tryCatch(
    {
      nb_fits <- abundances |>
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
      species_abund <- abundances |>
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
      pois_fits <- abundances |>
        group_by(species) |>
        summarize(fit = list(fitdist(abundance, "pois", method = "mme"))) |>
        mutate(lambda = map_dbl(fit, ~ .x$estimate["lambda"])) |>
        rowwise() |>
        mutate(prob_mass = list(dpois(0:200, lambda = lambda)))
      unlist(pois_fits$prob_mass) / nrow(pois_fits)
    }
  )
}

# Observed abundance of each species in the training or test data of a split
get_abundances <- function(data) {
  data.frame(species = data$series, abundance = data$y) |>
    filter(!is.na(abundance))
}

get_composition_distance <- function(data_split, test_start) {
  train <- get_abundances(data_split$train)
  test <- get_abundances(data_split$test)
  ddhellingerpar(
    get_probability_masses(train, test_start, "train"),
    get_probability_masses(test, test_start, "test")
  )
}

get_env_distance <- function(environment, window) {
  train <- environment |>
    filter(between(newmoonnumber, window$train_start, window$train_end)) |>
    dplyr::select(-newmoonnumber)
  test <- environment |>
    filter(between(newmoonnumber, window$test_start, window$test_end)) |>
    dplyr::select(-newmoonnumber)
  hellinger(train, test)
}

get_test_train_distances <- function(
  train_start,
  train_win_width,
  data_all,
  environment,
  dist_species
) {
  window <- get_window(train_start, train_win_width)
  env_distance <- data.frame(
    enviro_dist = get_env_distance(environment, window),
    test_start_newmoonnumber = window$test_start
  )

  data_split <- split_train_test(
    data_all,
    gap = 0,
    train_start = window$train_start,
    train_end = window$train_end,
    test_start = window$test_start,
    test_end = window$test_end,
    species = dist_species
  )
  if (length(data_split$species_list) == 0) {
    return(list(env = env_distance, composition = NULL))
  }
  composition_distance <- data.frame(
    comp_dist = get_composition_distance(data_split, window$test_start),
    test_start_newmoonnumber = window$test_start,
    species_list = paste(data_split$species_list, collapse = "_")
  )

  list(env = env_distance, composition = composition_distance)
}


calculate_distances <- function(data_all, analysis_config) {
  dist_species <- analysis_config$dist_species
  if (identical(dist_species, "all")) {
    dist_species <- levels(data_all$series)
  }
  environment <- data.frame(
    newmoonnumber = data_all$newmoonnumber,
    ndvi = as.vector(data_all$ndvi),
    meantemp_lag_1 = as.vector(data_all$meantemp_lag_1)
  ) |>
    distinct(newmoonnumber, .keep_all = TRUE)

  train_starts <- get_train_starts(data_all$newmoonnumber, analysis_config)
  safe_test_train_distances <- purrr::safely(get_test_train_distances)
  results <- purrr::map(
    train_starts,
    \(train_start) safe_test_train_distances(
      train_start,
      analysis_config$train_win_width,
      data_all,
      environment,
      dist_species
    ),
    .progress = TRUE
  ) |>
    handle_failed_windows(train_starts, "Distances")

  saveRDS(purrr::map_dfr(results, "env"), "env_distances.rds")
  saveRDS(purrr::map_dfr(results, "composition"), "composition_distances.rds")
}
