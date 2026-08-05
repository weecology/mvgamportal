library(dplyr)
library(cmdstanr)
library(dad)
library(fitdistrplus)
library(forecast)
library(mvgam)
library(portalr)
library(glue)
library(furrr)
library(purrr)
library(tidyr)
source("R/get_regime.R")
source("R/skill_scores.r")
source("R/models.R")
source("R/forecast_plots.r")

analysis_config <- config::get(file = "config.yaml")
model_config <- analysis_config$models
names(model_config) <- purrr::map_chr(model_config, "name")
check_model_definitions(names(model_config))

# Turn the configured species groups into one entry per set of species that
# gets fit together. A group with `separately` becomes one entry per species,
# labelled by species name rather than by group name.
expand_species_sets <- function(species_sets, all_species) {
  expanded <- purrr::flatten(purrr::map(species_sets, function(species_set) {
    species <- species_set$species %||% all_species
    if (!is.null(species_set$separately)) {
      purrr::map(species, \(focal_species) list(
        name = focal_species,
        group = species_set$name,
        species = focal_species,
        separately = TRUE
      ))
    } else {
      list(list(
        name = species_set$name,
        group = species_set$name,
        species = species,
        separately = FALSE
      ))
    }
  }))

  purrr::set_names(expanded, purrr::map_chr(expanded, "name"))
}

# Convert config info into models to be run.
# One model per model name (from config) per species group.
# Models name the groups from config.yaml and`separately` expands
# into one run per species.
build_run_plan <- function(model_config, species_sets) {
  groups <- purrr::map_chr(species_sets, "group")
  purrr::flatten(purrr::map(model_config, function(model) {
    set_names <- names(species_sets)[
      groups %in% (model$species_sets %||% groups)
    ]
    purrr::map(set_names, \(set_name) list(
      model = model$name,
      species_set = set_name,
      separately = species_sets[[set_name]]$separately,
      label = paste(model$name, set_name, sep = "_")
    ))
  }))
}

# Set number of workers. Each worker spawns 4 cmdstanr chains, so total
# cores ~= n_workers * 4. To run on HPC set MVGAM_N_WORKERS environmental
# variable based on number of requested cores
n_workers <- as.integer(Sys.getenv("MVGAM_N_WORKERS", unset = "8"))
plan(multisession, workers = n_workers / 4)

data_all <- readRDS("data_heteromyid.rds")

species_sets <- expand_species_sets(
  analysis_config$species_sets,
  levels(data_all$series)
)
run_plan <- build_run_plan(model_config, species_sets)

# Drop species group groups that no model in the config uses
# This prevents their splits and distances from being computed for no reason
species_sets <- species_sets[unique(purrr::map_chr(run_plan, "species_set"))]

split_train_test <- function(
  data_all,
  gap,
  train_start,
  train_end,
  test_start,
  test_end,
  species = NULL
) {
  species <- species %||% levels(data_all$series)
  species_list <- data.frame(
    newmoonnumber = data_all$newmoonnumber,
    series = data_all$series,
    y = data_all$y
  ) |>
    filter(
      newmoonnumber >= train_start,
      newmoonnumber <= (train_end),
      series %in% species
    ) |>
    group_by(newmoonnumber, series) |>
    summarise(abundance = sum(y, na.rm = TRUE), .groups = 'drop') |>
    group_by(series) |>
    summarise(occupancy = sum(abundance > 0) / n(), .groups = 'drop') |>
    filter(occupancy >= 0.30) |>
    pull(series) |>
    droplevels()

  train_inds <- which(
    data_all$newmoonnumber >= train_start &
      data_all$newmoonnumber <= (train_end) &
      data_all$series %in% species_list
  )
  data_train <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][train_inds, ]
    } else {
      data_all[[x]][train_inds]
    }
  })

  names(data_train) <- names(data_all)
  data_train$series <- droplevels(data_train$series)

  test_inds <- which(
    data_all$newmoonnumber >= (test_start) &
      data_all$newmoonnumber <= test_end &
      data_all$series %in% species_list
  )
  data_test <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][test_inds, ]
    } else {
      data_all[[x]][test_inds]
    }
  })

  names(data_test) <- names(data_all)
  data_test$series <- droplevels(data_test$series)

  return(list(
    train = data_train,
    test = data_test,
    species_list = species_list
  ))
}

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

newmoon_min <- min(data_all$newmoonnumber)
newmoon_max <- max(data_all$newmoonnumber)
train_win_width <- analysis_config$train_win_width

# Run all possible windows unless limited by config.yaml
if (is.null(analysis_config$test_starts)) {
  train_starts <- newmoon_min:(newmoon_max - train_win_width - 12 + 1)
} else {
  test_starts <- seq(
    from = analysis_config$test_starts$from,
    to = analysis_config$test_starts$to,
    by = analysis_config$test_starts$by %||% 1
  )
  train_starts <- test_starts - train_win_width
}

run_window <- function(train_start) {
  train_end <- train_start + train_win_width - 1
  test_start <- train_end + 1
  test_end <- test_start + 12 - 1
  print(glue("Training test start {test_start} (newmoon)"))

  data_splits <- purrr::map(species_sets, \(species_set) split_train_test(
    data_all,
    gap = 0,
    train_start = train_start,
    train_end = train_end,
    test_start = test_start,
    test_end = test_end,
    species = species_set$species
  ))

  # Fit, evaluate and plot one run at a time so that only a single model is
  # held in memory at once
  evaluations <- purrr::map(run_plan, function(run) {
    data_split <- data_splits[[run$species_set]]
    model <- fit_model(
      run$model,
      data_split$train,
      data_split$test,
      data_split$species_list
    )
    evaluation <- evaluate_model(
      model,
      run$model,
      run$species_set,
      test_start,
      data_split$species_list
    )
    plot_run_forecasts(
      model,
      run$label,
      tidy_score(evaluation$score, run$model, run$separately),
      data_split$species_list,
      test_start,
      trace_plots = analysis_config$trace_plot
    )
    rm(model)
    gc()
    evaluation
  }) |>
    purrr::set_names(purrr::map_chr(run_plan, "label"))

  distances <- purrr::imap(data_splits, function(data_split, set_name) {
    data_train <- data_split$train
    data_test <- data_split$test
    species_str <- paste(data_split$species_list, collapse = "_")

    comp_data_train <- bind_cols(
      species = data_train$series,
      abundance = data_train$y
    )
    comp_data_test <- bind_cols(
      species = data_test$series,
      abundance = data_test$y
    )
    composition_distance <- data.frame(
      comp_dist = get_composition_distance(
        comp_data_train,
        comp_data_test,
        test_start
      )
    )
    composition_distance$test_start_newmoonnumber <- test_start
    composition_distance$species_set <- set_name
    composition_distance$species_list <- species_str

    env_train <- data.frame(
      ndvi = data_train$ndvi,
      mintemp = data_train$meantemp_lag_1
    )
    env_test <- data.frame(
      ndvi = data_test$ndvi,
      mintemp = data_test$meantemp_lag_1
    )
    env_distance <- data.frame(enviro_dist = hellinger(env_train, env_test))
    env_distance$test_start_newmoonnumber <- test_start
    env_distance$species_set <- set_name
    env_distance$species_list <- species_str

    list(composition = composition_distance, env = env_distance)
  })

  # Keyed by model rather than run label so BASELINE is recognisable; runs are
  # kept apart by the species_set column carried on each score
  scores <- purrr::map(evaluations, "score") |>
    purrr::set_names(purrr::map_chr(run_plan, "model")) |>
    get_skill_scores(purrr::map_lgl(run_plan, "separately"))

  gc()

  list(
    scores = scores,
    summaries = purrr::map(evaluations, "summary"),
    loos = purrr::map(evaluations, "loo"),
    env_distance = purrr::map_dfr(distances, "env"),
    composition_distance = purrr::map_dfr(distances, "composition")
  )
}

safe_run_window <- purrr::safely(run_window)
results <- future_map(
  train_starts,
  safe_run_window,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

errored <- purrr::map_lgl(results, ~ !is.null(.x$error))
if (any(errored)) {
  warning(glue(
    "{sum(errored)} window(s) failed: train_starts {paste(train_starts[errored], collapse=', ')}"
  ))
}
results <- purrr::map(results, "result")

scores <- purrr::map_dfr(results, "scores")
summaries <- purrr::flatten(purrr::map(results, "summaries"))
loos <- purrr::flatten(purrr::map(results, "loos"))
env_distances <- purrr::map(results, "env_distance")
composition_distances <- purrr::map(results, "composition_distance")

saveRDS(scores, "scores.rds")
saveRDS(summaries, "summaries.rds")
saveRDS(loos, "loos.rds")
saveRDS(env_distances, "env_distances.rds")
saveRDS(composition_distances, "composition_distances.rds")
