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
source("R/split_train_test.R")

analysis_config <- config::get(file = "config.yaml")
model_config <- analysis_config$models
names(model_config) <- purrr::map_chr(model_config, "name")
# BASELINE is fit for every species set rather than configured, so check it too
check_model_definitions(c("BASELINE", names(model_config)))

# Turn the configured species groups into one entry per set of species that
# gets fit together. A group with `separately: yes` becomes one entry per
# species, labelled by species name rather than by group name.
expand_species_sets <- function(species_sets, all_species) {
  expanded <- purrr::flatten(purrr::map(species_sets, function(species_set) {
    species <- species_set$species %||% all_species
    if (species_set$separately %||% FALSE) {
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
# into one run per species. `existing_runs` sets whether to rerun
# results already in the rds files and replace them or skip and keep
# the old values.
build_run_plan <- function(model_config, species_sets, existing_runs) {
  groups <- purrr::map_chr(species_sets, "group")
  make_run <- function(model_name, set_name, existing_runs) list(
    model = model_name,
    species_set = set_name,
    separately = species_sets[[set_name]]$separately,
    label = paste(model_name, set_name, sep = "_"),
    existing_runs = existing_runs
  )

  runs <- purrr::flatten(purrr::map(model_config, function(model) {
    set_names <- names(species_sets)[
      groups %in% (model$species_sets %||% groups)
    ]
    purrr::map(set_names, \(set_name) make_run(
      model$name,
      set_name,
      model$existing_runs %||% existing_runs
    ))
  }))

  baseline_runs <- purrr::map(
    unique(purrr::map_chr(runs, "species_set")),
    \(set_name) make_run("BASELINE", set_name, existing_runs)
  )
  c(baseline_runs, runs)
}

# Set number of workers. Each worker spawns 4 cmdstanr chains, so total
# cores ~= n_workers * 4. To run on HPC set MVGAM_N_WORKERS environmental
# variable based on number of requested cores
n_workers <- as.integer(Sys.getenv("MVGAM_N_WORKERS", unset = "20"))
plan(multisession, workers = n_workers / 4)

data_all <- readRDS("data_heteromyid.rds")

species_sets <- expand_species_sets(
  analysis_config$species_sets,
  levels(data_all$series)
)
run_plan <- build_run_plan(
  model_config,
  species_sets,
  analysis_config$existing_runs %||% "skip"
)

# Drop species group groups that no model in the config uses
# This prevents their splits and distances from being computed for no reason
species_sets <- species_sets[unique(purrr::map_chr(run_plan, "species_set"))]

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
  train_starts <- newmoon_min:(newmoon_max - train_win_width - 13 + 1)
} else {
  test_starts <- seq(
    from = analysis_config$test_starts$from,
    to = analysis_config$test_starts$to,
    by = analysis_config$test_starts$by %||% 1
  )
  train_starts <- test_starts - train_win_width
}

output_files <- c(
  scores = "scores.rds",
  summaries = "summaries.rds",
  loos = "loos.rds",
  env_distance = "env_distances.rds",
  composition_distance = "composition_distances.rds"
)
read_results <- function(path, empty) {
  if (file.exists(path)) readRDS(path) else empty
}
no_distances <- tibble(
  species_set = character(),
  test_start_newmoonnumber = numeric()
)
no_scores <- mutate(no_distances, model = character())
existing <- list(
  # Skill scores are recomputed once everything is combined
  scores = read_results(output_files[["scores"]], no_scores) |>
    dplyr::select(-any_of(c("baseline_score", "skill_score"))),
  summaries = read_results(output_files[["summaries"]], list()),
  loos = read_results(output_files[["loos"]], list()),
  # bind_rows() also reads distances saved as one data frame per window
  env_distance = bind_rows(
    read_results(output_files[["env_distance"]], no_distances)
  ),
  composition_distance = bind_rows(
    read_results(output_files[["composition_distance"]], no_distances)
  )
)

# Identifies a run across windows, for matching runs to the saved results
run_id <- function(model, species_set, test_start) {
  paste(model, species_set, test_start, sep = "|")
}
completed_runs <- with(
  existing$scores,
  run_id(model, species_set, test_start_newmoonnumber)
)

# Only keep runs that are flagged to overwrite or are not present in RDS files
window_plans <- purrr::map(train_starts + train_win_width, function(test_start) {
  purrr::keep(run_plan, \(run) {
    run$existing_runs == "overwrite" ||
      !run_id(run$model, run$species_set, test_start) %in% completed_runs
  })
})
has_runs <- lengths(window_plans) > 0
train_starts <- train_starts[has_runs]
window_plans <- window_plans[has_runs]

run_window <- function(train_start, window_plan) {
  train_end <- train_start + train_win_width - 1
  test_start <- train_end + 1
  test_end <- test_start + 13 - 1
  print(glue("Training test start {test_start} (newmoon)"))

  window_sets <- species_sets[unique(purrr::map_chr(window_plan, "species_set"))]
  data_splits <- purrr::map(window_sets, \(species_set) split_train_test(
    data_all,
    gap = 0,
    train_start = train_start,
    train_end = train_end,
    test_start = test_start,
    test_end = test_end,
    species = species_set$species
  ))

  # Fit, evaluate and plot one run at a time so that only a single model is
  # held in memory at once. A run with nothing left after the occupancy filter
  # is dropped, or a run that fails, is scored as NAs, so the rest of the nb_fits
  # can complete successfully.
  evaluations <- purrr::map(window_plan, function(run) {
    data_split <- data_splits[[run$species_set]]
    if (length(data_split$species_list) == 0) {
      return(NULL)
    }

    evaluation <- tryCatch(
      {
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
        evaluation
      },
      error = function(e) {
        message(glue(
          "{run$label} failed for test start {test_start}: {conditionMessage(e)}"
        ))
        list(
          score = na_scores(
            data_split$species_list,
            run$species_set,
            test_start
          ),
          summary = NULL,
          loo = NULL
        )
      }
    )
    gc()
    evaluation
  })

  # Keep the run plan aligned with the runs that produced something
  fitted_runs <- window_plan[!purrr::map_lgl(evaluations, is.null)]
  evaluations <- purrr::compact(evaluations) |>
    purrr::set_names(purrr::map_chr(fitted_runs, "label"))

  distances <- purrr::imap(data_splits, function(data_split, set_name) {
    if (length(data_split$species_list) == 0) {
      return(NULL)
    }
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

  scores <- purrr::map2_dfr(
    evaluations,
    fitted_runs,
    \(evaluation, run) tidy_score(evaluation$score, run$model, run$separately)
  )

  gc()

  list(
    # List all runs attempted. Ensures replacement if run produced no rows
    attempted = purrr::map_chr(
      window_plan,
      \(run) run_id(run$model, run$species_set, test_start)
    ),
    scores = scores,
    summaries = purrr::map(evaluations, "summary"),
    loos = purrr::map(evaluations, "loo"),
    env_distance = purrr::map_dfr(distances, "env"),
    composition_distance = purrr::map_dfr(distances, "composition")
  )
}

safe_run_window <- purrr::safely(run_window)
results <- future_map2(
  train_starts,
  window_plans,
  safe_run_window,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

errored <- purrr::map_lgl(results, ~ !is.null(.x$error))
if (any(errored)) {
  purrr::walk2(
    train_starts[errored],
    purrr::map(results[errored], "error"),
    \(train_start, error) message(glue(
      "train_start {train_start} failed: {conditionMessage(error)}"
    ))
  )
  warning(glue(
    "{sum(errored)} window(s) failed: train_starts {paste(train_starts[errored], collapse=', ')}"
  ))
}
results <- purrr::map(results, "result")

# Update results with new runs
attempted <- unlist(purrr::map(results, "attempted"))
keep_unreplaced <- function(results) {
  purrr::discard(
    purrr::compact(results),
    \(result) run_id(
      result$model,
      result$species_set,
      result$test_start_newmoonnumber
    ) %in% attempted
  )
}
# Distances depend only on the window and species set, so new ones replace old
update_distances <- function(existing, new) {
  bind_rows(
    filter(
      existing,
      !paste(species_set, test_start_newmoonnumber) %in%
        paste(new[["species_set"]], new[["test_start_newmoonnumber"]])
    ),
    new
  )
}

scores <- bind_rows(
  filter(
    existing$scores,
    !run_id(model, species_set, test_start_newmoonnumber) %in% attempted
  ),
  purrr::map_dfr(results, "scores")
) |>
  get_skill_scores()
summaries <- c(
  keep_unreplaced(existing$summaries),
  purrr::compact(purrr::flatten(purrr::map(results, "summaries")))
)
loos <- c(
  keep_unreplaced(existing$loos),
  purrr::compact(purrr::flatten(purrr::map(results, "loos")))
)
env_distances <- update_distances(
  existing$env_distance,
  purrr::map_dfr(results, "env_distance")
)
composition_distances <- update_distances(
  existing$composition_distance,
  purrr::map_dfr(results, "composition_distance")
)

saveRDS(scores, output_files[["scores"]])
saveRDS(summaries, output_files[["summaries"]])
saveRDS(loos, output_files[["loos"]])
saveRDS(env_distances, output_files[["env_distance"]])
saveRDS(composition_distances, output_files[["composition_distance"]])
