#!/usr/bin/env Rscript
#
# Fit one model to one window, so a new model can be checked on its own before
# it is run across every window by sliding_window_fits.R.
#
# It goes through the same steps a real run does — split_train_test(),
# fit_model(), evaluate_model(), tidy_score(), plot_run_forecasts() — and
# reports the window, sampler diagnostics and CRPS, so a model that samples
# badly or cannot forecast shows up here rather than after a full run.
#
# From the shell:
#   Rscript test_model.R AR_SINGLE_SPECIES --species DM --test-start 240
#   Rscript test_model.R GAM_AR --train-win-width 120 --trace
#
# From an R session, which also gives back the fitted model to poke at:
#   source("test_model.R")
#   fit <- test_model("AR_SINGLE_SPECIES", species = "DM")
#
# A model that has an R/models/<NAME>.R file but is not yet wired into
# fit_model() or config.yaml can still be tested: it is called directly and
# the report says what is left to do to integrate it.

library(dplyr)
library(glue)
library(mvgam)
library(purrr)
library(tidyr)
source("R/models.R")
source("R/skill_scores.r")
source("R/forecast_plots.r")
source("R/split_train_test.R")

# Fit and score `model` on the single window ending at `test_start`. Defaults
# for the window come from config.yaml, so a test matches how the model would
# be run for real. Returns the fitted model, its scores and the split, invisibly.
test_model <- function(
  model,
  species = NULL,
  test_start = NULL,
  train_win_width = NULL,
  data_file = "data_heteromyid.rds",
  plot = TRUE,
  trace_plot = FALSE,
  horizon = 12
) {
  analysis_config <- config::get(file = "config.yaml")
  train_win_width <- train_win_width %||% analysis_config$train_win_width
  data_all <- readRDS(data_file)
  test_start <- test_start %||%
    analysis_config$test_starts$from %||%
    (min(data_all$newmoonnumber) + train_win_width)

  train_start <- test_start - train_win_width
  train_end <- test_start - 1
  test_end <- test_start + horizon - 1

  split <- split_train_test(
    data_all,
    gap = 0,
    train_start = train_start,
    train_end = train_end,
    test_start = test_start,
    test_end = test_end,
    species = species
  )
  species_list <- split$species_list
  if (length(species_list) == 0) {
    stop(glue(
      "No species clear the 30% occupancy filter in newmoons ",
      "{train_start}-{train_end}; sliding_window_fits.R would skip this run."
    ))
  }

  # Labelled apart from the real runs so test plots do not overwrite theirs
  run_label <- paste0(model, "_test")
  # With one species score()'s all_series element is a copy of that species'
  # score, so drop it the way a `separately` run does
  single_species <- length(species_list) == 1

  cat(sprintf(
    "\nFitting %s to %s, test start %d\n",
    model,
    paste(species_list, collapse = ", "),
    test_start
  ))
  elapsed <- system.time(
    fit <- dispatch_fit(model, split$train, split$test, species_list)
  )[["elapsed"]]

  evaluation <- evaluate_model(
    fit,
    model,
    species_set = "test",
    test_start = test_start,
    species_list = species_list
  )
  scores <- tidy_score(evaluation$score, model, single_species)

  if (plot) {
    plot_run_forecasts(
      fit,
      run_label,
      scores,
      species_list,
      test_start,
      trace_plots = trace_plot
    )
  }

  report_test(
    model,
    run_label,
    data_file,
    species,
    species_list,
    train_start,
    train_end,
    test_start,
    test_end,
    scores,
    elapsed,
    plot,
    trace_plot
  )

  invisible(list(
    model = fit,
    scores = scores,
    summary = evaluation$summary,
    loo = evaluation$loo,
    split = split
  ))
}

# Prefer fit_model(), so the branch a real run goes through is what gets
# tested, but fall back to the fit function itself for a model that has not
# been added to fit_model() yet.
dispatch_fit <- function(model, data_train, data_test, species_list) {
  if (has_fit_branch(model)) {
    fit_model(model, data_train, data_test, species_list)
  } else {
    get(paste0("fit_", model), mode = "function")(
      data_train,
      data_test,
      species_list
    )
  }
}

# fit_model() is a single switch() over model names (see CLAUDE.md), so its
# branches are the names of that call's arguments.
has_fit_branch <- function(model) {
  model %in% names(as.list(body(fit_model)[[2]]))
}

report_test <- function(
  model,
  run_label,
  data_file,
  requested_species,
  species_list,
  train_start,
  train_end,
  test_start,
  test_end,
  scores,
  elapsed,
  plot,
  trace_plot
) {
  requested <- requested_species %||% "all"
  dropped <- setdiff(requested_species, as.character(species_list))
  # Newmoons with no observation are scored NA, so say how many of the horizon
  # each mean is actually over
  crps <- scores |>
    filter(species != "all_series") |>
    group_by(species) |>
    summarise(
      mean_crps = mean(score, na.rm = TRUE),
      n_scored = sum(!is.na(score)),
      n = n(),
      .groups = "drop"
    )

  field <- function(name, ...) cat(sprintf("  %-14s%s\n", name, paste0(...)))
  cat(sprintf("\n=== %s ===\n", model))
  field("data", data_file)
  field("train", sprintf(
    "newmoon %d-%d (%d newmoons)",
    train_start,
    train_end,
    train_end - train_start + 1
  ))
  field("test", sprintf("newmoon %d-%d", test_start, test_end))
  field(
    "species",
    paste(species_list, collapse = ", "),
    sprintf(" (requested %s", paste(requested, collapse = ", ")),
    if (length(dropped)) {
      sprintf("; %s below 30%% occupancy", paste(dropped, collapse = ", "))
    },
    ")"
  )
  field("fit time", sprintf("%.1f min", elapsed / 60))

  cat("\n  Sampler\n")
  field("mean rhat", round(scores$rhat[1], 4))
  field("rhat > 1.05", sprintf("%.1f%% of parameters", 100 * scores$prhat_high[1]))
  field("divergences", scores$n_divergences[1])

  cat("\n  Mean CRPS over the forecast horizon\n")
  pwalk(crps, \(species, mean_crps, n_scored, n) field(
    species,
    sprintf("%.3f (%d of %d newmoons observed)", mean_crps, n_scored, n)
  ))

  if (plot) {
    cat("\n  Plots\n")
    field("forecasts", sprintf("figures/%s_%d.png", run_label, test_start))
    if (trace_plot) {
      field("trace", sprintf("figures/%s_trace_%d.png", run_label, test_start))
    }
  }

  report_integration(model)
  cat("\n")
}

# What is still needed before sliding_window_fits.R can run this model.
report_integration <- function(model) {
  configured <- map_chr(config::get(file = "config.yaml")$models, "name")
  todo <- c(
    if (!has_fit_branch(model)) {
      glue(
        "add `{model} = fit_{model}(data_train, data_test, species_list),` ",
        "to fit_model() in R/models.R"
      )
    },
    if (!model %in% configured) {
      glue("add `- name: {model}` to the models list in config.yaml")
    }
  )
  if (length(todo) == 0) {
    return(invisible(NULL))
  }
  cat("\n  To run this model in sliding_window_fits.R\n")
  walk(todo, \(step) cat(sprintf("    - %s\n", step)))
}

parse_args <- function(args) {
  opts <- list()
  i <- 1
  while (i <= length(args)) {
    switch(
      args[i],
      "--trace" = {
        opts$trace_plot <- TRUE
        i <- i + 1
      },
      "--no-plot" = {
        opts$plot <- FALSE
        i <- i + 1
      },
      if (grepl("^--", args[i])) {
        opts[[gsub("-", "_", sub("^--", "", args[i]))]] <- args[i + 1]
        i <- i + 2
      } else {
        opts$model <- args[i]
        i <- i + 1
      }
    )
  }
  opts$species <- if (!is.null(opts$species)) strsplit(opts$species, ",")[[1]]
  for (n in intersect(names(opts), c("test_start", "train_win_width", "horizon"))) {
    opts[[n]] <- as.integer(opts[[n]])
  }
  opts
}

args <- commandArgs(trailingOnly = TRUE)
if (!interactive() && length(args) > 0) {
  do.call(test_model, parse_args(args))
}
