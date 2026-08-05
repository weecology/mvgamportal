library(mvgam)

# Shared model settings and output handling for the sliding window analysis.
# Which models are run is set in config.yaml.
#
# Each model lives in its own file in R/models/, which defines a single
# fit_<NAME>() function. To add a model: write R/models/<NAME>.R defining
# fit_<NAME>(data_train, data_test, species_list), add a branch to fit_model()
# below, and add <NAME> to the models list in config.yaml.
invisible(lapply(
  list.files("R/models", pattern = "\\.R$", full.names = TRUE),
  source
))

# Priors

# sigma_prior <- prior(beta(10, 10), class = sigma, lb = 0.2, ub = 1)

sigma_prior <- c(prior(exponential(5), class = sigma))

ndvi_random_slopes_prior <- prior(
  inv_gamma(2.3693353, 0.7311319),
  class = sigma_raw_trend
)

ar_priors <- c(sigma_prior)
gam_ar_priors <- c(sigma_prior)
gam_var_priors <- c(sigma_prior)

# Species given their own temperature response on top of the shared one. Only
# those also present in the window get a term.
weighted_species <- c("DM", "DO", "PP")

# Built per window because which species are present varies. A single-species
# window gets no by-species terms, since the shared term already describes it.
get_trend_formula <- function(species_list) {
  species <- levels(species_list)
  weighted <- if (length(species) > 1) {
    intersect(weighted_species, species)
  } else {
    character(0)
  }

  reformulate(c(
    's(ndvi_ma12, trend, bs = "re")',
    'te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr"))',
    glue::glue(
      'te(mintemp, lag, by = weights_{tolower(weighted)}, ',
      'k = c(3, 4), bs = c("tp", "cr"))'
    )
  ))
}

fit_model <- function(model_name, data_train, data_test, species_list) {
  switch(
    model_name,
    BASELINE = fit_BASELINE(data_train, data_test, species_list),
    AR = fit_AR(data_train, data_test, species_list),
    AR_SINGLE_SPECIES = fit_AR_SINGLE_SPECIES(
      data_train,
      data_test,
      species_list
    ),
    GAM_AR = fit_GAM_AR(data_train, data_test, species_list),
    GAM_VAR = fit_GAM_VAR(data_train, data_test, species_list),
    SIMPLE = fit_SIMPLE(data_train, data_test, species_list),
    stop(glue::glue("No fit_model() branch for model '{model_name}'"))
  )
}

# Fail at startup rather than mid-run when a configured model has no
# R/models/<NAME>.R file defining fit_<NAME>().
check_model_definitions <- function(model_names) {
  undefined <- model_names[
    !purrr::map_lgl(paste0("fit_", model_names), exists, mode = "function")
  ]
  if (length(undefined) > 0) {
    stop(glue::glue(
      "No fit function for model(s) {paste(undefined, collapse = ', ')}. ",
      "Each needs R/models/<NAME>.R defining fit_<NAME>() ",
      "and a branch in fit_model()."
    ))
  }
  invisible(model_names)
}

# Score, summarize and LOO a fitted model, tagging each output with the window
# and species it came from so results can be combined across windows.
evaluate_model <- function(
  model,
  model_name,
  species_set,
  test_start,
  species_list
) {
  species_str <- paste(species_list, collapse = "_")
  rhats <- rhat(model)

  model_score <- score(forecast(model), score = "crps")
  model_score$test_start_newmoonnumber <- test_start
  model_score$species_set <- species_set
  model_score$species_list <- species_str
  model_score$rhat <- mean(rhats, na.rm = TRUE)
  model_score$prhat_high <- mean(rhats > 1.05, na.rm = TRUE)
  model_score$n_divergences <- sum(sapply(
    rstan::get_sampler_params(model$model_output, inc_warmup = FALSE),
    function(x) sum(x[, 'divergent__'])
  ))

  model_summary <- summary(model)
  model_summary$test_start_newmoonnumber <- test_start
  model_summary$species_set <- species_set
  model_summary$species_list <- species_str
  model_summary$model <- model_name

  model_loo <- loo(model)
  model_loo$test_start_newmoonnumber <- test_start
  model_loo$species_set <- species_set

  list(score = model_score, summary = model_summary, loo = model_loo)
}
