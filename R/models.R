library(mvgam)

# Model definitions and output handling for the sliding window analysis.
# Which models are run is set in config.yaml.

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

trend_formula_PP <- ~ s(ndvi_ma12, trend, bs = "re") +
  te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr")) +
  te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c("tp", "cr")) +
  te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c("tp", "cr")) +
  te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c("tp", "cr"))

trend_formula_DX <- ~ s(ndvi_ma12, trend, bs = "re") +
  te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr")) +
  te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c("tp", "cr")) +
  te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c("tp", "cr"))

get_trend_formula <- function(species_list) {
  if ("PP" %in% levels(species_list)) {
    trend_formula_PP
  } else {
    trend_formula_DX
  }
}

fit_model <- function(model_name, data_train, data_test, species_list) {
  trend_formula <- get_trend_formula(species_list)

  switch(
    model_name,
    BASELINE = mvgam(
      formula = y ~ -1 + series, #remove global intercept and allow species-specific intercepts
      data = data_train,
      newdata = data_test,
      family = poisson(),
      silent = 2,
      refresh = 0
    ),
    AR = mvgam(
      formula = y ~ -1 + series,
      data = data_train,
      newdata = data_test,
      family = nb(),
      trend_model = AR(),
      priors = ar_priors,
      burnin = 5000,
      samples = 2000,
      silent = 2,
      refresh = 0
    ),
    GAM_AR = mvgam(
      formula = y ~ -1,
      trend_formula = trend_formula,
      data = data_train,
      newdata = data_test,
      family = nb(),
      trend_model = AR(),
      priors = gam_ar_priors,
      burnin = 5000,
      samples = 2000,
      silent = 2,
      refresh = 0,
      control = list(adapt_delta = 0.95) # Increase from default 0.8 to decrease divergences
    ),
    GAM_VAR = mvgam(
      formula = y ~ -1,
      trend_formula = trend_formula,
      data = data_train,
      newdata = data_test,
      family = nb(),
      trend_model = VAR(),
      priors = gam_var_priors,
      burnin = 5000,
      samples = 2000,
      silent = 2,
      refresh = 0,
      control = list(adapt_delta = 0.95)
    ),
    SIMPLE = mvgam(
      formula = y ~ -1,
      trend_formula = ~ te(
        mintemp_lag_0,
        delta_mintemp,
        by = trend,
        k = 4,
        bs = "sz"
      ) +
        s(summer_ndvi, by = trend, k = 4) +
        s(winter_ndvi, by = trend, k = 4),
      data = data_train,
      newdata = data_test,
      family = nb(),
      trend_model = AR(),
      noncentred = TRUE,
      burnin = 5000,
      samples = 2000,
      silent = 2,
      refresh = 0,
      control = list(adapt_delta = 0.95) # Increase from default 0.8 to decrease divergences
    ),
    stop(glue::glue("No model definition for '{model_name}' in fit_model()"))
  )
}

# Score, summarize and LOO a fitted model, tagging each output with the window
# and species it came from so results can be combined across windows.
evaluate_model <- function(model, model_name, test_start, species_list) {
  species_str <- paste(species_list, collapse = "_")
  rhats <- rhat(model)

  model_score <- score(forecast(model), score = "crps")
  model_score$test_start_newmoonnumber <- test_start
  model_score$species_list <- species_str
  model_score$rhat <- mean(rhats, na.rm = TRUE)
  model_score$prhat_high <- mean(rhats > 1.05, na.rm = TRUE)
  model_score$n_divergences <- sum(sapply(
    rstan::get_sampler_params(model$model_output, inc_warmup = FALSE),
    function(x) sum(x[, 'divergent__'])
  ))

  model_summary <- summary(model)
  model_summary$test_start_newmoonnumber <- test_start
  model_summary$species_list <- species_str
  model_summary$model <- model_name

  model_loo <- loo(model)
  model_loo$test_start_newmoonnumber <- test_start

  list(score = model_score, summary = model_summary, loo = model_loo)
}
