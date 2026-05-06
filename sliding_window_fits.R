library(dplyr)
library(cmdstanr)
library(forecast)
library(mvgam)
library(portalr)
library(glue)
library(statip)
library(furrr)
source("R/get_regime.R")

# Workers for the outer (train_start) loop. Each worker also spawns
# cmdstanr chains, so total cores in use ~= n_workers * 4.
n_workers <- as.integer(Sys.getenv("MVGAM_N_WORKERS", unset = "4"))
plan(multicore, workers = n_workers)

data_all <- readRDS("data_heteromyid.rds")

split_train_test <- function(data_all, gap, train_start, train_end, test_start, test_end) {
  species_list <- data.frame(newmoonnumber=data_all$newmoonnumber,series=data_all$series,y=data_all$y) |>
    filter(newmoonnumber >= train_start, newmoonnumber <= (train_end)) |>
    group_by(newmoonnumber, series) |>
    summarise(abundance = sum(y, na.rm = TRUE), .groups = 'drop') |>
    group_by(series) |>
    summarise(occupancy = sum(abundance > 0) / n(), .groups = 'drop') |>
    filter(occupancy >= 0.30) |>
    pull(series) |>
    droplevels()

  train_inds <- which(data_all$newmoonnumber >= train_start &
    data_all$newmoonnumber <= (train_end) &     # Old +1 here keeps the first observation after the gap to serve as the initial condition, also cut +gap since no gaps
    data_all$series %in% species_list)
  data_train <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][train_inds, ]
    } else {
      data_all[[x]][train_inds]
    }
  })

  names(data_train) <- names(data_all)
  data_train$series <- droplevels(data_train$series)
  # gap_inds <- which(data_train$newmoonnumber >= (train_end + 1) &
  #   data_train$newmoonnumber <= (train_end + gap))
  # data_train$y[gap_inds] <- NA

  test_inds <- which(data_all$newmoonnumber >= (test_start) & # Old +1 here was to skip the initial condition point
    data_all$newmoonnumber <= test_end &
    data_all$series %in% species_list)
  data_test <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][test_inds, ]
    } else {
      data_all[[x]][test_inds]
    }
  })

  names(data_test) <- names(data_all)
  data_test$series <- droplevels(data_test$series)

  return(list(train = data_train, test = data_test, species_list=species_list))
}

# Priors

sigma_prior <- prior(beta(10, 10), class = sigma, lb = 0.2, ub = 1)

ndvi_random_slopes_prior <- prior(
  inv_gamma(2.3693353, 0.7311319),
  class = sigma_raw_trend
)
# AR model prior in Clark et al. 2025 https://github.com/nicholasjclark/portal_VAR/blob/main/2.%20models.R
# ar_sp_intercept_prior <- prior(std_normal(), class = b)

ar_priors <- c(sigma_prior)
gam_ar_priors <- c(sigma_prior, ndvi_random_slopes_prior)
gam_var_priors <- c(sigma_prior, ndvi_random_slopes_prior)

trend_formula_PP = ~ s(ndvi_ma12, trend, bs = "re") +
  te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr")) +
  te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c("tp", "cr")) +
  te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c("tp", "cr")) +
  te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c("tp", "cr"))

trend_formula_DX = ~ s(ndvi_ma12, trend, bs = "re") +
  te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr")) +
  te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c("tp", "cr")) +
  te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c("tp", "cr"))

newmoon_min <- min(data_all$newmoonnumber)
newmoon_max <- max(data_all$newmoonnumber)
train_win_width <- 60
train_starts <- newmoon_min:(newmoon_max - train_win_width - 12 + 1)

# For non-full runs uncomment the lines below and specify desired
# test starts as newmoonnumbers.
# test_starts = seq(from = 200, to = 400, by = 20)
# train_starts = test_starts - train_win_width

run_window <- function(train_start) {
  train_end <- train_start + train_win_width - 1
  test_start <- train_end + 1
  test_end <- test_start + 12 - 1
  print(glue("Training test start {test_start} (newmoon)"))
  data_split <- split_train_test(
    data_all,
    gap = 0,
    train_start = train_start,
    train_end = train_end,
    test_start = test_start,
    test_end = test_end
  )
  data_train <- data_split$train
  data_test <- data_split$test

  if("PP" %in% levels(data_split$species_list)) {
    trend_formula = trend_formula_PP
  } else {
    trend_formula = trend_formula_DX
  }

  baseline_model <- mvgam(
    formula = y ~ -1 + series, #remove global intercept and allow species-specific intercepts
    data = data_train,
    newdata = data_test,
    family = poisson(),
    silent = 2,
    refresh = 0
    #remove prior specification, since it is only sigma and not used
    #solves Warning message: no match found in model_file for parameter: sigma
  )

  ar_model <- mvgam(
    formula = y ~ -1 + series,
    data = data_train,
    newdata = data_test,
    family = poisson(),
    trend_model = AR(),
    priors = ar_priors,
    burnin = 5000,
    samples = 2000,
    silent = 2,
    refresh = 0
  )

  gam_ar_model <- mvgam(
    formula = y ~ -1,
    trend_formula = trend_formula,
    data = data_train,
    newdata = data_test,
    family = poisson(),
    trend_model = AR(),
    priors = gam_ar_priors,
    burnin = 5000,
    samples = 2000,
    silent = 2,
    refresh = 0,
    control = list(adapt_delta = 0.95) # Increase from default 0.8 to decrease divergences
  )

  gam_var_model <- mvgam(
    formula = y ~ -1,
    trend_formula = trend_formula,
    data = data_train,
    newdata = data_test,
    family = poisson(),
    trend_model = VAR(),
    priors = gam_var_priors,
    burnin = 5000,
    samples = 2000,
    silent = 2,
    refresh = 0,
    control = list(adapt_delta = 0.95)
  )

  simple_model <- mvgam(
    formula = y ~ -1,
    trend_formula = ~ te(mintemp_lag_0, delta_mintemp, by = trend, k = 4, bs = "sz") +
      s(summer_ndvi, by = trend, k = 4) +
      s(winter_ndvi, by = trend, k = 4),
   data = data_train,
   newdata = data_test,
   family = poisson(),
   trend_model = AR(),
   noncentred = TRUE,
   burnin = 5000,
   samples = 2000,
   silent = 2,
   refresh = 0,
   control = list(adapt_delta = 0.95) # Increase from default 0.8 to decrease divergences
 )

  baseline_score <- score(forecast(baseline_model), score = "crps")
  baseline_score$test_start_newmoonnumber <- test_start
  baseline_score$species_list <- paste(data_split$species_list,collapse="_")
  baseline_score$rhat <- mean(rhat(baseline_model),na.rm=TRUE)
  baseline_score$prhat_high <- mean(rhat(baseline_model)>1.05,na.rm=TRUE)
  baseline_score$n_divergences <- sum(sapply(
    rstan::get_sampler_params(baseline_model$model_output, inc_warmup = FALSE),
    function(x) sum(x[, 'divergent__'])))
  baseline_summary <- summary(baseline_model)
  baseline_summary$test_start_newmoonnumber <- test_start
  baseline_summary$species_list <- paste(data_split$species_list,collapse="_")

  ar_score <- score(forecast(ar_model), score = "crps")
  ar_score$test_start_newmoonnumber <- test_start
  ar_score$species_list <- paste(data_split$species_list,collapse="_")
  ar_score$rhat <- mean(rhat(ar_model),na.rm=TRUE)
  ar_score$prhat_high <- mean(rhat(ar_model)>1.05,na.rm=TRUE)
  ar_score$n_divergences <- sum(sapply(
    rstan::get_sampler_params(ar_model$model_output, inc_warmup = FALSE),
    function(x) sum(x[, 'divergent__'])))
  ar_summary <- summary(ar_model)
  ar_summary$test_start_newmoonnumber <- test_start
  ar_summary$species_list <- paste(data_split$species_list,collapse="_")

  gam_ar_score <- score(forecast(gam_ar_model), score = "crps")
  gam_ar_score$test_start_newmoonnumber <- test_start
  gam_ar_score$species_list <- paste(data_split$species_list,collapse="_")
  gam_ar_score$rhat <- mean(rhat(gam_ar_model),na.rm=TRUE)
  gam_ar_score$prhat_high <- mean(rhat(gam_ar_model)>1.05,na.rm=TRUE)
  gam_ar_score$n_divergences <- sum(sapply(
    rstan::get_sampler_params(gam_ar_model$model_output, inc_warmup = FALSE),
    function(x) sum(x[, 'divergent__'])))
  gam_ar_summary <- summary(gam_ar_model)
  gam_ar_summary$test_start_newmoonnumber <- test_start
  gam_ar_summary$species_list <- paste(data_split$species_list,collapse="_")

  gam_var_score <- score(forecast(gam_var_model), score = "crps")
  gam_var_score$test_start_newmoonnumber <- test_start
  gam_var_score$species_list <- paste(data_split$species_list,collapse="_")
  gam_var_score$rhat <- mean(rhat(gam_var_model),na.rm=TRUE)
  gam_var_score$prhat_high <- mean(rhat(gam_var_model)>1.05,na.rm=TRUE)
  gam_var_score$n_divergences <- sum(sapply(
    rstan::get_sampler_params(gam_var_model$model_output, inc_warmup = FALSE),
    function(x) sum(x[, 'divergent__'])))
  gam_var_summary <- summary(gam_var_model)
  gam_var_summary$test_start_newmoonnumber <- test_start
  gam_var_summary$species_list <- paste(data_split$species_list,collapse="_")

  simple_score <- score(forecast(simple_model), score = "crps")
  simple_score$test_start_newmoonnumber <- test_start
  simple_score$species_list <- paste(data_split$species_list,collapse="_")
  simple_score$rhat <- mean(rhat(simple_model),na.rm=TRUE)
  simple_score$prhat_high <- mean(rhat(simple_model)>1.05,na.rm=TRUE)
  simple_score$n_divergences <- sum(sapply(
    rstan::get_sampler_params(simple_model$model_output, inc_warmup = FALSE),
    function(x) sum(x[, 'divergent__'])))
  simple_summary <- summary(simple_model)
  simple_summary$test_start_newmoonnumber <- test_start
  simple_summary$species_list <- paste(data_split$species_list,collapse="_")

  env_distance <- data.frame(ndvi=hellinger(data_train$ndvi_ma12, data_test$ndvi_ma12),
                             mintemp=hellinger(data_train$mintemp, data_test$mintemp))
  env_distance$test_start_newmoonnumber <- test_start
  env_distance$species_list <- paste(data_split$species_list,collapse="_")

  source("R/skill_scores.r", local = TRUE)
  source("R/forecast_plots.r", local = TRUE, echo = TRUE)

  list(
    baseline_score = baseline_score,
    ar_score = ar_score,
    gam_ar_score = gam_ar_score,
    gam_var_score = gam_var_score,
    simple_score = simple_score,
    baseline_summary = baseline_summary,
    ar_summary = ar_summary,
    gam_ar_summary = gam_ar_summary,
    gam_var_summary = gam_var_summary,
    simple_summary = simple_summary,
    env_distance = env_distance
  )
}

safe_run_window <- purrr::safely(run_window)
results <- future_map(
  train_starts,
  safe_run_window,
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

errored <- purrr::map_lgl(results, ~ !is.null(.x$error))
if (any(errored)) {
  warning(glue("{sum(errored)} window(s) failed: train_starts {paste(train_starts[errored], collapse=', ')}"))
}
results <- purrr::map(results, "result")

baseline_scores <- purrr::map(results, "baseline_score")
ar_scores <- purrr::map(results, "ar_score")
gam_ar_scores <- purrr::map(results, "gam_ar_score")
gam_var_scores <- purrr::map(results, "gam_var_score")
simple_scores <- purrr::map(results, "simple_score")

baseline_summaries <- purrr::map(results, "baseline_summary")
ar_summaries <- purrr::map(results, "ar_summary")
gam_ar_summaries <- purrr::map(results, "gam_ar_summary")
gam_var_summaries <- purrr::map(results, "gam_var_summary")
simple_summaries <- purrr::map(results, "simple_summary")
env_distances <- purrr::map(results, "env_distance")

saveRDS(baseline_scores, "baseline_scores.rds")
saveRDS(ar_scores, "ar_scores.rds")
saveRDS(gam_ar_scores, "gam_ar_scores.rds")
saveRDS(gam_var_scores, "gam_var_scores.rds")
saveRDS(simple_scores, "simple_scores.rds")
saveRDS(ar_summaries, "ar_summaries.rds")
saveRDS(gam_ar_summaries, "gam_ar_summaries.rds")
saveRDS(gam_var_summaries, "gam_var_summaries.rds")
saveRDS(simple_summaries, "simple_summaries.rds")
saveRDS(env_distances, "env_distances.rds")

# baseline_scores <- readRDS("baseline_scores.rds")
# ar_scores <- readRDS("ar_scores.rds")
# gam_ar_scores <- readRDS("gam_ar_scores.rds")
# gam_var_scores <- readRDS("gam_var_scores.rds")
# simple_scores <- readRDS("simple_scores.rds")
# ar_summaries <- readRDS("ar_summaries.rds")
# gam_ar_summaries <- readRDS("gam_ar_summaries.rds")
# gam_var_summaries <- readRDS("gam_var_summaries.rds")
# simple_summaries <- readRDS("simple_summaries.rds")
# env_distances <- readRDS("env_distances.rds")
