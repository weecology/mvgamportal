library(dplyr)
library(cmdstanr)
library(forecast)
library(mvgam)
library(portalr)
library(glue)
library(statip)
source("R/get_regime.R")

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

baseline_scores <- vector(mode = "list", length = length(train_starts))
ar_scores <- vector(mode = "list", length = length(train_starts))
gam_ar_scores <- vector(mode = "list", length = length(train_starts))
gam_var_scores <- vector(mode = "list", length = length(train_starts))

baseline_summaries <- vector(mode = "list", length = length(train_starts))
ar_summaries <- vector(mode = "list", length = length(train_starts))
gam_ar_summaries <- vector(mode = "list", length = length(train_starts))
gam_var_summaries <- vector(mode = "list", length = length(train_starts))
env_distances <- vector(mode = "list", length = length(train_starts))

targets = c(158,180,200,220,240,248,260,280,285,300,320,340,360,380,400,420,428,440,460,480)
initial = targets - 60

for (i in initial) { #seq_along(train_starts)
  target = i + 60
  train_start <- train_starts[i - 95] # remove the -95
  print(glue("Starting training: {i}"))
  train_end <- train_start + 60 - 1
  test_start <- train_end + 1
  test_end <- test_start + 12 - 1
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
    family = poisson()         
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
    samples = 2000
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
    control = list(adapt_delta = 0.95)
  )

  baseline_score <- score(forecast(baseline_model), score = "crps")
  baseline_score$test_start_newmoonnumber <- test_start
  baseline_score$species_list <- paste(data_split$species_list,collapse="_")
  baseline_score$rhat <- mean(rhat(baseline_model),na.rm=TRUE)
  baseline_score$prhat_high <- mean(rhat(baseline_model)>1.05,na.rm=TRUE)
  baseline_score$n_divergences <- sum(sapply(
    rstan::get_sampler_params(baseline_model$model_output, inc_warmup = FALSE), 
    function(x) sum(x[, 'divergent__'])))
  baseline_scores[[i]] <- baseline_score
  baseline_summary <- summary(baseline_model)
  baseline_summary$test_start_newmoonnumber <- test_start
  baseline_summary$species_list <- paste(data_split$species_list,collapse="_")
  baseline_summaries[[i]] <- baseline_summary
  
  ar_score <- score(forecast(ar_model), score = "crps")
  ar_score$test_start_newmoonnumber <- test_start
  ar_score$species_list <- paste(data_split$species_list,collapse="_")
  ar_score$rhat <- mean(rhat(ar_model),na.rm=TRUE)
  ar_score$prhat_high <- mean(rhat(ar_model)>1.05,na.rm=TRUE)
  ar_score$n_divergences <- sum(sapply(
    rstan::get_sampler_params(ar_model$model_output, inc_warmup = FALSE), 
    function(x) sum(x[, 'divergent__'])))
  ar_scores[[i]] <- ar_score
  ar_summary <- summary(ar_model)
  ar_summary$test_start_newmoonnumber <- test_start
  ar_summary$species_list <- paste(data_split$species_list,collapse="_")
  ar_summaries[[i]] <- ar_summary
  
  gam_ar_score <- score(forecast(gam_ar_model), score = "crps")
  gam_ar_score$test_start_newmoonnumber <- test_start
  gam_ar_score$species_list <- paste(data_split$species_list,collapse="_")
  gam_ar_score$rhat <- mean(rhat(gam_ar_model),na.rm=TRUE)
  gam_ar_score$prhat_high <- mean(rhat(gam_ar_model)>1.05,na.rm=TRUE)
  gam_ar_score$n_divergences <- sum(sapply(
    rstan::get_sampler_params(gam_ar_model$model_output, inc_warmup = FALSE), 
    function(x) sum(x[, 'divergent__'])))
  gam_ar_scores[[i]] <- gam_ar_score
  gam_ar_summary <- summary(gam_ar_model)
  gam_ar_summary$test_start_newmoonnumber <- test_start
  gam_ar_summary$species_list <- paste(data_split$species_list,collapse="_")
  gam_ar_summaries[[i]] <- gam_ar_summary
  
  gam_var_score <- score(forecast(gam_var_model), score = "crps")
  gam_var_score$test_start_newmoonnumber <- test_start
  gam_var_score$species_list <- paste(data_split$species_list,collapse="_")
  gam_var_score$rhat <- mean(rhat(gam_var_model),na.rm=TRUE)
  gam_var_score$prhat_high <- mean(rhat(gam_var_model)>1.05,na.rm=TRUE)
  gam_var_score$n_divergences <- sum(sapply(
    rstan::get_sampler_params(gam_var_model$model_output, inc_warmup = FALSE), 
    function(x) sum(x[, 'divergent__'])))
  gam_var_scores[[i]] <- gam_var_score
  gam_var_summary <- summary(gam_var_model)
  gam_var_summary$test_start_newmoonnumber <- test_start
  gam_var_summary$species_list <- paste(data_split$species_list,collapse="_")
  gam_var_summaries[[i]] <- gam_var_summary
  
  env_distance <- data.frame(ndvi=hellinger(data_train$ndvi_ma12, data_test$ndvi_ma12),
                             mintemp=hellinger(data_train$mintemp, data_test$mintemp))
  env_distance$test_start_newmoonnumber <- test_start
  env_distance$species_list <- paste(data_split$species_list,collapse="_")
  env_distances[[i]] <- env_distance

  source("./R/skill_scores.r")
  source("./R/forecast_plots.r", echo = TRUE)
}
