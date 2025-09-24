library(dplyr)
library(cmdstanr)
library(forecast)
library(mvgam)
library(portalr)
source("R/get_regime.R")

data_all <- readRDS("data_all.rds")

split_train_test <- function(data_all, gap, train_start, train_end, test_start, test_end) {
  train_inds <- which(data_all$newmoonnumber >= train_start &
    data_all$newmoonnumber <= (train_end)) # Old +1 here keeps the first observation after the gap to serve as the initial condition, also cut +gap since no gaps
  data_train <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][train_inds, ]
    } else {
      data_all[[x]][train_inds]
    }
  })

  names(data_train) <- names(data_all)
  # gap_inds <- which(data_train$newmoonnumber >= (train_end + 1) &
  #   data_train$newmoonnumber <= (train_end + gap))
  # data_train$y[gap_inds] <- NA

  test_inds <- which(data_all$newmoonnumber >= (test_start) & # Old +1 here was to skip the initial condition point
    data_all$newmoonnumber <= test_end)
  data_test <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][test_inds, ]
    } else {
      data_all[[x]][test_inds]
    }
  })

  names(data_test) <- names(data_all)

  return(list(train = data_train, test = data_test))
}

# Priors

sigma_prior <- prior(beta(10, 10), class = sigma, lb = 0.2, ub = 1)
# The observation-level intercept is a nuisance parameter in this model that we
# don't really need (and cannot adequately estimate anyway).
# At present there is no way to drop this term using mvgam, but we can
# heavily regularize it to zero with a strong prior
intercept_prior <- prior(normal(0, 0.001), class = Intercept)
ndvi_random_slopes_prior <- prior(
  inv_gamma(2.3693353, 0.7311319),
  class = sigma_raw_trend
)
# AR model prior in Clark et al. 2025 https://github.com/nicholasjclark/portal_VAR/blob/main/2.%20models.R
ar_sp_intercept_prior <- prior(std_normal(), class = b)

ar_priors <- c(sigma_prior, intercept_prior, ar_sp_intercept_prior)
gam_ar_priors <- c(sigma_prior, intercept_prior, ndvi_random_slopes_prior)
gam_var_priors <- c(sigma_prior, intercept_prior, ndvi_random_slopes_prior)

newmoon_min <- min(data_all$newmoonnumber)
newmoon_max <- max(data_all$newmoonnumber)
train_win_width <- 60
train_starts <- newmoon_min:97 # (newmoon_max - train_win_width - 12 + 1)

ar_scores <- vector(mode = "list", length = length(train_starts))
gam_ar_scores <- vector(mode = "list", length = length(train_starts))
gam_var_scores <- vector(mode = "list", length = length(train_starts))

ar_summaries <- vector(mode = "list", length = length(train_starts))
gam_ar_summaries <- vector(mode = "list", length = length(train_starts))
gam_var_summaries <- vector(mode = "list", length = length(train_starts))

for (i in seq_along(train_starts)) {
  train_start <- train_starts[i]
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

  ar_model <- mvgam(
    formula = y ~ series,
    data = data_train,
    newdata = data_test,
    family = poisson(),
    trend_model = AR(),
    priors = ar_priors,
    samples = 1600
  )

  gam_ar_model <- mvgam(
    formula = y ~ -1,
    trend_formula = ~ s(ndvi_ma12, trend, bs = "re") +
      te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr")) +
      te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c("tp", "cr")) +
      te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c("tp", "cr")) +
      te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c("tp", "cr")),
    data = data_train,
    newdata = data_test,
    family = poisson(),
    trend_model = AR(),
    priors = gam_ar_priors,
    samples = 1600
  )

  gam_var_model <- mvgam(
    formula = y ~ -1,
    trend_formula = ~ s(ndvi_ma12, trend, bs = "re") +
      te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr")) +
      te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c("tp", "cr")) +
      te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c("tp", "cr")) +
      te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c("tp", "cr")),
    data = data_train,
    newdata = data_test,
    family = poisson(),
    trend_model = "VAR1",
    priors = gam_var_priors,
    samples = 1600
  )
  # TODO
  # 3. Store convergence and other model fit info
  # 3. Make cool graphs
  # 4. ?
  # 5. Profit

  ar_score <- score(forecast(ar_model), score = "crps")
  ar_score$test_start_newmoonnumber <- test_start
  ar_scores[[i]] <- ar_score
  ar_summary <- summary(ar_model)
  ar_summary$test_start_newmoonnumber <- test_start
  ar_summaries[[i]] <- ar_summary
  gam_ar_score <- score(forecast(gam_ar_model), score = "crps")
  gam_ar_score$test_start_newmoonnumber <- test_start
  gam_ar_scores[[i]] <- gam_ar_score
  gam_ar_summary <- summary(gam_ar_model)
  gam_ar_summary$test_start_newmoonnumber <- test_start
  gam_ar_summaries[[i]] <- gam_ar_summary
  gam_var_score <- score(forecast(gam_var_model), score = "crps")
  gam_var_score$test_start_newmoonnumber <- test_start
  gam_var_scores[[i]] <- gam_var_score
  gam_var_summary <- summary(gam_var_model)
  gam_var_summary$test_start_newmoonnumber <- test_start
  gam_var_summaries[[i]] <- gam_var_summary
}
