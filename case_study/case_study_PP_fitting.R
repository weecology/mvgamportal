path <- "./case_study/"
data_all <- readRDS(paste0(path, "data_pb_regime.rds"))

species <- "PP"
pp_rows <- which(data_all$series == species)   # Identify rows for PP
species_index <- which(levels(data_all$series) == species)
species_index
# Subset predictors and weights for PP

lag_pp     <- data_all$lag[pp_rows, , drop = FALSE]
mintemp_pp <- data_all$mintemp[pp_rows, , drop = FALSE]
maxtemp_pp <- data_all$maxtemp[pp_rows, , drop = FALSE]
weights_pp <- data_all$weights_pp[pp_rows, , drop = FALSE]

mintemp_ma3 <- data_all$mintemp_ma3[pp_rows, , drop = FALSE]
maxtemp_ma3 <- data_all$maxtemp_ma3[pp_rows, , drop = FALSE]
ndvi_ma12   <- data_all$ndvi_ma12[pp_rows, , drop = FALSE]

time_pp          <- data_all$time[pp_rows]
newmoonnumber_pp <- data_all$newmoonnumber[pp_rows]
meantemp_pp      <- data_all$meantemp[pp_rows]
meantemp_lag1_pp <- data_all$meantemp_lag_1[pp_rows]

data_pp <- list(
  y = data_all$y[pp_rows],
  series = factor(rep(species, length(pp_rows)), levels = species),
  lag = lag_pp,
  mintemp = mintemp_pp,
  maxtemp = maxtemp_pp,
  mintemp_ma3 = mintemp_ma3,
  maxtemp_ma3 = maxtemp_ma3,
  ndvi_ma12 = ndvi_ma12,
  weights_pp = weights_pp,
  time = time_pp,
  newmoonnumber = newmoonnumber_pp,
  meantemp = meantemp_pp,
  meantemp_lag_1 = meantemp_lag1_pp
)

data_pp$lag

stopifnot(
  length(data_pp$y) == nrow(data_pp$lag),
  length(data_pp$y) == nrow(data_pp$mintemp),
  length(data_pp$y) == nrow(data_pp$maxtemp),
  length(data_pp$y) == nrow(data_pp$weights_pp),
  length(levels(data_pp$series)) == 1
)
message("Data cleaned for PP: rows = ", length(data_pp$y),
        ", series levels = ", paste(levels(data_pp$series), collapse = ", "))



split_train_test <- function(data_pp, gap, train_start, train_end, test_start, test_end) {
  train_inds <- which(data_pp$newmoonnumber >= train_start &
                        data_pp$newmoonnumber <= (train_end + gap + 1))
  data_train <- lapply(seq_along(data_pp), function(x) {
    if (is.matrix(data_pp[[x]])) {
      data_pp[[x]][train_inds, , drop = FALSE]
    } else {
      data_pp[[x]][train_inds]
    }
  })
  names(data_train) <- names(data_pp)
  

  if("series" %in% names(data_train)) {
    data_train$series <- droplevels(data_train$series)
  }
  
  gap_inds <- which(data_train$newmoonnumber >= (train_end + 1) &
                      data_train$newmoonnumber <= (train_end + gap))
  data_train$y[gap_inds] <- NA
  
  test_inds <- which(data_pp$newmoonnumber >= (test_start + 1) &
                       data_pp$newmoonnumber <= test_end)
  data_test <- lapply(seq_along(data_pp), function(x) {
    if (is.matrix(data_pp[[x]])) {
      data_pp[[x]][test_inds, , drop = FALSE]
    } else {
      data_pp[[x]][test_inds]
    }
  })
  names(data_test) <- names(data_pp)
  

  if("series" %in% names(data_test)) {
    data_test$series <- droplevels(data_test$series)
  }
  
  return(list(train = data_train, test = data_test))
}


fit_models_single_species <- function(data_pp, train_starts, type) {
  train_win_width <- 60
  gap <- 16
  horizon <- 7
  
  for (i in seq_along(train_starts)) {
    train_start <- train_starts[i]
    train_stop  <- train_start + train_win_width - 1
    test_start  <- train_stop + gap + 1
    test_stop   <- test_start + horizon - 1
    
    # Split data
    data_split <- split_train_test(data_pp, gap, train_start, train_stop, test_start, test_stop)
    data_train <- data_split$train
    data_test  <- data_split$test
    
    # Priors
    sigma_prior <- prior(beta(10, 10), class = sigma, lb = 0.2, ub = 1)
    intercept_prior <- prior(normal(0, 0.001), class = Intercept)
    ndvi_random_slopes_prior <- prior(inv_gamma(2.3693353, 0.7311319), class = sigma_raw_trend)
    ar_sp_intercept_prior <- prior(std_normal(), class = b)
    
    ar_priors <- c(intercept_prior, ar_sp_intercept_prior)
    gam_ar_priors <- c(intercept_prior, ndvi_random_slopes_prior)
    
    
    print(paste("start baseline model",i, sep=": "))
    baseline_model <- mvgam(
      formula = y ~ 1,
      data = data_train,
      newdata = data_test,
      family = poisson(),
      trend_model = "None",
      #priors = ar_priors,
      samples = 1600
    )
    
    linear_mintemp = mvgam(
      y ~ mintemp[,1],
      trend_model = AR(p = 1),
      family = poisson(),
      data = data_train,
      newdata = data_test
    )
    # print(paste("start gam_ar",i, sep=": "))
    # model_gam_ar <- mvgam(
    #   formula = y ~ -1,
    #   trend_formula = ~ s(ndvi_ma12, trend, bs = "re") +
    #     te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr")),
    #   data = data_train,
    #   newdata = data_test,
    #   family = poisson(),
    #   trend_model = AR(),
    #   priors = gam_ar_priors,
    #   samples = 1600
    # )
    # print(paste("start gam_var", i, sep=": "))
    # model_gam_var <- mvgam(
    #   formula = y ~ -1,
    #   trend_formula = ~ s(ndvi_ma12, trend, bs = "re") +
    #     te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr")),
    #   data = data_train,
    #   newdata = data_test,
    #   family = poisson(),
    #   trend_model = "VAR1",
    #   priors = gam_var_priors,
    #   samples = 1600
    # )
    # print(paste("start ar_model", i, sep=": "))
    # model_ar <- mvgam(
    #   formula = y ~ series,
    #   data = data_train,
    #   newdata = data_test,
    #   family = poisson(),
    #   trend_model = AR(),
    #   priors = ar_priors,
    #   samples = 1600
    # )
# 
#     saveRDS(model_gam_ar,
#             paste0(path, "model_outputs/gam_ar_", type, "_PP_regime_output_", test_start, ".rds"))
#     saveRDS(model_ar,
#             paste0(path, "model_outputs/ar_", type, "_PP_regime_output_", test_start, ".rds"))
    saveRDS(baseline_model,
            paste0(path, "model_outputs/baseline_", type, "regime_output", test_start, ".rds"))
    saveRDS(linear_mintemp,
            paste0(path, "model_outputs/linear_mintemp_", type, "regime_output", test_start, ".rds"))
    
  }
}


in_train_starts <- c(287, 299, 311)
out_train_start <- c(336)

fit_models_single_species(data_pp, in_train_starts, "in")
fit_models_single_species(data_pp, out_train_start, "out")

linear_mintemp <- readRDS("case_study/model_outputs/linear_mintemp_inregime_output387.rds")
baseline_model <-readRDS("case_study/model_outputs/baseline_inregime_output387.rds")
mcmc_plot(linear_mintemp, type = 'trace')

plot(linear_mintemp, type = 'forecast')
plot(baseline_model, type = 'forecast')



# From class
pp_data = read.csv("pp_abundance_timeseries.csv") |>
  mutate(time = newmoonnumber) |>
  mutate(series = as.factor("PP"))

data_train = filter(pp_data, time <=max(time) - 24) 
data_test = filter(pp_data, time > max(time) -24)

plot_mvgam_series(data = pp_data, y = "abundance")

baseline_model = mvgam(
  abundance ~ mintemp,
  trend_model = AR(p = 1),
  family = gaussian(),
  data = data_train
)

