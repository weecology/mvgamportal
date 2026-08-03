fit_SIMPLE <- function(data_train, data_test, species_list) {
  mvgam(
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
  )
}
