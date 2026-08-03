fit_AR <- function(data_train, data_test, species_list) {
  mvgam(
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
  )
}
