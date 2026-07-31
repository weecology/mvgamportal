fit_BASELINE <- function(data_train, data_test, species_list) {
  mvgam(
    formula = y ~ -1 + series, #remove global intercept and allow species-specific intercepts
    data = data_train,
    newdata = data_test,
    family = poisson(),
    silent = 2,
    refresh = 0
  )
}
