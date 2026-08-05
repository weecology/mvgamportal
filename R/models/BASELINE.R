fit_BASELINE <- function(data_train, data_test, species_list) {
  #remove global intercept and allow species-specific intercepts, except with a
  #single species where `series` is a one-level factor that cannot be contrast
  #coded and the species-specific intercepts collapse to one intercept anyway
  baseline_formula <- if (length(species_list) > 1) {
    y ~ -1 + series
  } else {
    y ~ 1
  }

  mvgam(
    formula = baseline_formula,
    data = data_train,
    newdata = data_test,
    family = poisson(),
    silent = 2,
    refresh = 0
  )
}
