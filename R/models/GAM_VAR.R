fit_GAM_VAR <- function(data_train, data_test, species_list) {
  mvgam(
    formula = y ~ -1,
    trend_formula = get_trend_formula(species_list),
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
  )
}
