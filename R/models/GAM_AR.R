fit_GAM_AR <- function(data_train, data_test, species_list) {
  mvgam(
    formula = y ~ -1,
    trend_formula = get_trend_formula(species_list),
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
  )
}
