library(dplyr)
library(lubridate)
library(portalr)
library(tidyr)

#' Create distributed lag matrices for environmental covariates
#' Code modified from Clark et al. 2024
#' df should have a newmoonnumber column and a value column
#' the value column will be lagged
lagard <- function(df, n_lag = 6) {
  covar <- names(df)[2]
  n <- nrow(df)
  time_series <- covars[[2]] # Get just values
  lags <- matrix(NA, n, n_lag + 1)
  for (i in 1:(n_lag + 1)) {
    lags[i:n, i] <- time_series[i:n - i + 1]
  }
  lags_df <- data.frame(lags)
  colnames(lags_df) <- paste0(covar, "_lag_", 0:(n_lag))
  bind_cols(df["newmoonnumber"], lags_df)
}

download_observations()

rodent_data <- summarize_rodent_data(
  path = get_default_data_path(),
  clean = FALSE,
  level = "Treatment",
  type = "Rodents",
  plots = "Longterm",
  unknowns = FALSE,
  shape = "long",
  time = "all",
  output = "abundance",
  na_drop = FALSE,
  zero_drop = FALSE,
  min_traps = 1,
  min_plots = 1,
  effort = TRUE,
  download_if_missing = TRUE,
  quiet = FALSE
)

target_species <- c("DM", "PP", "OT", "PE", "PB", "RM", "DO")

# In comparison to Clarke et al. 2025
# we've added the correction for under sampling of controls
rodent_data <- rodent_data |>
  filter(treatment == "control") |>
  filter(species %in% target_species) |>
  mutate(species = factor(species, levels = unique(species))) |>
  mutate(abundance = as.integer(round(abundance * 4 / nplots, 0))) |>
  select(-ntraps, -nplots)
max_newmoon <- max(rodent_data$newmoonnumber)

weather <- weather(
  level = "newmoon",
  fill = TRUE,
  horizon = 365
)

# NDVI not filled in Clarke et al. 2025
ndvi <- ndvi(level = "newmoon", fill = TRUE) |>
  filter(newmoonnumber <= max_newmoon)
min_newmoon <- min(ndvi$newmoonnumber)

# Different from Clarke et al. 2025 we are scaling the ma12 of ndvi not
# the raw ndvi which should give us better values for back transforming
covars <- weather |>
  left_join(ndvi) |>
  filter(newmoonnumber >= min_newmoon, newmoonnumber <= max_newmoon) |>
  mutate(
    ndvi_ma12 = scale(zoo::rollmean(
      ndvi,
      k = 12,
      align = "right",
      fill = NA
    )),
    mintemp_ma3 = scale(zoo::rollmean(
      mintemp,
      k = 3,
      align = "right",
      na.pad = TRUE
    )),
    maxtemp_ma3 = scale(zoo::rollmean(
      maxtemp,
      k = 3,
      align = "right",
      na.pad = TRUE
    )),
    mintemp = scale(mintemp),
    maxtemp = scale(maxtemp),
    # From Pat's paper
    meantemp_lag_1 = scale(lag(meantemp, order_by = newmoonnumber)),
  ) |>
  select(
    newmoonnumber,
    meantemp,
    mintemp,
    maxtemp,
    precipitation,
    warm_precip,
    cool_precip,
    ndvi,
    ndvi_ma12,
    meantemp_lag_1,
    mintemp_ma3,
    maxtemp_ma3
  )

# Temperature 6-month lag matrices
mintemp_df <- lagard(select(covars, newmoonnumber, mintemp), 6)
maxtemp_df <- lagard(select(covars, newmoonnumber, maxtemp), 6)

covars_w_lags <- covars |>
  inner_join(mintemp_df, by = c("newmoonnumber")) |>
  inner_join(maxtemp_df, by = c("newmoonnumber"))

model_dat <- rodent_data |>
  right_join(covars_w_lags, by = "newmoonnumber") |>
  drop_na(meantemp:maxtemp_lag_5) |>
  # add mvgam columns
  mutate(
    y = abundance,
    series = species,
    time = newmoonnumber - (min(newmoonnumber) - 1)
  ) |>
  select(time, y, series, newmoonnumber, meantemp:maxtemp_lag_6) |>
  arrange(time, series)

# Compared to Clarke et al. 2025 we don't need to filter species
# by occurance because we already filtered by name

lag <- matrix(0:5, nrow(model_dat), 6, byrow = TRUE)

# Create weight matrices for hierarchical distributed lag terms
# Directly from Clarke et al. 2025's code
weights_dm <- weights_do <-
  weights_pp <- weights_ol <-
  weights_ot <- weights_pf <-
  weights_pb <- weights_pe <- weights_rm <-
  matrix(1, ncol = ncol(lag), nrow = nrow(lag))

weights_dm[!(model_dat$series == "DM"), ] <- 0
weights_do[!(model_dat$series == "DO"), ] <- 0
weights_ol[!(model_dat$series == "OL"), ] <- 0
weights_ot[!(model_dat$series == "OT"), ] <- 0
weights_pb[!(model_dat$series == "PB"), ] <- 0
weights_pe[!(model_dat$series == "PE"), ] <- 0
weights_pf[!(model_dat$series == "PF"), ] <- 0
weights_pp[!(model_dat$series == "PP"), ] <- 0
weights_rm[!(model_dat$series == "RM"), ] <- 0

data_all <- list(
  lag = lag,
  meantemp = model_dat$meantemp,
  meantemp_lag_1 = model_dat$meantemp_lag_1,
  mintemp = as.matrix(select(model_dat, mintemp_lag_1:mintemp_lag_6)),
  mintemp_ma3 = model_dat$mintemp_ma3,
  maxtemp = as.matrix(select(model_dat, maxtemp_lag_1:maxtemp_lag_6)),
  maxtemp_ma3 = model_dat$maxtemp_ma3,
  warm_precip = model_dat$warm_precip,
  cool_precip = model_dat$cool_precip,
  ndvi_ma12 = model_dat$ndvi_ma12,
  weights_dm = weights_dm,
  weights_do = weights_do,
  weights_ol = weights_ol,
  weights_ot = weights_ot,
  weights_pb = weights_pb,
  weights_pe = weights_pe,
  weights_pf = weights_pf,
  weights_pp = weights_pp,
  weights_rm = weights_rm,
  y = model_dat$y,
  series = model_dat$series,
  time = model_dat$time
)

saveRDS(data_all, file = "data_all.rds")
