library(dplyr)
library(lubridate)
library(portalr)
library(tidyr)

#' Create distributed lag matrices for environmental covariates
#' Code modified from Clark et al. 2024
#' df should have a time column, a series column, and a value column
#' the value column will be lagged
lagard <- function(df, n_lag = 6) {
  n <- nrow(df)
  time_series <- pull(select(df, -time, -series)) # Get just values
  lags <- matrix(NA, n, n_lag)
  for (i in 1:n_lag) {
    lags[i:n, i] <- time_series[i:n - i + 1]
  }
  lags_df <- data.frame(lags)
  bind_cols(select(df, time, series), lags_df)
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
    meantemp_lag1 = scale(lag(meantemp, order_by = newmoonnumber)), # From Pat's paper
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
    meantemp_lag1,
    mintemp_ma3,
    maxtemp_ma3
  )

# TODO: See if we need to keep the month and year columns created here
#       Think we don't need them but will keep for the moment
model_dat <- rodent_data |>
  right_join(covars, by = "newmoonnumber") |>
  # add mvgam columns + month and year
  mutate(
    y = abundance,
    series = species,
    time = newmoonnumber - (min(newmoonnumber) - 1),
    month = lubridate::month(censusdate),
    year = lubridate::year(censusdate)
  ) |>
  select(time, y, series, month, year, newmoonnumber, mintemp:maxtemp_ma3) |>
  drop_na(mintemp:maxtemp_ma3) |>
  arrange(time, series)

# Compared to Clarke et al. 2025 we don't need to filter species
# by occurance because we already filtered by name

# Mintemp 6-month lag matrix
mintemp_df <- lagard(select(model_dat, time, series, mintemp), 6)
maxtemp_df <- lagard(select(model_dat, time, series, maxtemp), 6)

lag <- matrix(0:5, nrow(model_dat),
  6,
  byrow = TRUE
)
