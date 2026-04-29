library(dplyr)
library(lubridate)
library(portalr)
library(tidyr)
library(slider)

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
  type = "Granivores",
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

# Keep species present in >=29% of sampling periods
target_species <- rodent_data |>
  select(-c(censusdate, period, nplots, ntraps)) |>
  group_by(newmoonnumber, species) |>
  summarise(abundance = sum(abundance, na.rm = TRUE)) |>
  group_by(species) |>
  summarise(occupancy = sum(abundance > 0) / max(newmoonnumber)) |>
  filter(occupancy >= 0.29) |>
  select(species)

# In comparison to Clarke et al. 2025
# we've added the correction for under sampling of controls
rodent_data <- rodent_data |>
  filter(treatment == "control") |>
  filter(species %in% target_species$species) |>
  mutate(species = factor(species, levels = unique(species))) |>
  mutate(abundance = as.integer(round(abundance * 4 / nplots, 0))) |>
  select(-ntraps, -nplots)
max_newmoon <- max(rodent_data$newmoonnumber)

moon_dates <- load_datafile("Rodents/moon_dates.csv", na.strings = c("NA")) |>
  mutate(censusdate = as.Date(coalesce(censusdate, newmoondate)))

weather <- weather(
  level = "daily",
  fill = TRUE,
  horizon = 365) |> 
  left_join(moon_dates, by = c("date"="censusdate")) |>
 mutate(
  mintemp = slide_index_dbl(
   .x = mintemp,
   .i = date,
   .f = ~mean(.x, na.rm = TRUE),
   .before = days(28)
   ),
  maxtemp = slide_index_dbl(
   .x = maxtemp,
   .i = date,
   .f = ~mean(.x, na.rm = TRUE),
   .before = days(28)
   ),
  meantemp = slide_index_dbl(
   .x = meantemp,
   .i = date,
   .f = ~mean(.x, na.rm = TRUE),
   .before = days(28)
   ),
  precipitation = slide_index_dbl(
   .x = precipitation,
   .i = date,
   .f = ~sum(.x, na.rm = TRUE),
   .before = days(28)
 )) |>
 filter(!is.na(newmoonnumber)) |>
 mutate(
  mintemp_lag_0 = mintemp,
  mintemp_lag_1 = lag(mintemp),
  delta_mintemp = mintemp_lag_0 - mintemp_lag_1
 ) |>
 select(-c(mintemp_lag_0, mintemp_lag_1))

# NDVI not filled in Clarke et al. 2025
ndvi <- ndvi(level = "newmoon", fill = TRUE) |>
  filter(newmoonnumber <= max_newmoon) |>
  left_join(moon_dates, by = join_by(newmoonnumber)) |>
  mutate(year=lubridate::year(newmoondate),month=lubridate::month(newmoondate))
min_newmoon <- min(ndvi$newmoonnumber)

# Add seasonal NDVI
ndvi_seasons <- ndvi |>
  group_by(year) |>
  summarise(winter_ndvi=mean(ndvi[month %in% 3:4]),
            summer_ndvi=mean(ndvi[month %in% 7:9])) 

ndvi <- ndvi |>
  mutate(winter_ndvi = ifelse(
    month < 5, 
    ndvi_seasons$winter_ndvi[match(year-1, ndvi_seasons$year)], 
    ndvi_seasons$winter_ndvi[match(year, ndvi_seasons$year)]
  ),
  summer_ndvi = ifelse(
    month < 10, 
    ndvi_seasons$summer_ndvi[match(year-1, ndvi_seasons$year)], 
    ndvi_seasons$summer_ndvi[match(year, ndvi_seasons$year)]
  ))

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
    delta_mintemp = scale(delta_mintemp),
    maxtemp = scale(maxtemp),
    # From Pat's paper
    meantemp_lag_1 = scale(lag(meantemp, order_by = newmoonnumber)),
    ndvi = scale(ndvi),
    winter_ndvi = scale(winter_ndvi),
    summer_ndvi = scale(summer_ndvi)
  ) |>
  select(
    newmoonnumber,
    meantemp,
    mintemp,
    delta_mintemp,
    maxtemp,
    precipitation,
    warm_precip,
    cool_precip,
    ndvi,
    ndvi_ma12,
    winter_ndvi,
    summer_ndvi,
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
  weights_ds <- weights_pp <-
  weights_pf <- weights_pb <-
  weights_pe <- weights_rm <-
  matrix(1, ncol = ncol(lag), nrow = nrow(lag))

weights_dm[!(model_dat$series == "DM"), ] <- 0
weights_do[!(model_dat$series == "DO"), ] <- 0
weights_ds[!(model_dat$series == "DS"), ] <- 0
weights_pb[!(model_dat$series == "PB"), ] <- 0
weights_pe[!(model_dat$series == "PE"), ] <- 0
weights_pf[!(model_dat$series == "PF"), ] <- 0
weights_pp[!(model_dat$series == "PP"), ] <- 0
weights_rm[!(model_dat$series == "RM"), ] <- 0

data_all <- list(
  lag = lag,
  meantemp = model_dat$meantemp,
  meantemp_lag_1 = model_dat$meantemp_lag_1,
  mintemp = as.matrix(select(model_dat, mintemp_lag_0:mintemp_lag_5)),
  delta_mintemp = model_dat$delta_mintemp,
  mintemp_ma3 = model_dat$mintemp_ma3,
  maxtemp = as.matrix(select(model_dat, maxtemp_lag_1:maxtemp_lag_6)),
  maxtemp_ma3 = model_dat$maxtemp_ma3,
  warm_precip = model_dat$warm_precip,
  cool_precip = model_dat$cool_precip,
  ndvi_ma12 = model_dat$ndvi_ma12,
  weights_dm = weights_dm,
  weights_do = weights_do,
  weights_ds = weights_ds,
  weights_pb = weights_pb,
  weights_pe = weights_pe,
  weights_pf = weights_pf,
  weights_pp = weights_pp,
  weights_rm = weights_rm,
  y = model_dat$y,
  series = model_dat$series,
  time = model_dat$time,
  newmoonnumber = model_dat$newmoonnumber
)

saveRDS(data_all, file = "data_all.rds")

model_dat_pb_regime <- model_dat |>
  filter(!(series %in% c("PE","RM","DS")))
filter_indices <- as.numeric(setdiff(rownames(model_dat), rownames(model_dat_pb_regime)))

data_pb_regime <- list(
  lag = lag[-filter_indices, , drop = FALSE],
  meantemp = model_dat_pb_regime$meantemp,
  meantemp_lag_1 = model_dat_pb_regime$meantemp_lag_1,
  mintemp = as.matrix(select(model_dat, mintemp_lag_0:mintemp_lag_5)),
  delta_mintemp = model_dat$delta_mintemp,
  mintemp_ma3 = model_dat_pb_regime$mintemp_ma3,
  maxtemp = as.matrix(select(model_dat_pb_regime, maxtemp_lag_1:maxtemp_lag_6)),
  maxtemp_ma3 = model_dat_pb_regime$maxtemp_ma3,
  warm_precip = model_dat_pb_regime$warm_precip,
  cool_precip = model_dat_pb_regime$cool_precip,
  ndvi_ma12 = model_dat_pb_regime$ndvi_ma12,
  weights_dm = weights_dm[-filter_indices, , drop = FALSE],
  weights_do = weights_do[-filter_indices, , drop = FALSE],
  weights_pb = weights_pb[-filter_indices, , drop = FALSE],
  weights_pf = weights_pf[-filter_indices, , drop = FALSE],
  weights_pp = weights_pp[-filter_indices, , drop = FALSE],
  y = model_dat_pb_regime$y,
  series = droplevels(model_dat_pb_regime$series),
  time = model_dat_pb_regime$time,
  newmoonnumber = model_dat_pb_regime$newmoonnumber
)

saveRDS(data_pb_regime, file = "data_pb_regime.rds")

model_dat_heteromyids <- model_dat |>
  filter(!(series %in% c("PE","RM")))
filter_indices2 <- as.numeric(setdiff(rownames(model_dat), rownames(model_dat_heteromyids)))

data_heteromyid <- list(
  lag = lag[-filter_indices2, , drop = FALSE],
  meantemp = model_dat_heteromyids$meantemp,
  meantemp_lag_1 = model_dat_heteromyids$meantemp_lag_1,
  mintemp = as.matrix(select(model_dat, mintemp_lag_0:mintemp_lag_5)),
  delta_mintemp = model_dat$delta_mintemp,
  mintemp_ma3 = model_dat_heteromyids$mintemp_ma3,
  maxtemp = as.matrix(select(model_dat_heteromyids, maxtemp_lag_1:maxtemp_lag_6)),
  maxtemp_ma3 = model_dat_heteromyids$maxtemp_ma3,
  warm_precip = model_dat_heteromyids$warm_precip,
  cool_precip = model_dat_heteromyids$cool_precip,
  ndvi = model_dat_heteromyids$ndvi,
  ndvi_ma12 = model_dat_heteromyids$ndvi_ma12,
  winter_ndvi = model_dat_heteromyids$winter_ndvi,
  summer_ndvi = model_dat_heteromyids$summer_ndvi,
  weights_dm = weights_dm[-filter_indices2, , drop = FALSE],
  weights_do = weights_do[-filter_indices2, , drop = FALSE],
  weights_ds = weights_ds[-filter_indices2, , drop = FALSE],
  weights_pb = weights_pb[-filter_indices2, , drop = FALSE],
  weights_pf = weights_pf[-filter_indices2, , drop = FALSE],
  weights_pp = weights_pp[-filter_indices2, , drop = FALSE],
  y = model_dat_heteromyids$y,
  series = droplevels(model_dat_heteromyids$series),
  time = model_dat_heteromyids$time,
  newmoonnumber = model_dat_heteromyids$newmoonnumber
)
saveRDS(data_heteromyid, file = "data_heteromyid.rds")
