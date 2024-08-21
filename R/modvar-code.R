#### 1. Prep data for modelling ####
library(portalcasting)
library(mgcv)
library(dplyr)
library(cmdstanr)
library(forecast)
library(mvgam)
library(portalr)

# Load the most recent rodents survey table for control plots
moon_dates=read.csv("https://raw.githubusercontent.com/weecology/PortalData/main/Rodents/moon_dates.csv")

covars=weather(level="newmoon", fill=TRUE, horizon=365, path=get_default_data_path())%>%
  select(newmoonnumber, meantemp, mintemp, maxtemp, precipitation, warm_precip, cool_precip)%>%
  mutate(meantemp_lag1=lag(meantemp,order_by=newmoonnumber))

ndvi_dat=ndvi(level="newmoon")

rodent_data=summarize_rodent_data(
  path = get_default_data_path(),
  clean = TRUE,
  level="Treatment",
  type = "Rodents",
  plots = "Longterm",
  unknowns = FALSE,
  shape = "crosstab",
  time = "newmoon",
  output = "abundance",
  na_drop = FALSE,
  zero_drop = FALSE,
  min_traps = 1,
  min_plots = 1,
  effort = TRUE,
  download_if_missing = TRUE,
  quiet = FALSE
)

rodents_table <- left_join(rodent_data,moon_dates)%>%
  left_join(covars)%>%
  left_join(ndvi_dat)%>%
  filter(!newmoonnumber<85, treatment=="control")

# Calculate means and sds of covariates for later unscaled plotting
rodents_table %>%
  dplyr::mutate(month = lubridate::month(newmoondate),
                year = lubridate::year(newmoondate)) %>%
  dplyr::filter(year > 1983) %>%
  dplyr::select(mintemp) %>%
  dplyr::summarise(mintemp_mean = mean(mintemp, na.rm = TRUE),
                   mintemp_sd = sd(mintemp, na.rm = TRUE)) -> mintemp_stats

rodents_table %>%
  dplyr::mutate(month = lubridate::month(newmoondate),
                year = lubridate::year(newmoondate)) %>%
  dplyr::filter(year > 1983) %>%
  dplyr::select(ndvi) %>%
  dplyr::summarise(ndvi_mean = mean(ndvi, na.rm = TRUE),
                   ndvi_sd = sd(ndvi, na.rm = TRUE)) -> ndvi_stats

# Prep for modelling
rodents_table %>%
  dplyr::mutate(month = lubridate::month(newmoondate),
                year = lubridate::year(newmoondate)) %>%
  dplyr::filter(year > 1983) %>%
  # Scale continuous variables for massively improved efficiency of
  # Stan sampling
  dplyr::mutate(ndvi = as.vector(scale(ndvi)),
                mintemp = as.vector(scale(mintemp)),
                maxtemp = as.vector(scale(maxtemp))) %>%
  dplyr::mutate(ndvi_ma12 = zoo::rollmean(ndvi, k = 12, align = 'right',
                                          na.pad = TRUE)) %>%
  # # Keep the first observation if multiple taken in the same month
  # dplyr::arrange(year, month) %>%
  # dplyr::group_by(month, year) %>%
  # dplyr::slice_head(n = 1) %>%
  tidyr::pivot_longer(cols = colnames(rodents_table)[5:25],
                      names_to = 'series', values_to = 'y') %>%
  dplyr::select(y, series, month, year,
                newmoonnumber, mintemp:ndvi_ma12) %>%
  dplyr::mutate(time = newmoonnumber - (min(newmoonnumber) - 1))-> model_dat

# Many models will fail if the series of observations is nearly all zeroes.
# Remove species with < 33 (out of 330) total unique observations (i.e. captures in at least
# 10% of unique trapping sessions)
# as forecasting these is not really useful anyway
model_dat %>%
  dplyr::group_by(series) %>%
  dplyr::summarise(total_obs = length(which(y >= 1))) %>%
  dplyr::filter(total_obs >= 33) %>%
  dplyr::pull(series) -> series_keep

model_dat %>%
  dplyr::filter(series %in% series_keep) %>%
  dplyr::filter(series != 'total') %>%
  dplyr::mutate(series = as.factor(series)) %>%
  dplyr::arrange(time, series) -> model_dat

#UNSURE WHY THIS IS FALSE FOR NOW
(max(model_dat$time) * length(unique(model_dat$series))) == NROW(model_dat)

# Feature engineering
#1. Distributed lag matrices for environmental covariates
# Function to set up a lag matrix for distributed lag nonlinear models
lagard <- function(x, n_lag = 6){
  n <- length(x)
  X <- matrix(NA, n, n_lag)
  for (i in 1:n_lag) X[i:n, i] <- x[i:n - i + 1]
  X
}

# Function to generate predictions for missing real-valued environmental variables
# using a GAM with seasonality and yearly components
approx_gam = function(df, family = gaussian()){
  require(mgcv)
  mod <- bam(y ~
               s(month, bs = 'cc', k = 10) +
               s(year, bs = 'bs', m = c(2,1,0),
                 k = 12), data = df, discrete = TRUE)
  preds <- predict(mod, newdata = df, type = 'response')

  # Replace any missing values with model-based predictions
  truth <- df$y
  truth[is.na(truth)] <- preds[is.na(truth)]
  truth
}

# Mintemp 6-month lag matrix
unique_times <- sort(unique(model_dat$time))
mintemp <- lagard(approx_gam(model_dat %>%
                               dplyr::select(mintemp, month, year, time) %>%
                               dplyr::arrange(time) %>%
                               dplyr::distinct() %>%
                               dplyr::mutate(y = mintemp)), 6)
mintemp_df <- data.frame(mintemp)
mintemp_df$time <- unique_times
model_dat %>%
  dplyr::select(time, year, month) %>%
  dplyr::left_join(mintemp_df) -> mintemp_df
dim(mintemp_df)[1] == NROW(model_dat)

# Maxtemp 6-month lag matrix
maxtemp <- lagard(approx_gam(model_dat %>%
                               dplyr::select(maxtemp, month, year, time) %>%
                               dplyr::arrange(time) %>%
                               dplyr::distinct() %>%
                               dplyr::mutate(y = maxtemp)), 6)
maxtemp_df <- data.frame(maxtemp)
maxtemp_df$time <- unique_times
model_dat %>%
  dplyr::select(time, year, month) %>%
  dplyr::left_join(maxtemp_df) -> maxtemp_df
dim(maxtemp_df)[1] == NROW(model_dat)

# The lag matrix
lag <- matrix(0:5, nrow(model_dat),
              6, byrow = TRUE)
dim(lag)[1] == NROW(model_dat)

#2. Create remaining moving average / anomaly versions of environmental covariates
model_dat %>%
  dplyr::left_join(model_dat %>%
                     dplyr::select(time, mintemp, maxtemp, ndvi) %>%
                     dplyr::distinct() %>%
                     dplyr::mutate(mintemp_ma3 = zoo::rollmean(mintemp, k = 3, align = 'right',
                                                               na.pad = TRUE),
                                   maxtemp_ma3 = zoo::rollmean(maxtemp, k = 3, align = 'right',
                                                               na.pad = TRUE))) -> model_dat

# As we now have NAs for the first 11 rows of observations for each lag matrix,
# as well as NAs for some rows of the moving average covariates,
# filter the data so that no NAs remain for covariates
model_dat %>%
  dplyr::filter(time > 11) %>%
  dplyr::mutate(time = time - 11) -> model_dat

# Impute ndvi_ma12
model_dat %>%
  dplyr::select(-ndvi_ma12) %>%
  dplyr::left_join(model_dat %>%
                     dplyr::select(time) %>%
                     dplyr::distinct() %>%
                     dplyr::arrange(time) %>%
                     dplyr::bind_cols(data.frame(ndvi_ma12 = approx_gam(model_dat %>%
                                                                          dplyr::select(ndvi_ma12, month, year, time) %>%
                                                                          dplyr::arrange(time) %>%
                                                                          dplyr::distinct() %>%
                                                                          dplyr::mutate(y = ndvi_ma12))))) -> model_dat


mintemp_df %>%
  dplyr::ungroup() %>%
  dplyr::filter(time > 11) %>%
  dplyr::select(-time, -year, -month) %>%
  as.matrix() -> mintemp
dim(mintemp)[1] == NROW(model_dat)

maxtemp_df %>%
  dplyr::ungroup() %>%
  dplyr::filter(time > 11) %>%
  dplyr::select(-time, -year, -month) %>%
  as.matrix() -> maxtemp
dim(maxtemp)[1] == NROW(model_dat)

lag <- tail(lag, NROW(model_dat))
dim(lag)[1] == NROW(model_dat)

# Now create weight matrices that can be used for setting up hierarchical
# distributed lag terms
weights_dm <- weights_do <-
  weights_pp <- weights_ol <-
  weights_ot <- weights_pf <-
  weights_pb <- weights_pe <- weights_rm <-
  matrix(1, ncol = ncol(lag), nrow = nrow(lag))

weights_dm[!(model_dat$series == 'DM'), ] <- 0
weights_do[!(model_dat$series == 'DO'), ] <- 0
weights_ol[!(model_dat$series == 'OL'), ] <- 0
weights_ot[!(model_dat$series == 'OT'), ] <- 0
weights_pb[!(model_dat$series == 'PB'), ] <- 0
weights_pe[!(model_dat$series == 'PE'), ] <- 0
weights_pf[!(model_dat$series == 'PF'), ] <- 0
weights_pp[!(model_dat$series == 'PP'), ] <- 0
weights_rm[!(model_dat$series == 'RM'), ] <- 0


# Create a list to store the full dataset, including lag matrices and
# moving averages for the environmental covariates
data_all <- list(lag = lag,
                 mintemp = mintemp,
                 mintemp_ma3 = model_dat$mintemp_ma3,
                 maxtemp = maxtemp,
                 maxtemp_ma3 = model_dat$maxtemp_ma3,
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
                 month = model_dat$month,
                 year = model_dat$year,
                 series = model_dat$series,
                 time = model_dat$time)

# Split data into training and testing; stop training
# at the end of 2018 so that 2019 can be evaluated. Conditions were
# challenging in COVID and post-COVID, so evaluation of models may not be
# as 'fair'
train_inds <- which(model_dat$year < 2019)

data_train <- lapply(seq_along(data_all), function(x){
  if(is.matrix(data_all[[x]])){
    data_all[[x]][train_inds,]
  } else {
    data_all[[x]][train_inds]
  }
})

test_inds <- which(model_dat$year %in% c(2019))

data_test <- lapply(seq_along(data_all), function(x){
  if(is.matrix(data_all[[x]])){
    data_all[[x]][test_inds,]
  } else {
    data_all[[x]][test_inds]
  }
})
names(data_train) <- names(data_test) <- names(data_all)

# Save all objects for forecast modelling
save(model_dat,
     data_all,
     data_train,
     data_test,
     file = 'rodents_data_tsobjects.rda')


# Load the pre-prepared modelling data
load('rodents_data_tsobjects.rda')

# Source the bespoke checking / graphical functions
source('D:/Dropbox (UFL)/PhD-stuff/mvgamportal/R/checking_functions.R')

# View some of the raw time series
plot_mvgam_series(data = data_train, series = 'all')
plot_mvgam_series(data = data_train, newdata = data_test, series = 8)

# With prior distributions derived, we can construct the mvgam model with a latent
# VAR1 temporal process
priors <- get_mvgam_priors(formula = y ~ 1,
                           trend_formula = ~ s(ndvi_ma12, trend, bs = 're') +
                             te(mintemp, lag, k = c(3, 4), bs = c('tp', 'cr')) +
                             te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
                             te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
                             te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
                             te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')),
                           data = data_train,
                           family = poisson(),
                           trend_model = 'VAR1cor')

# Update the prior for the species-level process model variances, with appropriate upper and
# lower bounds. Here we use brms functionality to easily set priors
priors <- prior(beta(10, 10), class = sigma, lb = 0.2, ub = 1)

# Update the prior for the NDVI random slopes
priors <- c(priors,
            prior(inv_gamma(2.3693353, 0.7311319), class = sigma_raw_trend))

# The observation-level intercept is a nuisance parameter in this model that we
# don't really need (and cannot adequately estimate anyway).
# At present there is no way to drop this term using mvgam, but we can
# heavily regularize it to zero with a strong prior
priors <- c(priors,
            prior(normal(0, 0.001), class = Intercept))

priors

# Ensure the priors are properly incorporated into the model's Stan code
code(mvgam(formula = y ~ 1,
           trend_formula = ~ s(ndvi_ma12, trend, bs = 're') +
             te(mintemp, lag, k = c(3, 4), bs = c('tp', 'cr')) +
             te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
             te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
             te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
             te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')),
           data = data_train,
           newdata = data_test,
           family = poisson(),
           trend_model = 'VAR1cor',
           priors = priors,
           run_model = FALSE))

# Fit the model
modvar <- mvgam(formula = y ~ 1,
                trend_formula = ~ s(ndvi_ma12, trend, bs = 're') +
                  te(mintemp, lag, k = c(3, 4), bs = c('tp', 'cr')) +
                  te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c('tp', 'cr')) +
                  te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c('tp', 'cr')) +
                  te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c('tp', 'cr')) +
                  te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c('tp', 'cr')),
                data = data_train,
                newdata = data_test,
                family = poisson(),
                trend_model = 'VAR1cor',
                priors = priors,
                samples = 1600)

par(mfrow=c(3,1))

plot(modvar, forecast, series=1)
plot(modvar, forecast, series=2)
plot(modvar, forecast, series=3)
plot(modvar, forecast, series=4)
plot(modvar, forecast, series=5)
plot(modvar, forecast, series=6)
plot(modvar, forecast, series=7)
plot(modvar, forecast, series=8)
plot(modvar, forecast, series=9)

#### Fit full models to the entire series of data ####
modvar_all <- update(modvar, data = data_all, samples = 1600)
save(modvar_all, file = 'modvar_all.rda')

modvar_scores <- exp(log(score(forecast(modvar), score = 'variogram',
                               log = TRUE)$all_series$score[1:12] *
                           score(forecast(modvar), score = 'energy',
                                 log = TRUE)$all_series$score[1:12]))
