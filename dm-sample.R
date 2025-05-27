#preliminary analyses for timeseries length and forecasting project

#load packages####

#for mvgam installation refer to: https://course.naturecast.org/lessons/r-complex-time-series-models-1/material/

require(portalr)
require(mvgam)
require(tidyr)
require(ggpubr)
require(tidyverse)
require(dplyr)
require(vctrs)
require(lubridate)
require(rsample)

#options(mc.cores = parallel::detectCores())

#generate data subsets####

rodent_data=summarize_rodent_data(
  path = get_default_data_path(),
  clean = FALSE,
  level="Treatment",
  type = "Rodents",
  plots = "Longterm",
  unknowns = FALSE,
  shape = "crosstab",
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

dmcont_dat=rodent_data%>%
  filter(treatment%in%c("control", NA))%>%
  select(censusdate,newmoonnumber,DM)

covars=weather(level="newmoon", fill=TRUE, horizon=365, path=get_default_data_path())%>%
  select(newmoonnumber,meantemp, mintemp, maxtemp, precipitation, warm_precip, cool_precip)

dmcont_covs=right_join(covars,dmcont_dat, by="newmoonnumber")

dmdat=dmcont_covs%>%
  rename("abundance"="DM")%>%filter(!newmoonnumber>526)%>%
  mutate(time=newmoonnumber - min(newmoonnumber) + 1)%>%
  mutate(series=as.factor('DM'))%>%
  select(time, censusdate, newmoonnumber, series, abundance, meantemp,warm_precip, cool_precip)%>%
  arrange(time)%>%
  filter(!time<32)

# Define the model-fitting function
fit_mod_AR <- function(train_data, test_data) {
  mvgam::mvgam(
    abundance ~ 1,  # model abundance with intercept only
    trend_formula = ~ -1,  # remove default smooth trend
    trend_model = "AR1",  # AR1 temporal correlation
    family = poisson(link = "log"),
    data = train_data,
    newdata = test_data,
    priors = prior(normal(0, 2), class = "Intercept"),
    chains = 4,
    samples = 2000
  )
}

fit_mod_ARlin=function(train_data, test_data) {
  mvgam::mvgam(formula= abundance~1,
               trend_model='VAR1cor',
               family= poisson(link = "log"),
               data=train_data,
               newdata= test_data,
               trend_formula=  ~ s(warm_precip, trend, bs="re")+
                 s(cool_precip, trend, bs="re")+
                 s(meantemp_lag1, trend, bs="re"), # process model formula, which includes the smooth functions,
               #trend is a latent process
               chains = 4,
               samples = 2000)
}

step_size <- 5                       # How much to grow training set each step
sampling_mos <- 12                  # Test set size in months
training_increment <- sampling_mos * step_size
forecast_horizon <- 12

max_train_length <- 480             # Maximum training data length

max_index <- nrow(dmdat)

# Prepare result list
mod1_results <- list()
split_id <- 1

# Start test set 12 months before end, and move backward
for (test_end_index in seq(max_index, sampling_mos + training_increment, by = -1)) {

  test_start_index <- test_end_index - sampling_mos + 1

  # Only proceed if there's room for a training set before test set
  if ((test_start_index - 1) >= training_increment) {

    # Loop through training sizes, ending at test_start - 1
    for (train_length in seq(from = training_increment,
                             to = min(test_start_index - 1, max_train_length),
                             by = training_increment)) {

      train_start_index <- test_start_index - train_length
      train_end_index <- test_start_index - 1

      # Check bounds
      if (train_start_index >= 1) {

        train_data <- dmdat[train_start_index:train_end_index, ]
        test_data <- dmdat[test_start_index:test_end_index, ]

        # model, prediction, scoring
        model <- fit_mod_AR(train_data, test_data)
        preds <- forecast(model, newdata = test_data, type = 'response')
        forecast_scores <- mvgam::score(preds, score = 'crps')

        # Store results
        mod1_results[[split_id]] <- list(
          split_id = split_id,
          train_start = train_start_index,
          train_end = train_end_index,
          test_start = test_start_index,
          test_end = test_end_index,
          model = model,
          preds = preds,
          scores = forecast_scores
        )

        split_id <- split_id + 1
      }
    }
  }
}

saveRDS(mod1_results, "mod_AR_results.RDS")

#check list if it looks correct
results_df <- do.call(rbind, lapply(mod1_results, function(x) {
  data.frame(
    split_id = x$split_id,
    train_start = x$train_start,
    train_end = x$train_end,
    test_start = x$test_start,
    test_end = x$test_end,
    training_length = x$train_end - x$train_start + 1
  )
}))

# plot
par(mfrow = c(1, 3))

# Histogram of train_end indices
hist(results_df$train_end,
     breaks = 30,
     main = "Distribution of Train End Indices",
     xlab = "Train End Index",
     col = "skyblue", border = "white")
abline(v = 60, col = "red", lty = 2, lwd = 2)
text(x = 60, y = par("usr")[4] * 0.9, labels = "60", col = "red", pos = 4)
abline(v = 483, col = "red", lty = 2, lwd = 2)
text(x = 483, y = par("usr")[4] * 0.9, labels = "483", col = "red", pos = 4)

# Histogram of test start indices
hist(results_df$test_start,
     breaks = 30,
     main = "Distribution of Test Start Indices",
     xlab = "Test Start Index",
     col = "skyblue", border = "white")
abline(v = 61, col = "red", lty = 2, lwd = 2)
text(x = 61, y = par("usr")[4] * 0.9, labels = "61", col = "red", pos = 4)
abline(v = 484, col = "red", lty = 2, lwd = 2)
text(x = 484, y = par("usr")[4] * 0.9, labels = "484", col = "red", pos = 4)


# Histogram of training data lengths
hist(results_df$training_length,
     breaks = 20,
     main = "Distribution of Training Set Lengths",
     xlab = "Training Length",
     col = "lightgreen", border = "white")

# Extract training lengths
training_lengths <- sapply(mod1_results, function(x) {
  x$train_end - x$train_start + 1
})

# Get shortest and longest lengths
min(training_lengths)
max(training_lengths)
# Compile CRPS scores into a dataframe

crps_df <- map_dfr(mod1_results, ~{
  n_scores <- length(.x$scores$DM$score)  # number of forecast horizons scored
  train_len <- .x$train_end - .x$train_start + 1
  test_range <- .x$test_start:.x$test_end

  tibble(
    split_id = .x$split_id,
    train_len = train_len,
    time = test_range[1:n_scores],
    date = dmdat$censusdate[test_range][1:n_scores],
    horizon = seq_len(n_scores),
    crps = .x$scores$DM$score
  )
})

#plot
p1=ggplot(crps_df, aes(x = as.factor(horizon), y = crps, color =train_len)) +
  geom_violin(alpha = 1, color="black") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1, aes(color = train_len)) +
  scale_colour_viridis_c()+
  labs(
    title = "CRPS Distribution by Forecast Horizon",
    x = "Forecast Horizon",
    y = "CRPS",
    color = "TS length"
  ) +
  theme_classic()

p2=ggplot(crps_df, aes(x=train_len, y=crps, col=horizon))+
  #geom_line(aes(x=train_len, y=crps, color = as.factor(horizon)))+
   geom_point(alpha = 1, size = 2, aes(x=train_len, y=crps, color = as.factor(horizon)))+
  theme_classic()+scale_color_viridis_d()+
  xlab("training data length")+ggtitle("DM (AR1 model)")+
  facet_wrap(~horizon)

ggarrange(p1,p2)
