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
fit_mod <- function(train_data, test_data) {
  mvgam::mvgam(
    abundance ~ 1,  # model abundance with intercept only
    trend_formula = ~ -1,  # remove default smooth trend
    trend_model = "AR1",  # AR1 temporal correlation
    family = poisson(link = "log"),
    data = train_data,
    newdata = test_data,
    priors = prior(normal(0, 2), class = "Intercept"),
    chains = 4,
    samples = 200
  )
}

# Define parameters

sampling_mos=12
step_size=5 #make datasets increase by 5 years at a time
training_increment <- sampling_mos * step_size          # training increment in time steps (e.g., months)
forecast_horizon <- 12     # number of steps ahead to forecast (3 for short, 12 for long)

# Prepare result list

mod1_results <- list()

# Loop through incremental training windows
split_id <- 1
for (train_end in seq(from = training_increment,
                      to = nrow(dmdat) - forecast_horizon,
                      by = training_increment)) {

  # Define training and test datasets
  train_data <- dmdat[1:train_end, ]
  test_data <- dmdat[(train_end + 1):(train_end + forecast_horizon), ]


  # Ensure test data exists and is complete
  if (nrow(test_data) == forecast_horizon) {

    # Fit the model
    model <- fit_mod(train_data, test_data)

    # forecasts
    preds<- forecast(model, newdata = test_data, type = 'response')

    # scores
    forecast_scores <- mvgam::score(preds, score = 'crps')

    # Store results
    mod1_results[[split_id]] <- list(
      split_id = split_id,
      train_end = train_end,
      test_start = train_end + 1,
      test_end = train_end + forecast_horizon,
      model = model,
      preds = preds,
      scores = forecast_scores
    )

    split_id <- split_id + 1
  }
}

# Compile CRPS scores into a dataframe
crps_df <- map_dfr(mod1_results, ~{
  tibble(
    split_id = .x$split_id,
    time = .x$test_start:.x$test_end,
    date = dmdat$censusdate[.x$test_start:.x$test_end],
    horizon = 1:12,  # step-ahead horizon from 1 to 12
    crps = .x$scores$DM$score
  )
})

ggplot(crps_df, aes(x = as.factor(horizon), y = crps, col=split_id)) +
  geom_violin(color = "black", alpha = 1) +
  geom_jitter(width = 0.2, alpha = 1, size = 2, aes(color = split_id)) +
 # scale_color_viridis_d() +
  labs(title = "CRPS Distribution by Forecast Horizon",
       x = "Forecast Horizon (months)",
       y = "CRPS") +
  theme_classic()

ggplot(crps_df, aes(x = date, y = crps)) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue") +
  labs(title = "CRPS Over Time",
       x = "Time (newmoonnumber)",
       y = "CRPS") +
  theme_classic()
