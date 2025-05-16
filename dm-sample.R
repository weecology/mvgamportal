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
               samples = 200)
}
step_size <- 5
sampling_mos <- 12  # 1 year
training_increment <- sampling_mos * step_size  # How much to grow training set each time

forecast_horizon <- 12
test_start_index <- 484 #will need to make this dynamic if we were to slide
test_end_index <- 495

# Prepare result list
mod1_results <- list()

# Loop through different training window sizes, ending at the point before test
split_id <- 1
for (train_length in seq(from = training_increment,
                         to = test_start_index - 1,
                         by = training_increment)) {

  # Define start and end of training data to make sure it immediately precedes test data
  train_start_index <- test_start_index - train_length
  train_end_index <- test_start_index - 1

  # Ensure valid range
  if (train_start_index >= 1) {

    # Subset training and test data
    train_data <- dmdat[train_start_index:train_end_index, ]
    test_data <- dmdat[test_start_index:test_end_index, ]

    # Fit model
    model <- fit_mod_AR(train_data, test_data)

    # Forecast
    preds <- forecast(model, newdata = test_data, type = 'response')

    # Score
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
  geom_jitter(width = 0.2, alpha = 1, size = 2, aes(color = train_len)) +
  scale_colour_viridis_c()+
  labs(
    title = "CRPS Distribution by Forecast Horizon",
    x = "Forecast Horizon",
    y = "CRPS",
    color = "TS length"
  ) +
  theme_classic()

p2=ggplot(crps_df, aes(x=train_len, y=crps, col=horizon))+
  geom_point(width = 0.2, alpha = 1, size = 2, aes(color = as.factor(horizon)))+
  geom_line(aes(color = as.factor(horizon)))+
  theme_classic()+scale_color_viridis_d()+
  xlab("training data length")+ggtitle("DM (AR1 model)")

ggarrange(p1,p2)
