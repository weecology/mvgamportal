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

#up to Dec 2019 only
dmcont_covs=right_join(covars,dmcont_dat, by="newmoonnumber")

dmdat=dmcont_covs%>%
  rename("abundance"="DM")%>%filter(!newmoonnumber>526)%>%
  mutate(time=newmoonnumber - min(newmoonnumber) + 1)%>%
  mutate(series=as.factor('DM'))%>%
  select(time, censusdate, newmoonnumber, series, abundance, meantemp,warm_precip, cool_precip)%>%
  arrange(time)

#apply sliding-index to create subsets of training data at different windows

dmdat20=sliding_index(
  data= dmdat, #all DM control data
  newmoonnumber,
  lookback=240,
  assess_stop=12,
  complete = TRUE
)

#function for fitting AR1 model, getting forecasts, and scoring them
#output is a list
#for testing purposes, set fewer iters
#trend formula as suggested by Nick

fitmod1_cast_score=function(split) {

  data_train= analysis(split) #training data

  data_test= assessment(split) # test data

  model= mvgam(abundance~1,
               trend_formula = ~ -1,
               trend_model="AR1",
               family= poisson(link = "log"),
               data=data_train,
               newdata= data_test,
               priors = prior(normal(0, 2), class = Intercept),
               chains = 4,
               samples = 200)

  preds= as.vector(forecast(model, data_test))

  get_score= score(preds)

  return(list(model, preds, get_score))

}

#shorter subset moving window for different lengths of training data
dmdat20v1=dmdat20[1:5,]

#fit, predict, score
#dmdat20v1$output=map(dmdat20v1$splits, fitmod1_cast_score)

saveRDS(dmdat20v1, "mvgam-sample-output.RDS")
