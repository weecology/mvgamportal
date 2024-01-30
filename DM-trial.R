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


#generate data subsets####

rodent_data=summarize_rodent_data(
  path = get_default_data_path(),
  clean = TRUE,
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
  filter(treatment=="control")%>%
  select(censusdate,newmoonnumber,DM)

covars=weather(level="newmoon", fill=TRUE, horizon=365, path=get_default_data_path())%>%
  select(newmoonnumber,meantemp, mintemp, maxtemp, precipitation, warm_precip, cool_precip)

#up to Dec 2019 only
dmcont_covs=right_join(covars,dmcont_dat, by="newmoonnumber")

dmdat=dmcont_covs%>%
  rename("abundance"="DM")%>%filter(!newmoonnumber>526)%>%
  mutate(time=newmoonnumber - min(newmoonnumber) + 1)%>%
  mutate(series=as.factor('DM'))%>%
  select(time, series, abundance, meantemp,warm_precip, cool_precip)


#shorter subset
dmdat_sample=dmdat%>%filter(!time<1, !time>46)

dmdat1=rolling_origin(
  data       = dmdat_sample, #all DM control data
  initial    = 12, #samples used for modelling (training)
  assess     = 24, # number of samples used for each assessment resample (horizon)
  cumulative = TRUE #length of analysis set is growing; 
)

#function for fitting AR1 model, getting forecasts, and scoring them
#output is a list
#for testing purposes, set a short burn-in period

fit_cast_score=function(split) {
  
  data_train= analysis(split) #training data
  
  data_test= assessment(split) # test data
  
  model= mvgam(abundance~1, 
               trend_model="AR1",
               family= poisson(link = "log"), 
               data=data_train,
               newdat= data_test,
               chains = 4,
               burnin = 100)
  
  preds= as.vector(forecast(model, data_test))
  
  get_score= score(preds)
  
  return(list(model, preds, get_score))
  
}

#fit, predict, score

dmdat1$output=map(dmdat1$splits, fit_cast_score)

