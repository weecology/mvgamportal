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

options(mc.cores = parallel::detectCores())

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

dmdat2=sliding_index(
  data= dmdat, #all DM control data
  newmoonnumber,
  lookback=24,
  assess_stop=12,
  complete = TRUE
)

dmdat5=sliding_index(
  data= dmdat, #all DM control data
  newmoonnumber,
  lookback=60,
  assess_stop=12,
  complete = TRUE
)

dmdat10=sliding_index(
  data= dmdat, #all DM control data
  newmoonnumber,
  lookback=120,
  assess_stop=12,
  complete = TRUE
)

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

fit_cast_score=function(split) {

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
dmdat20v1$output=map(dmdat20v1$splits, fitmod1_cast_score)
#dm1=lapply(dmdat20v1$splits, fitmod1_cast_score)

#access results

#predictions
N=length(dmdat20v1$id)

preds= c()

for (i in 1:N) {

  item = paste("forecast_", i)
  preds[[item]]= dmdat20v1$output[[i]][[2]]$forecasts$DM
}

#get scores
score= c()

for (i in 1:N) {

  item = paste("score_", i)
  score[[item]]= dmdat20v1$output[[i]][[3]]$DM
}

###################
#this part not needed; just for manual checking and stuff
#predictions:

d1=as.data.frame(dmdat20v1[[3]][[1]][[2]]$forecasts$DM)
d2=as.data.frame(dmdat20v1[[3]][[2]][[2]]$forecasts$DM)
d3=as.data.frame(dmdat20v1[[3]][[3]][[2]]$forecasts$DM)

#scores:

s1=as.data.frame(dmdat20v1[[3]][[1]][[3]]$DM)
s2=as.data.frame(dmdat20v1[[3]][[2]][[3]]$DM)
s3=as.data.frame(dmdat20v1[[3]][[3]][[3]]$DM)


#assess model performance/convergence

m1=dmdat20v1[[3]][[1]][[1]]$model_output
m2=dmdat20v1[[3]][[2]][[1]]$model_output
m3=dmdat20v1[[3]][[3]][[1]]$model_output
