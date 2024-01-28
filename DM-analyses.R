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
  rename("abundance"="DM")%>%
  mutate(time= seq(1:length(dmcont_covs$censusdate)), series=as.factor('DM'))%>%
  select(time, series, abundance, meantemp,warm_precip, cool_precip)

#create rolling origin object for analysis

#generate short-term forecasts for now (1-3 time steps ahead)

n_moons_yr=12
n_yrs=2 #min. TS length
n_moons_train=n_moons_yr*n_yrs
n_moons_test=3

DMcontrol_dat <- 
  rolling_origin(
    data       = dmdat, #all DM control data
    initial    = n_moons_train, #samples used for modelling (training)
    assess     = n_moons_test, # number of samples used for each assessment resample (horizon)
    cumulative = TRUE #length of analysis set is growing; 
  )

plot_mvgam_series(data = dmdat, y = 'abundance')

#fit AR1 model####

fit_ar1mod=function(split) {
  
  data_train= analysis(split) #get dataframe
  
  fit_model= mvgam(abundance~1, 
                   trend_model="AR1",
                   family= poisson(link = "log"), 
                   data=data_train)
}

DMcontrol_dat$ar1mod=map(DMcontrol_dat$splits, fit_ar1mod)

#save(DMcontrol_dat$ar1mod, file="DMcontrol_dat_ar1mod.RData") #save output

#check model convergence####
mcmc_plot(DMcontrol_dat$ar1mod[[1]], type = "trace", variable =  "ar1[1]")

# generate forecasts####

get_preds= function(split, model) {
  
  data_train=analysis(split)
  data_test=assessment(split)
  
  preds=forecast(model, newdata=data_test)
}

DMcontrol_dat$ar1preds=pmap(list(DMcontrol_dat$splits, DMcontrol_dat$ar1mod), get_preds)

#saveRDS(DMcontrol_dat$ar1pred, file="DMcontrol_dat_ar1preds.RDS")
