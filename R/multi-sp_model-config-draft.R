#preliminary analyses
#model configuration for timeseries length project

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

multi_sp_cont_dat=rodent_data%>%
  filter(treatment%in%c("control", NA))%>%
  select(censusdate,newmoonnumber,DM, PP, OT, PE, PB, RM, DO)

covars=weather(level="newmoon", fill=TRUE, horizon=365, path=get_default_data_path())%>%
  select(newmoonnumber,meantemp, mintemp, maxtemp, precipitation, warm_precip, cool_precip)%>%
  mutate(meantemp_lag1=lag(meantemp,order_by=newmoonnumber))

#select PP period:Sept 2010- Dec 2019
#zeros in cool_precip causing break
multi_sp_cont_covs=right_join(covars,multi_sp_cont_dat, by="newmoonnumber")%>%
  pivot_longer(cols=10:16, names_to="species")

multi_sp_dat=multi_sp_cont_covs%>%
  rename("abundance"="value")%>%filter(!newmoonnumber<411, !newmoonnumber>526)%>%
  mutate(month = lubridate::month(censusdate),
         year = lubridate::year(censusdate)) %>%
  mutate(time=newmoonnumber - min(newmoonnumber) + 1, series=as.factor(species))%>%
  select(time, censusdate,month, year, newmoonnumber, series, abundance, meantemp_lag1,warm_precip, cool_precip)%>%
  arrange(time)

n_moons_yr=12
n_yrs=5
nsp=length(unique(multi_sp_dat$series))
n_moons_train=n_moons_yr*n_yrs*nsp
n_moons_test=n_moons_yr*1*nsp

data_train=multi_sp_dat[1:n_moons_train,]
data_test=multi_sp_dat[n_moons_train+1:n_moons_test,]

#not sure why but it is important to treat species as factors

multisp_modvar=mvgam(formula= abundance~1,
               trend_model='VAR1cor',
               family= poisson(link = "log"),
               data=data_train,
               newdata= data_test,
               trend_formula=  ~ s(warm_precip, trend, bs="re")+
                 s(cool_precip, trend, bs="re")+
                 s(meantemp_lag1, trend, bs="re"), # process model formula, which includes the smooth functions,
               #trend is a latent process
               chains = 4,
               samples = 2000)

#inspect results

summary(multisp_modvar)
plot_mvgam_fc(multisp_modvar, series=1)
