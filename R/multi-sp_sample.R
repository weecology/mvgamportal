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

multi_sp_cont_dat=rodent_data%>%
  filter(treatment%in%c("control", NA))%>%
  select(censusdate,newmoonnumber,DM, PP, OT, PE, PB, RM, DO)

covars=weather(level="newmoon", fill=TRUE, horizon=365, path=get_default_data_path())%>%
  select(newmoonnumber,meantemp, mintemp, maxtemp, precipitation, warm_precip, cool_precip)

#Jan 1980-Dec 2019 only
multi_sp_cont_covs=right_join(covars,multi_sp_cont_dat, by="newmoonnumber")%>%
  pivot_longer(cols=9:15, names_to="species")

multi_sp_cont_covs$series=as.factor(multi_sp_cont_covs$species)

multi_sp_dat=multi_sp_cont_covs%>%
  rename("abundance"="value")%>%filter(!newmoonnumber>526)%>%
  mutate(time=newmoonnumber - min(newmoonnumber) + 1)%>%
  select(time, censusdate, newmoonnumber, series, abundance, meantemp,warm_precip, cool_precip)%>%
  arrange(time)%>%filter(!newmoonnumber<33)

#apply sliding-window to create subsets of training data at different windows
#set 5 years of training data, rolling over every year (every 12 months)
# set 1 year forward of test data

multi_sp_dat5=sliding_index(
  data= multi_sp_dat, #all DM, DO, PP, PB, OT, PE, RM control data
  index=newmoonnumber,
  lookback=60,
  assess_stop=12,
  complete = TRUE,
  step=12
)

#function for fitting VAR model, getting forecasts, and scoring them
#output is a list
#for testing purposes, set fewer iters

fitmod2_cast_score=function(split) {

  data_train= analysis(split) #training data

  data_test= assessment(split) # test data

  model= mvgam(formula= abundance~1,
               trend_model='VAR1cor',
               family= poisson(link = "log"),
               data=data_train,
               newdata= data_test,
               trend_formula=  ~ s(warm_precip, trend, bs="re"), #what is trend doing? where is it from?
              # priors = prior(normal(0, 2), class = Intercept),
               chains = 3,
               samples = 100)

  preds= as.vector(forecast(model, data_test))

  get_score= score(preds)

  return(list(model, preds, get_score))
}

#shorter subset moving window for different lengths of training data
multi_sp_dat5v1=multi_sp_dat5[1:3,]

#fit, predict, score
multi_sp_dat5v1$output=map(multi_sp_dat5v1$splits, fitmod2_cast_score)

#make plots of species-level foreceasts

#DM

#very manual way of getting forecasts and test data at each window

obs1=unlist(multi_sp_dat5v1[[3]][[1]][[2]]$test_observations$DM) #get test data
fors1=as.matrix(multi_sp_dat5v1[[3]][[1]][[2]]$forecasts$DM) # get forecasts
time1=unlist(multi_sp_dat5v1[[3]][[1]][[2]]$test_times) # get test

obs2=unlist(multi_sp_dat5v1[[3]][[2]][[2]]$test_observations$DM) #get test data
fors2=as.matrix(multi_sp_dat5v1[[3]][[2]][[2]]$forecasts$DM) # get forecasts
time2=unlist(multi_sp_dat5v1[[3]][[2]][[2]]$test_times) # get test

obs3=unlist(multi_sp_dat5v1[[3]][[3]][[2]]$test_observations$DM) #get test data
fors3=as.matrix(multi_sp_dat5v1[[3]][[3]][[2]]$forecasts$DM) # get forecasts
time3=unlist(multi_sp_dat5v1[[3]][[3]][[2]]$test_times) # get test


par(mfrow=c(3,1))

matplot(t(fors1), type="l", col="grey", ylab="count", xaxt="n", main="DM")
axis(1,time1,at=c(1:12))
points(obs1, col="red", pch=21, cex=0.5, bg="red")
lines(obs1, col="black")

matplot(t(fors2), type="l", col="grey", ylab="count", xaxt="n", main="DM")
axis(1,time2,at=c(1:12))
points(obs2, col="red", pch=21, cex=0.5, bg="red")
lines(obs2, col="black")

matplot(t(fors3), type="l", col="grey", ylab="count", xaxt="n", main="DM")
axis(1,time3,at=c(1:12))
points(obs3, col="red", pch=21, cex=0.5, bg="red")
lines(obs3, col="black")

#DO

#very manual way of getting forecasts and test data at each window

obs1=unlist(multi_sp_dat5v1[[3]][[1]][[2]]$test_observations$DO) #get test data
fors1=as.matrix(multi_sp_dat5v1[[3]][[1]][[2]]$forecasts$DO) # get forecasts
time1=unlist(multi_sp_dat5v1[[3]][[1]][[2]]$test_times) # get test

obs2=unlist(multi_sp_dat5v1[[3]][[2]][[2]]$test_observations$DO) #get test data
fors2=as.matrix(multi_sp_dat5v1[[3]][[2]][[2]]$forecasts$DO) # get forecasts
time2=unlist(multi_sp_dat5v1[[3]][[2]][[2]]$test_times) # get test

obs3=unlist(multi_sp_dat5v1[[3]][[3]][[2]]$test_observations$DO) #get test data
fors3=as.matrix(multi_sp_dat5v1[[3]][[3]][[2]]$forecasts$DO) # get forecasts
time3=unlist(multi_sp_dat5v1[[3]][[3]][[2]]$test_times) # get test


par(mfrow=c(3,1))

matplot(t(fors1), type="l", col="grey", ylab="count", xaxt="n", main="DO")
axis(1,time1,at=c(1:12))
points(obs1, col="red", pch=21, cex=0.5, bg="red")
lines(obs1, col="black")

matplot(t(fors2), type="l", col="grey", ylab="count", xaxt="n", main="DO")
axis(1,time2,at=c(1:12))
points(obs2, col="red", pch=21, cex=0.5, bg="red")
lines(obs2, col="black")

matplot(t(fors3), type="l", col="grey", ylab="count", xaxt="n", main="DO")
axis(1,time3,at=c(1:12))
points(obs3, col="red", pch=21, cex=0.5, bg="red")
lines(obs3, col="black")

#PP

#very manual way of getting forecasts and test data at each window

obs1=unlist(multi_sp_dat5v1[[3]][[1]][[2]]$test_observations$PP) #get test data
fors1=as.matrix(multi_sp_dat5v1[[3]][[1]][[2]]$forecasts$PP) # get forecasts
time1=unlist(multi_sp_dat5v1[[3]][[1]][[2]]$test_times) # get test

obs2=unlist(multi_sp_dat5v1[[3]][[2]][[2]]$test_observations$PP) #get test data
fors2=as.matrix(multi_sp_dat5v1[[3]][[2]][[2]]$forecasts$PP) # get forecasts
time2=unlist(multi_sp_dat5v1[[3]][[2]][[2]]$test_times) # get test

obs3=unlist(multi_sp_dat5v1[[3]][[3]][[2]]$test_observations$PP) #get test data
fors3=as.matrix(multi_sp_dat5v1[[3]][[3]][[2]]$forecasts$PP) # get forecasts
time3=unlist(multi_sp_dat5v1[[3]][[3]][[2]]$test_times) # get test


par(mfrow=c(3,1))

matplot(t(fors1), type="l", col="grey", ylab="count", xaxt="n", main="PP")
axis(1,time1,at=c(1:12))
points(obs1, col="red", pch=21, cex=0.5, bg="red")
lines(obs1, col="black")

matplot(t(fors2), type="l", col="grey", ylab="count", xaxt="n", main="PP")
axis(1,time2,at=c(1:12))
points(obs2, col="red", pch=21, cex=0.5, bg="red")
lines(obs2, col="black")

matplot(t(fors3), type="l", col="grey", ylab="count", xaxt="n", main="PP")
axis(1,time3,at=c(1:12))
points(obs3, col="red", pch=21, cex=0.5, bg="red")
lines(obs3, col="black")

#PB

#very manual way of getting forecasts and test data at each window

obs1=unlist(multi_sp_dat5v1[[3]][[1]][[2]]$test_observations$PB) #get test data
fors1=as.matrix(multi_sp_dat5v1[[3]][[1]][[2]]$forecasts$PB) # get forecasts
time1=unlist(multi_sp_dat5v1[[3]][[1]][[2]]$test_times) # get test

obs2=unlist(multi_sp_dat5v1[[3]][[2]][[2]]$test_observations$PB) #get test data
fors2=as.matrix(multi_sp_dat5v1[[3]][[2]][[2]]$forecasts$PB) # get forecasts
time2=unlist(multi_sp_dat5v1[[3]][[2]][[2]]$test_times) # get test

obs3=unlist(multi_sp_dat5v1[[3]][[3]][[2]]$test_observations$PB) #get test data
fors3=as.matrix(multi_sp_dat5v1[[3]][[3]][[2]]$forecasts$PB) # get forecasts
time3=unlist(multi_sp_dat5v1[[3]][[3]][[2]]$test_times) # get test


par(mfrow=c(3,1))

matplot(t(fors1), type="l", col="grey", ylab="count", xaxt="n", main="PB")
axis(1,time1,at=c(1:12))
points(obs1, col="red", pch=21, cex=0.5, bg="red")
lines(obs1, col="black")

matplot(t(fors2), type="l", col="grey", ylab="count", xaxt="n", main="PB")
axis(1,time2,at=c(1:12))
points(obs2, col="red", pch=21, cex=0.5, bg="red")
lines(obs2, col="black")

matplot(t(fors3), type="l", col="grey", ylab="count", xaxt="n", main="PB")
axis(1,time3,at=c(1:12))
points(obs3, col="red", pch=21, cex=0.5, bg="red")
lines(obs3, col="black")

#RM

#very manual way of getting forecasts and test data at each window

obs1=unlist(multi_sp_dat5v1[[3]][[1]][[2]]$test_observations$RM) #get test data
fors1=as.matrix(multi_sp_dat5v1[[3]][[1]][[2]]$forecasts$RM) # get forecasts
time1=unlist(multi_sp_dat5v1[[3]][[1]][[2]]$test_times) # get test

obs2=unlist(multi_sp_dat5v1[[3]][[2]][[2]]$test_observations$RM) #get test data
fors2=as.matrix(multi_sp_dat5v1[[3]][[2]][[2]]$forecasts$RM) # get forecasts
time2=unlist(multi_sp_dat5v1[[3]][[2]][[2]]$test_times) # get test

obs3=unlist(multi_sp_dat5v1[[3]][[3]][[2]]$test_observations$RM) #get test data
fors3=as.matrix(multi_sp_dat5v1[[3]][[3]][[2]]$forecasts$RM) # get forecasts
time3=unlist(multi_sp_dat5v1[[3]][[3]][[2]]$test_times) # get test


par(mfrow=c(3,1))

matplot(t(fors1), type="l", col="grey", ylab="count", xaxt="n", main="RM")
axis(1,time1,at=c(1:12))
points(obs1, col="red", pch=21, cex=0.5, bg="red")
lines(obs1, col="black")

matplot(t(fors2), type="l", col="grey", ylab="count", xaxt="n", main="RM")
axis(1,time2,at=c(1:12))
points(obs2, col="red", pch=21, cex=0.5, bg="red")
lines(obs2, col="black")

matplot(t(fors3), type="l", col="grey", ylab="count", xaxt="n", main="RM")
axis(1,time3,at=c(1:12))
points(obs3, col="red", pch=21, cex=0.5, bg="red")
lines(obs3, col="black")

#PE

#very manual way of getting forecasts and test data at each window

obs1=unlist(multi_sp_dat5v1[[3]][[1]][[2]]$test_observations$PE) #get test data
fors1=as.matrix(multi_sp_dat5v1[[3]][[1]][[2]]$forecasts$PE) # get forecasts
time1=unlist(multi_sp_dat5v1[[3]][[1]][[2]]$test_times) # get test

obs2=unlist(multi_sp_dat5v1[[3]][[2]][[2]]$test_observations$PE) #get test data
fors2=as.matrix(multi_sp_dat5v1[[3]][[2]][[2]]$forecasts$PE) # get forecasts
time2=unlist(multi_sp_dat5v1[[3]][[2]][[2]]$test_times) # get test

obs3=unlist(multi_sp_dat5v1[[3]][[3]][[2]]$test_observations$PE) #get test data
fors3=as.matrix(multi_sp_dat5v1[[3]][[3]][[2]]$forecasts$PE) # get forecasts
time3=unlist(multi_sp_dat5v1[[3]][[3]][[2]]$test_times) # get test


par(mfrow=c(3,1))

matplot(t(fors1), type="l", col="grey", ylab="count", xaxt="n", main="PE")
axis(1,time1,at=c(1:12))
points(obs1, col="red", pch=21, cex=0.5, bg="red")
lines(obs1, col="black")

matplot(t(fors2), type="l", col="grey", ylab="count", xaxt="n", main="PE")
axis(1,time2,at=c(1:12))
points(obs2, col="red", pch=21, cex=0.5, bg="red")
lines(obs2, col="black")

matplot(t(fors3), type="l", col="grey", ylab="count", xaxt="n", main="PE")
axis(1,time3,at=c(1:12))
points(obs3, col="red", pch=21, cex=0.5, bg="red")
lines(obs3, col="black")


