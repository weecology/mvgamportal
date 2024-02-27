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

#function to fit AR1 model, predict, and score

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
  
  pred_score= score(preds)
  
  return(list(model, preds, pred_score))
  
}

# access data and perform iterative forecasting
#run on only a subset; for checking

dmdat20=dmdat%>%
  sliding_index(
    newmoonnumber,
    lookback=240,
    assess_stop=12,
    complete = TRUE
  )%>%
  head(n=5)%>%
  mutate(output=map(splits, fitmod1_cast_score))

