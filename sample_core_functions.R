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


