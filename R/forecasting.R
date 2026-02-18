# from Feb 11th 
ar_outregime <- readRDS("ar_outregime_output412.rds")
ar_outregime
library(mvgam)
plot(ar_outregime, series = 1)

ar_out_forecast = forecast(ar_outregime)
plot(ar_out_forecast, series = 1)
