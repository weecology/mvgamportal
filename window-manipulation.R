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
  arrange(time)

# Create sliding windows

#parameters

window_size <- 6  # size of the sliding window
step_size <- 12    # step size for the sliding window

# Initialize a list to store the results
dm_slide <- list()

for (nyr in seq(5,15, by=1)) {

  window_size <- window_size  # size of the sliding window

  for (start in seq(1, nrow(dmdat) - window_size + 1, by = step_size)) {

    end <- start + window_size - 1
    dm_slide[[length(dm_slide) + 1]] <- dmdat[start:end, ]

  }
}
