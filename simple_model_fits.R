library(dplyr)
library(glue)
library(statip)
library(mvgam)

### Data

data_all <- readRDS("data_heteromyid.rds")
data_all$mintemp_lag_1 <- data_all$mintemp[,'mintemp_lag_1']

### Split train test
split_train_test <- function(data_all, gap, train_start, train_end, test_start, test_end) {
  species_list <- data.frame(newmoonnumber=data_all$newmoonnumber,series=data_all$series,y=data_all$y) |>
    filter(newmoonnumber >= train_start, newmoonnumber <= (train_end)) |>
    group_by(newmoonnumber, series) |>
    summarise(abundance = sum(y, na.rm = TRUE), .groups = 'drop') |>
    group_by(series) |>
    summarise(occupancy = sum(abundance > 0) / n(), .groups = 'drop') |>
    filter(occupancy >= 0.30) |>
    pull(series) |>
    droplevels()
  
  train_inds <- which(data_all$newmoonnumber >= train_start &
                        data_all$newmoonnumber <= (train_end) &     # Old +1 here keeps the first observation after the gap to serve as the initial condition, also cut +gap since no gaps
                        data_all$series %in% species_list)
  data_train <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][train_inds, ]
    } else {
      data_all[[x]][train_inds]
    }
  })
  
  names(data_train) <- names(data_all)
  data_train$series <- droplevels(data_train$series)
  # gap_inds <- which(data_train$newmoonnumber >= (train_end + 1) &
  #   data_train$newmoonnumber <= (train_end + gap))
  # data_train$y[gap_inds] <- NA
  
  test_inds <- which(data_all$newmoonnumber >= (test_start) & # Old +1 here was to skip the initial condition point
                       data_all$newmoonnumber <= test_end &
                       data_all$series %in% species_list)
  data_test <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][test_inds, ]
    } else {
      data_all[[x]][test_inds]
    }
  })
  
  names(data_test) <- names(data_all)
  data_test$series <- droplevels(data_test$series)
  
  return(list(train = data_train, test = data_test, species_list=species_list))
}

newmoon_min <- min(data_all$newmoonnumber)
newmoon_max <- max(data_all$newmoonnumber)
train_win_width <- 60
train_starts <- newmoon_min:(newmoon_max - train_win_width - 12 + 1)

### Loop setup
simple_scores <- vector(mode = "list", length = length(train_starts))
simple_summaries <- vector(mode = "list", length = length(train_starts))
env_distances <- vector(mode = "list", length = length(train_starts))

targets = c(158,180,200,220,240,248,260,280,285,300,320,340,360,380,400,420,428,440,460,480)
initial = targets - 60

### Fit all windows

for (i in initial) {
  target = i + 60
  train_start <- train_starts[i - 95]
  print(glue("Starting training: {i}"))
  train_end <- train_start + 60 - 1
  test_start <- train_end + 1
  test_end <- test_start + 12 - 1
  data_split <- split_train_test(
    data_all,
    gap = 0,
    train_start = train_start,
    train_end = train_end,
    test_start = test_start,
    test_end = test_end
  )
  data_train <- data_split$train
  data_test <- data_split$test

  baseline_model <- mvgam(
    formula = y ~ -1 + series,
    data = data_train,
    newdata = data_test,
    family = poisson()         
  )  

  simple_model <- mvgam(
    formula = y ~ -1,
    trend_formula = ~ -1 +
     s(mintemp_lag_1, by = trend, k = 4) +
     s(ndvi, k = 4) +
     s(winter_ndvi, by = trend, k = 4) +
     s(summer_ndvi, by = trend, k = 4),
   data = data_train,
   newdata = data_test,
   family = poisson(),
   trend_model = AR()
 )
  
 baseline_score <- score(forecast(baseline_model), score = "crps")
 baseline_score$test_start_newmoonnumber <- test_start
 simple_score <- score(forecast(simple_model), score = "crps")
 simple_score$test_start_newmoonnumber <- test_start
 simple_score$species_list <- paste(data_split$species_list,collapse="_")
 simple_score$rhat <- mean(rhat(simple_model),na.rm=TRUE)
 simple_score$prhat_high <- mean(rhat(simple_model)>1.05,na.rm=TRUE)
 simple_score$n_divergences <- sum(sapply(
   rstan::get_sampler_params(simple_model$model_output, inc_warmup = FALSE), 
   function(x) sum(x[, 'divergent__'])))
 simple_scores[[i]] <- simple_score
 simple_summary <- summary(simple_model)
 simple_summary$test_start_newmoonnumber <- test_start
 simple_summary$species_list <- paste(data_split$species_list,collapse="_")
 simple_summaries[[i]] <- simple_summary

 env_distance <- data.frame(ndvi=hellinger(data_train$ndvi, data_test$ndvi),
                            mintemp=hellinger(data_train$mintemp, data_test$mintemp))
 env_distance$test_start_newmoonnumber <- test_start
 env_distance$species_list <- paste(data_split$species_list,collapse="_")
 env_distances[[i]] <- env_distance
 
 baseline <- as.data.frame(baseline_score) %>%
   mutate(newmoonnumber = test_start_newmoonnumber + DM.eval_horizon - 1) %>%
   select(-contains("score_type")) %>%
   tidyr::pivot_longer(cols = !c("test_start_newmoonnumber","newmoonnumber"),
                       names_to = c("species", "type"),
                       names_sep = "\\.") %>%
   tidyr::pivot_wider(names_from = "type",
                      values_from = "value") %>%
   select(test_start_newmoonnumber, newmoonnumber, species, score,eval_horizon)

 simple <- as.data.frame(simple_score) %>%
   mutate(newmoonnumber = test_start_newmoonnumber + DM.eval_horizon - 1) %>%
   select(-contains("score_type"),-c("rhat","prhat_high","n_divergences")) %>%
   tidyr::pivot_longer(cols = !c("test_start_newmoonnumber","newmoonnumber","species_list"),
                       names_to = c("species", "type"),
                       names_sep = "\\.") %>%
   tidyr::pivot_wider(names_from = "type",
                      values_from = "value") %>%
   left_join(baseline, by = join_by(test_start_newmoonnumber, newmoonnumber, species, eval_horizon)) %>%
   rename(score = score.x, baseline_score = score.y) %>%
   mutate(skill_score = 1 - score/baseline_score,
          model = "Simple") 
 
 source("~/mvgamportal/R/simple_forcastplots.R", echo = TRUE)

}
