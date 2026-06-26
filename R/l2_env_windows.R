library(dplyr)
library(dad)
library(ggplot2)

data_all <- readRDS("data_heteromyid.rds")

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
                        data_all$newmoonnumber <= (train_end) & 
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
test_starts = seq(from = 200, to = 500, by = 20)
train_starts = test_starts - train_win_width

calc_window_l2 <- function(train_start) {
  train_end <- train_start + train_win_width - 1
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
  
  env_train = data.frame(ndvi=data_train$ndvi, mintemp = data_train$meantemp_lag_1)
  env_test = data.frame(ndvi=data_test$ndvi, mintemp = data_test$meantemp_lag_1)
  l2_train_downsampled <- sapply(1:100, function(i) 
    distl2d(env_train, env_train[sample(nrow(env_train), 12, replace=FALSE), ], method="kern"))
  env_distance <- data.frame(l2_test = distl2d(env_train, env_test, method="kern"),
    l2_train_downsampled = l2_train_downsampled)
  env_distance$test_start_newmoonnumber <- test_start
  env_distance$species_list <- paste(data_split$species_list,collapse="_")
  
  env_distance
}

l2_list <- lapply(train_starts, calc_window_l2)
l2_results <- do.call(rbind, l2_list)

ggplot(data = l2_results, aes(x = test_start_newmoonnumber, y = l2_train_downsampled, group = test_start_newmoonnumber)) +
  geom_boxplot() +
  geom_point(aes(x = test_start_newmoonnumber, y = l2_test, group = test_start_newmoonnumber),color= "blue",size=4) +
  ylim(0,.5) +
  theme_minimal()
