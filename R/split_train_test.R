library(dplyr)

# Cut one training/testing window out of the prepared data, keeping vector and
# matrix fields indexed together. `species` restricts the split to a species
# group; species below 30% occupancy in the training window are then dropped,
# so the species actually fit are returned as `species_list`.
split_train_test <- function(
  data_all,
  gap,
  train_start,
  train_end,
  test_start,
  test_end,
  species = NULL
) {
  species <- species %||% levels(data_all$series)
  species_list <- data.frame(
    newmoonnumber = data_all$newmoonnumber,
    series = data_all$series,
    y = data_all$y
  ) |>
    filter(
      newmoonnumber >= train_start,
      newmoonnumber <= (train_end),
      series %in% species
    ) |>
    group_by(newmoonnumber, series) |>
    summarise(abundance = sum(y, na.rm = TRUE), .groups = 'drop') |>
    group_by(series) |>
    summarise(occupancy = sum(abundance > 0) / n(), .groups = 'drop') |>
    filter(occupancy >= 0.30) |>
    pull(series) |>
    droplevels()

  train_inds <- which(
    data_all$newmoonnumber >= train_start &
      data_all$newmoonnumber <= (train_end) &
      data_all$series %in% species_list
  )
  data_train <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][train_inds, ]
    } else {
      data_all[[x]][train_inds]
    }
  })

  names(data_train) <- names(data_all)
  data_train$series <- droplevels(data_train$series)

  test_inds <- which(
    data_all$newmoonnumber >= (test_start) &
      data_all$newmoonnumber <= test_end &
      data_all$series %in% species_list
  )
  data_test <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][test_inds, ]
    } else {
      data_all[[x]][test_inds]
    }
  })

  names(data_test) <- names(data_all)
  data_test$series <- droplevels(data_test$series)

  return(list(
    train = data_train,
    test = data_test,
    species_list = species_list
  ))
}
