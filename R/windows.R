test_win_width <- 13

get_window <- function(train_start, train_win_width) {
  train_end <- train_start + train_win_width - 1
  list(
    train_start = train_start,
    train_end = train_end,
    test_start = train_end + 1,
    test_end = train_end + test_win_width
  )
}

get_train_starts <- function(newmoons, analysis_config) {
  train_win_width <- analysis_config$train_win_width
  if (is.null(analysis_config$test_starts)) {
    return(min(newmoons):(max(newmoons) - train_win_width - test_win_width + 1))
  }
  test_starts <- seq(
    from = analysis_config$test_starts$from,
    to = analysis_config$test_starts$to,
    by = analysis_config$test_starts$by %||% 1
  )
  test_starts - train_win_width
}

# Report out on failed windows
handle_failed_windows <- function(results, train_starts, task) {
  errored <- purrr::map_lgl(results, ~ !is.null(.x$error))
  if (any(errored)) {
    purrr::walk2(
      train_starts[errored],
      purrr::map(results[errored], "error"),
      \(train_start, error) message(glue::glue(
        "{task} failed for train_start {train_start}: {conditionMessage(error)}"
      ))
    )
    warning(
      glue::glue(
        "{task} failed for {sum(errored)} window(s): train_starts ",
        "{paste(train_starts[errored], collapse = ', ')}"
      ),
      call. = FALSE
    )
  }
  purrr::map(results, "result")
}
