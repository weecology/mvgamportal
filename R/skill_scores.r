library(dplyr)
library(purrr)
library(tidyr)

tidy_score <- function(score_df, model_name) {
  as.data.frame(score_df) |>
    # Any species' eval_horizon column works; which species are present varies
    # by window, so grab the first one rather than naming it
    mutate(
      newmoonnumber = test_start_newmoonnumber +
        pick(ends_with(".eval_horizon"))[[1]] -
        1
    ) |>
    # Explicit namespacing required - something is loading MASS & overwriting
    dplyr::select(-contains("score_type")) |>
    pivot_longer(
      cols = !c(
        "test_start_newmoonnumber",
        "newmoonnumber",
        "species_list",
        "rhat",
        "prhat_high",
        "n_divergences"
      ),
      names_to = c("species", "type"),
      names_sep = "\\."
    ) |>
    pivot_wider(names_from = "type", values_from = "value") |>
    mutate(model = model_name)
}

get_skill_scores <- function(score_dfs) {
  scores <- bind_rows(purrr::imap(score_dfs, \(df, name) tidy_score(df, name)))

  baseline_ref <- scores |>
    filter(model == "BASELINE") |>
    # Explicit namespacing required - something is loading MASS & overwriting
    dplyr::select(
      test_start_newmoonnumber,
      newmoonnumber,
      species,
      eval_horizon,
      baseline_score = score
    )

  scores |>
    left_join(
      baseline_ref,
      by = join_by(
        test_start_newmoonnumber,
        newmoonnumber,
        species,
        eval_horizon
      )
    ) |>
    mutate(skill_score = 1 - score / baseline_score)
}
