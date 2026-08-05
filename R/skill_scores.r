library(dplyr)
library(purrr)
library(tidyr)

tidy_score <- function(score_df, model_name, separately = FALSE) {
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
      cols = !any_of(c(
        "test_start_newmoonnumber",
        "newmoonnumber",
        "species_set",
        "species_list",
        "rhat",
        "prhat_high",
        "n_divergences"
      )),
      names_to = c("species", "type"),
      names_sep = "\\."
    ) |>
    pivot_wider(names_from = "type", values_from = "value") |>
    # score() adds an all_series element summing across species,
    # for separately fit species this is confusing because it
    # is just a copy of the single species score. When comparing
    # single species models to multi-species models we want to do
    # this as the species level, so remove all_series when fitting
    # separately to avoid confusion/mistakes
    filter(!separately | species != "all_series") |>
    mutate(model = model_name)
}

get_skill_scores <- function(score_dfs, separately) {
  scores <- bind_rows(purrr::pmap(
    list(score_dfs, names(score_dfs), separately),
    tidy_score
  ))

  # Matched within species_set so each model is compared against the baseline
  # fit to the same group of species
  baseline_ref <- scores |>
    filter(model == "BASELINE") |>
    # Explicit namespacing required - something is loading MASS & overwriting
    dplyr::select(
      test_start_newmoonnumber,
      newmoonnumber,
      species_set,
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
        species_set,
        species,
        eval_horizon
      )
    ) |>
    mutate(skill_score = 1 - score / baseline_score)
}
