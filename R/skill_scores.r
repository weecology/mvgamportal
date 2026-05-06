library(dplyr)
library(tidyr)

tidy_score <- function(score_df, model_name) {
  as.data.frame(score_df) %>%
    mutate(newmoonnumber = test_start_newmoonnumber + DM.eval_horizon - 1) %>%
    select(-contains("score_type")) %>%
    pivot_longer(cols = !c("test_start_newmoonnumber", "newmoonnumber",
                           "species_list", "rhat", "prhat_high", "n_divergences"),
                 names_to = c("species", "type"),
                 names_sep = "\\.") %>%
    pivot_wider(names_from = "type", values_from = "value") %>%
    mutate(model = model_name)
}

scores <- bind_rows(
  tidy_score(baseline_score, "BASELINE"),
  tidy_score(ar_score,       "AR"),
  tidy_score(gam_ar_score,   "GAM_AR"),
  tidy_score(gam_var_score,  "GAM_VAR"),
  tidy_score(simple_score,   "SIMPLE")
)

baseline_ref <- scores %>%
  filter(model == "BASELINE") %>%
  select(test_start_newmoonnumber, newmoonnumber, species, eval_horizon,
         baseline_score = score)

scores <- scores %>%
  left_join(baseline_ref,
            by = join_by(test_start_newmoonnumber, newmoonnumber, species, eval_horizon)) %>%
  mutate(skill_score = 1 - score / baseline_score)
