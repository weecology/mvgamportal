baseline <- as.data.frame(baseline_score) %>%
  mutate(newmoonnumber = test_start_newmoonnumber + DM.eval_horizon - 1) %>%
  select(-contains("score_type")) %>%
  tidyr::pivot_longer(cols = !c("test_start_newmoonnumber","newmoonnumber"),
                      names_to = c("species", "type"),
                      names_sep = "\\.") %>%
  tidyr::pivot_wider(names_from = "type",
                     values_from = "value") %>%
  select(test_start_newmoonnumber, newmoonnumber, species, baseline_score=score)

ar <- as.data.frame(ar_score) %>%
  mutate(newmoonnumber = test_start_newmoonnumber + DM.eval_horizon - 1) %>%
  select(-contains("score_type")) %>%
  tidyr::pivot_longer(cols = !c("test_start_newmoonnumber","newmoonnumber"),
                      names_to = c("species", "type"),
                      names_sep = "\\.") %>%
  tidyr::pivot_wider(names_from = "type",
                     values_from = "value") %>%
  left_join(baseline, by = join_by(test_start_newmoonnumber, newmoonnumber, species)) %>%
  mutate(skill_score = 1 - score/baseline_score,
         model = "AR") 

gam_ar <- as.data.frame(gam_ar_score) %>%
  mutate(newmoonnumber = test_start_newmoonnumber + DM.eval_horizon - 1) %>%
  select(-contains("score_type")) %>%
  tidyr::pivot_longer(cols = !c("test_start_newmoonnumber","newmoonnumber"),
                      names_to = c("species", "type"),
                      names_sep = "\\.") %>%
  tidyr::pivot_wider(names_from = "type",
                     values_from = "value") %>%
  left_join(baseline, by = join_by(test_start_newmoonnumber, newmoonnumber, species)) %>%
  mutate(skill_score = 1 - score/baseline_score,
         model = "GAM_AR") 

gam_var <- as.data.frame(gam_var_score) %>%
  mutate(newmoonnumber = test_start_newmoonnumber + DM.eval_horizon - 1) %>%
  select(-contains("score_type")) %>%
  tidyr::pivot_longer(cols = !c("test_start_newmoonnumber","newmoonnumber"),
                      names_to = c("species", "type"),
                      names_sep = "\\.") %>%
  tidyr::pivot_wider(names_from = "type",
                     values_from = "value") %>%
  left_join(baseline, by = join_by(test_start_newmoonnumber, newmoonnumber, species)) %>%
  mutate(skill_score = 1 - score/baseline_score,
         model = "GAM_VAR") 
