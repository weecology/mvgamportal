library(dplyr)
library(ggplot2)

baseline_dfs <- lapply(baseline_scores, data.frame, stringsAsFactors = FALSE)
baseline <- bind_rows(baseline_dfs) %>%
  mutate(newmoonnumber = test_start_newmoonnumber + DM.eval_horizon - 1) %>%
  select(-contains("score_type")) %>%
  tidyr::pivot_longer(cols = !c("test_start_newmoonnumber","newmoonnumber"),
                      names_to = c("species", "type"),
                      names_sep = "\\.") %>%
  tidyr::pivot_wider(names_from = "type",
                     values_from = "value") %>%
  select(test_start_newmoonnumber, newmoonnumber, species, baseline_score=score)

ar_dfs <- lapply(ar_scores, data.frame, stringsAsFactors = FALSE)
ar <- dplyr::bind_rows(ar_dfs) %>%
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

gam_ar_dfs <- lapply(gam_ar_scores, data.frame, stringsAsFactors = FALSE)
gam_ar <- dplyr::bind_rows(gam_ar_dfs) %>%
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

gam_var_dfs <- lapply(gam_var_scores, data.frame, stringsAsFactors = FALSE)
gam_var <- dplyr::bind_rows(gam_var_dfs) %>%
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

skill_scores <- rbind(ar,gam_ar,gam_var)
overall_skill_scores = skill_scores %>%
                       filter(species=="all_series")
species_skill_scores = skill_scores %>%
                       filter(species!="all_series")
dm_skill_scores = skill_scores[skill_scores$species=="DM",]
pp_skill_scores = skill_scores[skill_scores$species=="PP",]
pb_skill_scores = skill_scores[skill_scores$species=="PB",]
do_skill_scores = skill_scores[skill_scores$species=="DO",]

ggplot(data=overall_skill_scores, aes(x=newmoonnumber, y=skill_score, color=model)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_line() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none") 

ggplot(data=dm_skill_scores, aes(x=newmoonnumber, y=skill_score, color=model)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_line() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(data=pp_skill_scores, aes(x=newmoonnumber, y=skill_score, color=model)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_line() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(data=pb_skill_scores, aes(x=newmoonnumber, y=skill_score, color=model)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_line() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(data=do_skill_scores, aes(x=newmoonnumber, y=score, color=model)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_line() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(data=species_skill_scores, aes(x=newmoonnumber, y=score, color=model)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_line() +
  facet_wrap(~species, ncol = 4, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none")
