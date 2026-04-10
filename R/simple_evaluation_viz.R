library(dplyr)
library(ggplot2)

simple_dfs <- lapply(simple_scores, data.frame, stringsAsFactors = FALSE)
simple <- dplyr::bind_rows(simple_dfs) %>%
  mutate(newmoonnumber = test_start_newmoonnumber + DM.eval_horizon - 1) %>%
  select(-contains("score_type")) %>%
  tidyr::pivot_longer(cols = !c("test_start_newmoonnumber","newmoonnumber","species_list","rhat","prhat_high","n_divergences"),
                      names_to = c("species", "type"),
                      names_sep = "\\.") %>%
  tidyr::pivot_wider(names_from = "type",
                     values_from = "value") %>%
  left_join(baseline, by = join_by(test_start_newmoonnumber, newmoonnumber, species)) %>%
  mutate(skill_score = 1 - score/baseline_score,
         model = "Simple") 

env_dfs <- lapply(env_distances, data.frame, stringsAsFactors = FALSE)
environment <- bind_rows(env_dfs)

skill_scores <- simple
overall_skill_scores = skill_scores %>%
  filter(species=="all_series")
species_skill_scores = skill_scores %>%
  filter(species!="all_series")
dm_skill_scores = skill_scores[skill_scores$species=="DM",]
pp_skill_scores = skill_scores[skill_scores$species=="PP",]
pb_skill_scores = skill_scores[skill_scores$species=="PB",]
do_skill_scores = skill_scores[skill_scores$species=="DO",]

ggplot(data=overall_skill_scores, aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_line() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  #  ylim(-1,1) +
  theme_minimal() +
  theme(legend.position = "none") 

ggplot(data=dm_skill_scores, aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_line() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(data=pp_skill_scores, aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_line() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(data=pb_skill_scores, aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_line() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(data=do_skill_scores, aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_line() +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(data=subset(species_skill_scores,subset = model=="Simple"), aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point() +
  facet_wrap(~species, ncol = 3, scales = "free_y") +
  #  ylim(-1,1) +
  theme_minimal() +
  theme(legend.position = "none")

p1 <- ggplot(data=subset(overall_skill_scores,subset = model=="Simple")) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  # ylim(-1,1) +
  theme_minimal() +
  theme(legend.position = "none") 
p2 <- ggplot(data=subset(overall_skill_scores,subset = model=="Simple" & eval_horizon==1)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=newmoonnumber, y=rhat), color="red2") +
  ylim(1,1.02) +
  theme_minimal() +
  theme(legend.position = "none") 
p3 <- ggplot(data=subset(overall_skill_scores,subset = model=="Simple" & eval_horizon==1)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=newmoonnumber, y=n_divergences), color="red2") +
  theme_minimal() +
  theme(legend.position = "none")
p4 <- ggplot(environment) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=test_start_newmoonnumber, y=ndvi), color="green4") +
  geom_point(aes(x=test_start_newmoonnumber, y=mintemp), color="red4") +
  ylab("dist") +
  theme_minimal() +
  theme(legend.position = "none")

p1 / p2 / p3 / p4
