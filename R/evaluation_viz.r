library(dplyr)
library(ggplot2)
library(dplyr)
library(ggpp)

# scores <- readRDS("scores.rds")
# env_distances <- readRDS("env_distances.rds")

env_dfs <- lapply(env_distances, data.frame, stringsAsFactors = FALSE)
environment <- bind_rows(env_dfs)

skill_scores <- scores %>% filter(model != "BASELINE")
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

ggplot(data=subset(species_skill_scores,subset = model=="AR"), aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point() +
  facet_wrap(~species, ncol = 3, scales = "free_y") +
#  ylim(-1,1) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(data=subset(species_skill_scores,subset = model=="GAM_AR"), aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point() +
  facet_wrap(~species, ncol = 3, scales = "free_y") +
#  ylim(-1,1) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(data=subset(species_skill_scores,subset = model=="GAM_VAR"), aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point() +
  facet_wrap(~species, ncol = 3, scales = "free_y") +
#  ylim(-1,1) +
  theme_minimal() +
  theme(legend.position = "none")

p1 <- ggplot(data=subset(overall_skill_scores,subset = model=="AR")) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
# ylim(-1,1) +
  theme_minimal() +
  theme(legend.position = "none")
p2 <- ggplot(data=subset(overall_skill_scores,subset = model=="AR" & eval_horizon==1)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=newmoonnumber, y=rhat), color="red2") +
  ylim(1,1.02) +
  theme_minimal() +
  theme(legend.position = "none")
p3 <- ggplot(data=subset(overall_skill_scores,subset = model=="AR" & eval_horizon==1)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=newmoonnumber, y=n_divergences), color="red2") +
  theme_minimal() +
  theme(legend.position = "none")

p1 / p2 / p3

p1 <- ggplot(data=subset(overall_skill_scores,subset = model=="GAM_AR")) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  # ylim(-1,1) +
  theme_minimal() +
  theme(legend.position = "none")
p2 <- ggplot(data=subset(overall_skill_scores,subset = model=="GAM_AR" & eval_horizon==1)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=newmoonnumber, y=rhat), color="red2") +
  ylim(1,1.02) +
  theme_minimal() +
  theme(legend.position = "none")
p3 <- ggplot(data=subset(overall_skill_scores,subset = model=="GAM_AR" & eval_horizon==1)) +
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

p1 <- ggplot(data=subset(overall_skill_scores,subset = model=="GAM_VAR")) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  # ylim(-1,1) +
  theme_minimal() +
  theme(legend.position = "none")
p2 <- ggplot(data=subset(overall_skill_scores,subset = model=="GAM_VAR" & eval_horizon==1)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=newmoonnumber, y=rhat), color="red2") +
  ylim(1,1.02) +
  theme_minimal() +
  theme(legend.position = "none")
p3 <- ggplot(data=subset(overall_skill_scores,subset = model=="GAM_VAR" & eval_horizon==1)) +
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

p1 <- ggplot(data=subset(overall_skill_scores,subset = model=="SIMPLE")) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  theme_minimal() +
  theme(legend.position = "none")
p2 <- ggplot(data=subset(overall_skill_scores,subset = model=="SIMPLE" & eval_horizon==1)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point(aes(x=newmoonnumber, y=rhat), color="red2") +
  ylim(1,1.02) +
  theme_minimal() +
  theme(legend.position = "none")
p3 <- ggplot(data=subset(overall_skill_scores,subset = model=="SIMPLE" & eval_horizon==1)) +
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

ggplot(data=subset(species_skill_scores,subset = model=="SIMPLE"), aes(x=newmoonnumber, y=skill_score, color=eval_horizon)) +
  geom_rect(aes(xmin=140, xmax=230, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=263, xmax=278, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_rect(aes(xmin=396, xmax=411, ymin=0, ymax=Inf), fill="lightgrey",alpha=0.2, color=NA) +
  geom_point() +
  facet_wrap(~species, ncol = 3, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none")

skill_sign_overall <- overall_skill_scores %>%
  filter(!is.na(skill_score), is.finite(skill_score)) %>%
  mutate(skill_sign = ifelse(skill_score >= 0, "positive", "negative"))

ggplot(skill_sign_overall, aes(x=model, fill=skill_sign)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(positive = "steelblue", negative = "firebrick")) +
  labs(x = "Model", y = "Count", fill = "Skill score") +
  theme_minimal()

skill_sign_species <- species_skill_scores %>%
  filter(!is.na(skill_score), is.finite(skill_score)) %>%
  mutate(skill_sign = ifelse(skill_score >= 0, "positive", "negative"))

ggplot(skill_sign_species, aes(x=model, fill=skill_sign)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c(positive = "steelblue", negative = "firebrick")) +
  facet_wrap(~species, scales = "free_y") +
  labs(x = "Model", y = "Count", fill = "Skill score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

skill_binwidth <- 0.1
skill_clamp_value <- -1 - skill_binwidth / 2
inset_npcx <- skill_binwidth / (2 + skill_binwidth) + 0.05

skill_hist_overall <- overall_skill_scores %>%
  filter(!is.na(skill_score), is.finite(skill_score), skill_score <= 1) %>%
  mutate(skill_score_clamped = pmax(skill_score, skill_clamp_value),
         pooled = skill_score < -1)

get_hist_ymax <- function(x, binwidth, boundary = 0) {
  left  <- floor((min(x)  - boundary) / binwidth) * binwidth + boundary
  right <- ceiling((max(x) - boundary) / binwidth) * binwidth + boundary
  breaks <- seq(left, right + binwidth, by = binwidth)
  max(hist(x, breaks = breaks, plot = FALSE)$counts)
}

make_tail_inset <- function(d, ymax) {
  ggplot(d, aes(x = skill_score)) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "lightgrey", alpha = 0.6) +
    geom_histogram(binwidth = skill_binwidth, boundary = 0, fill = "steelblue4", color = "white") +
    coord_cartesian(ylim = c(0, ymax)) +
    theme_minimal(base_size = 6) +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 5),
          plot.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3))
}

hist_ymax_overall <- skill_hist_overall %>%
  group_by(model) %>%
  summarise(max_count = get_hist_ymax(skill_score_clamped, skill_binwidth), .groups = "drop")

tail_insets_overall <- skill_hist_overall %>%
  group_by(model) %>%
  group_modify(~tibble(tail_data = list(.x))) %>%
  left_join(hist_ymax_overall, by = "model") %>%
  mutate(inset = Map(make_tail_inset, tail_data, max_count), npcx = inset_npcx, npcy = 0.98) %>%
  select(-tail_data, -max_count)

ggplot(skill_hist_overall, aes(x=skill_score_clamped, fill=pooled)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "lightgrey", alpha = 0.6) +
  geom_histogram(binwidth = skill_binwidth, boundary = 0, color = "white") +
  scale_fill_manual(values = c(`FALSE` = "steelblue", `TRUE` = "steelblue4"),
                    labels = c(`FALSE` = "binned", `TRUE` = "pooled (< -1)")) +
  geom_plot_npc(data = tail_insets_overall,
                aes(npcx = npcx, npcy = npcy, label = inset),
                hjust = 0, vjust = 1, vp.width = 0.3, vp.height = 0.35) +
  facet_wrap(~model, scales = "free_y") +
  coord_cartesian(xlim = c(-1 - skill_binwidth, 1), expand = FALSE) +
  labs(x = "Skill score (values < -1 pooled)", y = "Count", fill = NULL) +
  theme_minimal()

skill_hist_species <- species_skill_scores %>%
  filter(!is.na(skill_score), is.finite(skill_score), skill_score <= 1) %>%
  mutate(skill_score_clamped = pmax(skill_score, skill_clamp_value),
         pooled = skill_score < -1)

hist_ymax_species <- skill_hist_species %>%
  group_by(species, model) %>%
  summarise(max_count = get_hist_ymax(skill_score_clamped, skill_binwidth), .groups = "drop")

tail_insets_species <- skill_hist_species %>%
  group_by(species, model) %>%
  group_modify(~tibble(tail_data = list(.x))) %>%
  left_join(hist_ymax_species, by = c("species", "model")) %>%
  mutate(inset = Map(make_tail_inset, tail_data, max_count), npcx = inset_npcx, npcy = 0.98) %>%
  select(-tail_data, -max_count)

ggplot(skill_hist_species, aes(x=skill_score_clamped, fill=pooled)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "lightgrey", alpha = 0.6) +
  geom_histogram(binwidth = skill_binwidth, boundary = 0, color = "white") +
  scale_fill_manual(values = c(`FALSE` = "steelblue", `TRUE` = "steelblue4"),
                    labels = c(`FALSE` = "binned", `TRUE` = "pooled (< -1)")) +
  geom_plot_npc(data = tail_insets_species,
                aes(npcx = npcx, npcy = npcy, label = inset),
                hjust = 0, vjust = 1, vp.width = 0.3, vp.height = 0.35) +
  facet_grid(species ~ model, scales = "free_y") +
  coord_cartesian(xlim = c(-1 - skill_binwidth, 1), expand = FALSE) +
  labs(x = "Skill score (values < -1 pooled)", y = "Count", fill = NULL) +
  theme_minimal()
