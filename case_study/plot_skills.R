library(dplyr)
library(ggplot2)
library(patchwork)

skills = read.csv("model_skills.csv", stringsAsFactors = FALSE)
skills$train_batch = as.factor(skills$train_batch)

plot_skills_by_start = function(data, sp){
p1 =  ggplot(data = filter(data, species == sp & model_type == "gam_var_scores"), aes(horizon, skill, color=train_batch)) + 
  geom_point() + geom_line() + ggtitle("GAM VAR")
p2 =  ggplot(data = filter(data, species == sp & model_type == "gam_ar_scores"), aes(horizon, skill, color=train_batch)) + 
  geom_point() + geom_line() + ggtitle("GAM AR")
p3 = ggplot(data = filter(data, species == sp & model_type == "ar_scores"), aes(horizon, skill, color=train_batch)) + 
  geom_point() + geom_line() + ggtitle("AR")
# p4 = ggplot(data = filter(data, species == sp & model_type == "linear_mintemp_scores"), aes(horizon, skill, color=train_batch)) + 
#   geom_point() + geom_line() + ggtitle("Linear Mintemp")

patchwork = (p1 | p2) / (p3)
patchwork + plot_annotation(title = sp)
}

plot_skills_by_start(skills, "DM")


data = skills |> group_by(model_type, species, horizon) |> summarise(mean_skill = mean(skill), 
                                                                   sd_value = sd(skill))
# plot_skills_by_species = 
sp = "DM"
ggplot(data = filter(data, species == sp), aes(horizon, mean_skill, color=model_type)) + 
  geom_point() + geom_line() + ggtitle(sp)
p1  

data = skills
sp = "DM"
species_gamvar = data |> filter(species == sp & model_type == "gam_var_scores")
ggplot(data = species_gamvar, aes(horizon, skill, color=train_batch)) + 
  geom_point() + geom_line() + ggtitle("GAM VAR")

p2 =  ggplot(data = filter(data, species == sp & model_type == "gam_ar_scores"), aes(horizon, skill, color=train_batch)) + 
  geom_point() + geom_line() + ggtitle("GAM AR")
p3 = ggplot(data = filter(data, species == sp & model_type == "ar_scores"), aes(horizon, skill, color=train_batch)) + 
  geom_point() + geom_line() + ggtitle("AR")
