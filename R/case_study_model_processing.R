library(dplyr)
library(mvgam)
library(ggplot2)
library(patchwork)

###############################
#   Extract CRPS scores & Plot model forecasts
###############################

grab_scores = function(model, test_start){
  score <- score(forecast(model), score = "crps")
  score$test_start_newmoonnumber <- test_start
  return(score)
}

test_starts = c(363, 375, 387)
baseline_scores = vector(mode = "list", length = length(test_starts))
ar_scores <- vector(mode = "list", length = length(test_starts))
gam_ar_scores <- vector(mode = "list", length = length(test_starts))
gam_var_scores <- vector(mode = "list", length = length(test_starts))

i = 1

for (test_start in test_starts){
  print("loading gam_var")
  model_gam_var <- readRDS(paste("gam_var_inregime_output",test_start,".rds", sep=""))
  print("loading gam_ar")
  model_gam_ar <- readRDS(paste("gam_ar_inregime_output",test_start,".rds", sep=""))
  print("loading ar")
  model_ar <- readRDS(paste("ar_inregime_output",test_start,".rds", sep=""))
  print("loading baseline")
  model_baseline = readRDS(paste("baseline_inregime_output",test_start,".rds", sep=""))
  
  species = c(1:8)
  
  for(sp in species){
 
  #   # Comment out this for-loop if don't want forecast graphs   
  #   par(mfrow = c(2, 2))
  #   
     plot(model_gam_var, "forecast", series = sp)
  #   plot(model_gam_ar, "forecast", series = sp) 
  #   plot(model_ar, "forecast", series = sp) 
  #   plot(model_baseline, "forecast", series = sp) 
  # 
  
  baseline_scores[[i]] = grab_scores(model_baseline, test_start)
  ar_scores[[i]] = grab_scores(model_ar, test_start)
  gam_ar_scores[[i]] = grab_scores(model_gam_ar, test_start)
  gam_var_scores[[i]] = grab_scores(model_gam_var, test_start)
i= i+1
}

saveRDS(gam_var_scores, "gam_var_inregime_scores.rds")
saveRDS(ar_scores, "ar_inregime_scores.rds")
saveRDS(gam_ar_scores, "gam_ar_inregime_scores.rds")
saveRDS(baseline_scores, "baseline_inregime_scores.rds")
}

#############################
#   Make skill scores
#############################

gam_var_scores <- readRDS("gam_var_inregime_scores.rds")
gam_ar_scores <- readRDS("gam_ar_inregime_scores.rds")
ar_scores <- readRDS("ar_inregime_scores.rds")
baseline_scores = readRDS("baseline_inregime_scores.rds")

series = c("all_series", "DO","DM","DS","PP","PB", "PF")

get_model_skill = function(scores, series){
  model_name = deparse(substitute(scores))
  output = data.frame(model_type = character(), 
                      train_set = factor(),
                      species = character(),
                      horizon = integer(),
                      skill = numeric())
  train_sets = c(1:3)
  for (train_set in train_sets){
    for (sp in series){
      baseline_crps = baseline_scores[[train_set]][[sp]][["score"]]
      model_crps = scores[[train_set]][[sp]][["score"]]
      horizon = baseline_scores[[train_set]][[sp]][["eval_horizon"]]
      skill = (model_crps - baseline_crps)/(0-baseline_crps)
      species = rep(sp, length(horizon))
      model_type = rep(model_name, length(horizon))
      train_batch = rep(train_set,length(horizon))
      temp_df = data.frame(model_type,
                           train_batch,
                           species,
                           horizon,
                           skill)
      output = rbind(output, temp_df)
    }
  }
  return(output)
}


gam_var_skills = get_model_skill(gam_var_scores, series)
gam_ar_skills = get_model_skill(gam_ar_scores, series)
ar_skills = get_model_skill(ar_scores, series)
skills_df = rbind(gam_var_skills,gam_ar_skills)
skills_df = rbind(skills_df, ar_skills)

write.csv(skills_df, "model_skills.csv", row.names = FALSE)

