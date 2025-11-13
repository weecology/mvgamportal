library(dplyr)
library(mvgam)

###############################
#   Functions to extract CRPS scores & Plot model forecasts
###############################

grab_scores = function(model, test_start){
  score <- score(forecast(model), score = "crps")
  score$test_start_newmoonnumber <- test_start
  return(score)
}

# note: need to change hard-coded series numbers to a number read from the file
process_modeloutput = function(within_starts, outside_starts) {
  test_starts = c(within_starts, outside_starts)
  
  baseline_scores = vector(mode = "list", length = length(test_starts))
  ar_scores <- vector(mode = "list", length = length(test_starts))
  gam_ar_scores <- vector(mode = "list", length = length(test_starts))
  gam_var_scores <- vector(mode = "list", length = length(test_starts))
  
  path = "./case_study/model_outputs/"
  i = 1
  for (test_start in test_starts) {
    if (test_start %in% outside_starts) {
      type = "out"
    } else {
      type = "in"
    }
    print(paste("loading gam_var", test_start, sep = ""))
    model_gam_var <- readRDS(paste(path,"gam_var_",type,"regime_output",test_start,".rds",
                                   sep = ""))
    
    print(paste("loading gam_ar", test_start, sep = ""))
    model_gam_ar <- readRDS(paste(path,"gam_ar_",type,"regime_output",test_start,".rds",
                                  sep = ""))
    
    print(paste("loading ar", test_start, sep = ""))
    model_ar <- readRDS(paste(path, "ar_", type, "regime_output", test_start, ".rds", 
                              sep = ""))
    
    print(paste("loading baseline", test_start, sep = ""))
    model_baseline = readRDS(paste(path,"baseline_",type,"regime_output", test_start,".rds",
                                   sep = "" ))
    
    species = c(1:6)
    
      baseline_scores[[i]] = grab_scores(model_baseline, test_start)
      print(paste("got baseline", test_start, sep=" "))
      ar_scores[[i]] = grab_scores(model_ar, test_start)
      print(paste("got ar", test_start, sep=" "))
      gam_ar_scores[[i]] = grab_scores(model_gam_ar, test_start)
      print(paste("got gam_ar", test_start,  sep=" "))
      gam_var_scores[[i]] = grab_scores(model_gam_var, test_start)
      print(paste("got gam var", test_start, sep=" "))
    i = i + 1
  }
  
  saveRDS(gam_var_scores, paste(path, "gam_var_scores.rds", sep = ""))
  saveRDS(ar_scores, paste(path, "ar_scores.rds", sep = ""))
  saveRDS(gam_ar_scores, paste(path, "gam_ar_scores.rds", sep = ""))
  saveRDS(baseline_scores, paste(path, "baseline_scores.rds", sep = ""))
  
}

#############################
#   Functions to calculate skill scores
#############################

calc_model_skill = function(scores, series){
  model_name = deparse(substitute(scores))
  output = data.frame(model_type = character(), 
                      train_set = factor(),
                      species = character(),
                      horizon = integer(),
                      skill = numeric())
  path = "./case_study/model_outputs/"
  baseline_scores = readRDS(paste(path,"baseline_scores.rds", sep=""))
  train_sets = c(1:4)
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

get_model_skill = function(series) {
  path = "./case_study/model_outputs/"
  gam_var_scores <- readRDS(paste(path, "gam_var_scores.rds", sep = ""))
  gam_ar_scores <- readRDS(paste(path, "gam_ar_scores.rds", sep = ""))
  ar_scores <- readRDS(paste(path, "ar_scores.rds", sep = ""))
  
  gam_var_skills = calc_model_skill(gam_var_scores, series)
  gam_ar_skills = calc_model_skill(gam_ar_scores, series)
  ar_skills = calc_model_skill(ar_scores, series)
  skills_df = rbind(gam_var_skills, gam_ar_skills)
  skills_df = rbind(skills_df, ar_skills)
  write.csv(skills_df,
            paste(path, "model_skills.csv", sep = ""),
            row.names = FALSE)
}

######## Executable code
in_starts = c(363, 375, 387)
out_starts = c(412)


# generates score files for within and out of regime forecasts

process_modeloutput(in_starts, out_starts)

# generates skills files for within and out of regime forecasts
# FYI: would be nice to have the code figure out what species are in the file so it
# can adjust dynamically to any changes.
series = c("all_series", "DO","DM", "PP","PB", "PF")
get_model_skill(series)



