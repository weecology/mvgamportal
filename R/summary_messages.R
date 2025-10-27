# Check samples

for(i in 1:length(gam_var_summaries)) {
if(!("✔ No issues with effective samples per iteration" %in%  gam_var_summaries[[i]]$diagnostics$sampler_message)){
  print(paste(gam_var_summaries[[i]]$test_start_newmoonnumber , gam_var_summaries[[i]]$diagnostics$sampler_message[1:2]))
  } 
if(!("✔ Rhat looks good for all parameters" %in%  gam_var_summaries[[i]]$diagnostics$sampler_message)){
  print(paste(gam_var_summaries[[i]]$test_start_newmoonnumber , gam_var_summaries[[i]]$diagnostics$sampler_message[2:3]))
  } 
if(!("✔ No issues with divergences" %in%  gam_var_summaries[[i]]$diagnostics$sampler_message)){
  print(paste(gam_var_summaries[[i]]$test_start_newmoonnumber , gam_var_summaries[[i]]$diagnostics$sampler_message[3:4]))
  } 
if(!("✔ No issues with maximum tree depth"  %in%  gam_var_summaries[[i]]$diagnostics$sampler_message)){
  print(paste(gam_var_summaries[[i]]$test_start_newmoonnumber , gam_var_summaries[[i]]$diagnostics$sampler_message[4:5]))
  }   
  
  }

for(i in 1:length(gam_ar_summaries)) {
  if(!("✔ No issues with effective samples per iteration" %in%  gam_ar_summaries[[i]]$diagnostics$sampler_message)){
    print(paste(gam_ar_summaries[[i]]$test_start_newmoonnumber , gam_ar_summaries[[i]]$diagnostics$sampler_message[1:2]))
  } 
  if(!("✔ Rhat looks good for all parameters" %in%  gam_ar_summaries[[i]]$diagnostics$sampler_message)){
    print(paste(gam_ar_summaries[[i]]$test_start_newmoonnumber , gam_ar_summaries[[i]]$diagnostics$sampler_message[2:3]))
  } 
  if(!("✔ No issues with divergences" %in%  gam_ar_summaries[[i]]$diagnostics$sampler_message)){
    print(paste(gam_ar_summaries[[i]]$test_start_newmoonnumber , gam_ar_summaries[[i]]$diagnostics$sampler_message[3:4]))
  } 
  if(!("✔ No issues with maximum tree depth"  %in%  gam_ar_summaries[[i]]$diagnostics$sampler_message)){
    print(paste(gam_ar_summaries[[i]]$test_start_newmoonnumber , gam_ar_summaries[[i]]$diagnostics$sampler_message[4:5]))
  }   
  
}

for(i in 1:length(ar_summaries)) {
  if(!("✔ No issues with effective samples per iteration" %in%  ar_summaries[[i]]$diagnostics$sampler_message)){
    print(paste(ar_summaries[[i]]$test_start_newmoonnumber , ar_summaries[[i]]$diagnostics$sampler_message[1:2]))
  } 
  if(!("✔ Rhat looks good for all parameters" %in%  ar_summaries[[i]]$diagnostics$sampler_message)){
    print(paste(ar_summaries[[i]]$test_start_newmoonnumber , ar_summaries[[i]]$diagnostics$sampler_message[2:3]))
  } 
  if(!("✔ No issues with divergences" %in%  ar_summaries[[i]]$diagnostics$sampler_message)){
    print(paste(ar_summaries[[i]]$test_start_newmoonnumber , ar_summaries[[i]]$diagnostics$sampler_message[3:4]))
  } 
  if(!("✔ No issues with maximum tree depth"  %in%  ar_summaries[[i]]$diagnostics$sampler_message)){
    print(paste(ar_summaries[[i]]$test_start_newmoonnumber , ar_summaries[[i]]$diagnostics$sampler_message[4:5]))
  }   
  
}
