library(ggplot2)
library(cowplot)

forecast_plot <- function(model,model_name,skill_scores,species_list,target) {
p1<-tryCatch({plot(forecast(model), series=1)}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p2<-tryCatch({plot(forecast(model), series=2)}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)}) + ggtitle(paste(model_name,target)) + theme(plot.title = element_text(hjust = 0.5))
p3<-tryCatch({plot(forecast(model), series=3)}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p4<-tryCatch({if(length(species_list)>3) {plot(forecast(model), series=4)} else {NULL}}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p5<-tryCatch({if(length(species_list)>4) {plot(forecast(model), series=5)} else {NULL}}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p6<-tryCatch({if(length(species_list)>5) {plot(forecast(model), series=6)} else {NULL}}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p7<-ggplot(skill_scores, aes(x=eval_horizon,y=score,color=species)) + geom_point() + theme_minimal()
p8<-mcmc_plot(model, type = 'rhat_hist')
topleft<-plot_grid(p1,p2,p3,p4)
top2rows<-plot_grid(topleft,p7,rel_widths = c(.67,.33))
bottomrow<-plot_grid(p5,p6,p8,nrow=1)
plot_grid(top2rows,bottomrow,ncol=1,nrow=2,rel_heights = c(.67,.33))
}

png(paste0("figures/baseline_",target,".png"), width = 1500, height = 1000) 
forecast_plot(baseline_model,"Baseline",baseline,data_split$species_list,target)
dev.off()

png(paste0("figures/ar_",target,".png"), width = 1500, height = 1000) 
forecast_plot(ar_model,"AR",ar,data_split$species_list,target)
dev.off()

png(paste0("figures/gamar_",target,".png"), width = 1500, height = 1000) 
forecast_plot(gam_ar_model,"GAM AR",gam_ar,data_split$species_list,target)
dev.off()

png(paste0("figures/gamvar_",target,".png"), width = 1500, height = 1000) 
forecast_plot(gam_var_model,"GAM VAR",gam_var,data_split$species_list,target)
dev.off()

png(paste0("figures/ar_trace_",target,".png"), width = 1500, height = 1500) 
mcmc_plot(ar_model, type = 'trace')  
dev.off()

png(paste0("figures/gamar_trace_",target,".png"), width = 1500, height = 1500) 
tryCatch({mcmc_plot(gam_ar_model, type = 'trace')}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
dev.off()

png(paste0("figures/gamvar_trace_",target,".png"), width = 2000, height = 2000) 
tryCatch({mcmc_plot(gam_var_model, type = 'trace', variable = 'trend_params')}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)}) 
dev.off()




# mcmc_plot(
#   object = gam_var_model,
#   variable = "trend",
#   regex=TRUE,
#   type = "areas"
# )
#
# mcmc_plot(gam_ar_model, type = 'intervals')
