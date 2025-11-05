library(ggplot2)
library(cowplot)

png(paste0("figures/ar_",target,".png"), width = 1500, height = 1000) 

p1<-plot(forecast(ar_model, level=c(.1,.25,.5)), series=1)
p2<-plot(forecast(ar_model, level=c(.1,.25,.5)), series=2) + ggtitle(paste("AR",target)) + theme(plot.title = element_text(hjust = 0.5))
p3<-plot(forecast(ar_model, level=c(.1,.25,.5)), series=3)
p4<-plot(forecast(ar_model, level=c(.1,.25,.5)), series=4)
p5<-plot(forecast(ar_model, level=c(.1,.25,.5)), series=5)
p6<-plot(forecast(ar_model, level=c(.1,.25,.5)), series=6)
p7<-ggplot(ar, aes(x=eval_horizon,y=skill_score,color=species)) + geom_point() + theme_minimal()
p8<-mcmc_plot(ar_model, type = 'rhat_hist')
topleft<-plot_grid(p1,p2,p3,p4)
top2rows<-plot_grid(topleft,p7,rel_widths = c(.67,.33))
bottomrow<-plot_grid(p5,p6,p8,nrow=1)
plot_grid(top2rows,bottomrow,ncol=1,nrow=2,rel_heights = c(.67,.33))

dev.off()

png(paste0("figures/gamar_",target,".png"), width = 1500, height = 1000) 

p1<-plot(forecast(gam_ar_model, level=c(.1,.25,.5)), series=1)
p2<-plot(forecast(gam_ar_model, level=c(.1,.25,.5)), series=2) + ggtitle(paste("GAM AR",target)) + theme(plot.title = element_text(hjust = 0.5))
p3<-plot(forecast(gam_ar_model, level=c(.1,.25,.5)), series=3)
p4<-plot(forecast(gam_ar_model, level=c(.1,.25,.5)), series=4)
p5<-plot(forecast(gam_ar_model, level=c(.1,.25,.5)), series=5)
p6<-plot(forecast(gam_ar_model, level=c(.1,.25,.5)), series=6)
p7<-ggplot(gam_ar, aes(x=eval_horizon,y=skill_score,color=species)) + geom_point() + theme_minimal()
p8<-mcmc_plot(gam_ar_model, type = 'rhat_hist')
topleft<-plot_grid(p1,p2,p3,p4)
top2rows<-plot_grid(topleft,p7,rel_widths = c(.67,.33))
bottomrow<-plot_grid(p5,p6,p8,nrow=1)
plot_grid(top2rows,bottomrow,ncol=1,nrow=2,rel_heights = c(.67,.33))

dev.off()


png(paste0("figures/gamvar_",target,".png"), width = 1500, height = 1000) 

p1<-tryCatch({plot(forecast(gam_var_model, level=c(.1,.25,.5)), series=1)}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p2<-tryCatch({plot(forecast(gam_var_model, level=c(.1,.25,.5)), series=2) + ggtitle(paste("GAM VAR",target)) + theme(plot.title = element_text(hjust = 0.5))}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p3<-tryCatch({plot(forecast(gam_var_model), series=3)}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p4<-tryCatch({plot(forecast(gam_var_model, level=c(.1,.25,.5)), series=4)}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p5<-tryCatch({plot(forecast(gam_var_model, level=c(.1,.25,.5)), series=5)}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p6<-tryCatch({plot(forecast(gam_var_model, level=c(.1,.25,.5)), series=6)}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p7<-ggplot(gam_var, aes(x=eval_horizon,y=skill_score,color=species)) + geom_point() + theme_minimal()
p8<-mcmc_plot(gam_var_model, type = 'rhat_hist')
topleft<-plot_grid(p1,p2,p3,p4)
top2rows<-plot_grid(topleft,p7,rel_widths = c(.67,.33))
bottomrow<-plot_grid(p5,p6,p8,nrow=1)
plot_grid(top2rows,bottomrow,ncol=1,nrow=2,rel_heights = c(.67,.33))

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
