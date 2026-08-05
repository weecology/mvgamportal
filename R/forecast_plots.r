library(ggplot2)
library(cowplot)
library(dplyr)

forecast_plot <- function(model, model_name, model_scores, species_list, test_start) {
p1<-tryCatch({plot(forecast(model), series=1)}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)}) + ggtitle(paste(model_name,test_start)) + theme(plot.title = element_text(hjust = 0.5))
p2<-tryCatch({if(length(species_list)>1) {plot(forecast(model), series=2)} else {NULL}}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p3<-tryCatch({if(length(species_list)>2) {plot(forecast(model), series=3)} else {NULL}}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p4<-tryCatch({if(length(species_list)>3) {plot(forecast(model), series=4)} else {NULL}}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p5<-tryCatch({if(length(species_list)>4) {plot(forecast(model), series=5)} else {NULL}}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p6<-tryCatch({if(length(species_list)>5) {plot(forecast(model), series=6)} else {NULL}}, error = function(e) {message("An error occurred during plotting or saving: ", e$message)})
p7<-ggplot(model_scores, aes(x=eval_horizon,y=score,color=species)) + geom_point() + theme_minimal()
p8<-mcmc_plot(model, type = 'rhat_hist')
topleft<-plot_grid(p1,p2,p3,p4)
top2rows<-plot_grid(topleft,p7,rel_widths = c(.67,.33))
bottomrow<-plot_grid(p5,p6,p8,nrow=1)
plot_grid(top2rows,bottomrow,ncol=1,nrow=2,rel_heights = c(.67,.33))
}

trace_plot <- function(model) {
  tryCatch(
    {
      mcmc_plot(model, type = 'trace')
    },
    error = function(e) {
      message("An error occurred during plotting or saving: ", e$message)
    }
  )
}

# Write the forecast panel for one fitted model, plus an MCMC trace plot when
# trace_plots is TRUE. run_label identifies the model and the species group it
# was fit to, and names the files.
plot_run_forecasts <- function(
  model,
  run_label,
  model_scores,
  species_list,
  test_start,
  trace_plots = TRUE
) {
  png(
    paste0("figures/", run_label, "_", test_start, ".png"),
    width = 1500,
    height = 1000
  )
  # print() is required: plots built inside a function are not auto-printed
  print(forecast_plot(
    model,
    run_label,
    model_scores,
    species_list,
    test_start
  ))
  dev.off()

  if (trace_plots) {
    png(
      paste0("figures/", run_label, "_trace_", test_start, ".png"),
      width = 1500,
      height = 1500
    )
    print(trace_plot(model))
    dev.off()
  }
}
