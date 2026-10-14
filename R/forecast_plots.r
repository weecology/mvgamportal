library(ggplot2)
library(dplyr)

forecast_plot <- function(model, model_name, model_scores, species, sp_series, test_start) {
  tryCatch(
    {
      plot(forecast(model), series=sp_series)
    }, 
    error = function(e) {
      message("An error occurred during plotting or saving: ", e$message)
    }) + 
    ggtitle(paste(model_name,test_start,species)) + 
    labs(subtitle = paste0("Score:", round(mean(model_scores$score,na.rm=TRUE),4),
      "  Mean rhat:", round(mean(model_scores$rhat,na.rm=TRUE),4),
      "  High Rhats:", sum(model_scores$prhat_high,na.rm=TRUE),
      "  Divergences:", sum(model_scores$n_divergences,na.rm=TRUE)))
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
  species_list[species_list == "all"] <- "all_series"
  for (sp_series in seq_along(species_list)) {
    species <- species_list[[sp_series]]
    png(
      paste0("figures/", run_label, "_", test_start, "_", species, ".png"),
      width = 1000,
      height = 500
    )
    # print() is required: plots built inside a function are not auto-printed
    print(forecast_plot(
      model,
      run_label,
      model_scores,
      species, 
      sp_series,
      test_start
    ))
  dev.off()
  
    if (trace_plots) {
      png(
        paste0("figures/", run_label, "_trace_", test_start, "_", species, ".png"),
        width = 1500,
        height = 1500
      )
      print(trace_plot(model))
      dev.off()
    }
  }
}