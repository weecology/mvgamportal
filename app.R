# Interactive exploration of the skill scores in scores.rds.
#
# Launch with shiny::runApp(), or the Run App button in RStudio. Sourcing this
# file instead defines the functions below without starting the app, so the
# plots can be reused for a one-off figure:
#
#   source("app.R")
#   ggsave("skill.png", plot_skill_timeseries(scores, -5, TRUE))
#
# Two plots — a time series of skill scores and a histogram that pools the
# negative tail into a single bin — over any combination of species, models,
# forecast horizons and species sets present in the scores. Clicking a point in
# the time series opens the forecast panel for the run behind it.

library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(tidyr)

# Ecological regime periods shaded in the time series, in newmoonnumbers
regime_periods <- data.frame(xmin = c(140, 263, 396), xmax = c(230, 278, 411))

# The data as the time series draws it. Shared with the click handler, which
# has to hit test against the same coordinates the points were plotted at.
prepare_timeseries <- function(scores, y_min) {
  scores |>
    # Newmoons a series has no score for are filled in as NA so that windows
    # that failed to fit are drawn as breaks in the line rather than
    # interpolated straight over
    group_by(species, model, species_set, eval_horizon) |>
    complete(newmoonnumber = full_seq(newmoonnumber, 1)) |>
    ungroup() |>
    # Scores past the bottom of the axis are pinned to it and drawn larger, so
    # an excursion off the panel is still visible as one. The unpinned value is
    # kept so a clicked point can be reported at its real score.
    mutate(below = !is.na(skill_score) & skill_score < y_min,
           skill_score_actual = skill_score,
           skill_score = pmax(skill_score, y_min))
}

plot_skill_timeseries <- function(scores, y_min, show_regimes) {
  d <- prepare_timeseries(scores, y_min)

  # Models share a set of axes per species, separated by line type and point
  # shape, so they can be compared directly against each other
  p <- ggplot(d, aes(x = newmoonnumber, y = skill_score,
                     color = factor(eval_horizon),
                     linetype = model, shape = model))
  if (show_regimes) {
    # Clipped to the selected windows so the bands don't stretch the x axis
    # out over newmoons with no scores in them
    window <- range(d$newmoonnumber)
    regimes <- regime_periods |>
      filter(xmax >= window[1], xmin <= window[2]) |>
      mutate(xmin = pmax(xmin, window[1]), xmax = pmin(xmax, window[2]))
    p <- p + geom_rect(data = regimes,
                       aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                       inherit.aes = FALSE, fill = "lightgrey", alpha = 0.4)
  }
  p <- p +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
    geom_line(alpha = 0.6, na.rm = TRUE) +
    geom_point(aes(size = below), na.rm = TRUE) +
    facet_wrap(~species, ncol = 1) +
    # clip = "off" keeps points sitting exactly on a limit from being sliced in
    # half by the panel edge, which expand = FALSE would otherwise do
    coord_cartesian(ylim = c(y_min, 1), expand = FALSE, clip = "off") +
    scale_color_viridis_d(end = 0.9) +
    scale_size_manual(values = c(`FALSE` = 0.9, `TRUE` = 2.4), guide = "none") +
    labs(x = "New moon number", y = "Skill score", color = "Horizon",
         linetype = "Model", shape = "Model") +
    theme_minimal()

  if (any(d$below)) {
    p <- p + labs(caption = paste("Enlarged points fall below", y_min))
  }
  p
}

# Everything below `pool_below` is drawn as one bin at the left edge, since the
# negative tail runs to -Inf and would otherwise flatten the rest of the plot.
plot_skill_histogram <- function(scores, binwidth, pool_below) {
  clamp_value <- pool_below - binwidth / 2
  d <- scores |>
    mutate(skill_score_clamped = pmax(skill_score, clamp_value),
           pooled = skill_score < pool_below)

  ggplot(d, aes(x = skill_score_clamped, fill = pooled)) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = "lightgrey", alpha = 0.6) +
    geom_histogram(binwidth = binwidth, boundary = 0, color = "white") +
    scale_fill_manual(values = c(`FALSE` = "steelblue", `TRUE` = "steelblue4"),
                      labels = c(`FALSE` = "binned",
                                 `TRUE` = paste0("pooled (< ", pool_below, ")"))) +
    facet_grid(species ~ model, scales = "free_y") +
    coord_cartesian(xlim = c(pool_below - binwidth, 1), expand = FALSE) +
    labs(x = paste0("Skill score (values < ", pool_below, " pooled)"),
         y = "Count", fill = NULL) +
    theme_minimal()
}

skill_score_app <- function(scores_path = "scores.rds", figures_path = "figures") {
  # BASELINE is the reference every skill score is computed against, so its
  # own skill score is always 0
  scores <- readRDS(scores_path) |>
    filter(model != "BASELINE", is.finite(skill_score))

  # Lets the forecast panels be served straight out of the figures directory,
  # which is otherwise outside anything the app can serve
  addResourcePath("forecast_figures", normalizePath(figures_path))

  species <- sort(unique(scores$species))
  models <- sort(unique(scores$model))
  horizons <- sort(unique(scores$eval_horizon))
  species_sets <- sort(unique(scores$species_set))
  newmoons <- range(scores$newmoonnumber)

  ui <- page_sidebar(
    title = "Portal forecast skill scores",
    sidebar = sidebar(
      width = 260,
      selectizeInput("species", "Species", species, selected = species,
                     multiple = TRUE),
      selectizeInput("models", "Models", models, selected = models,
                     multiple = TRUE),
      selectizeInput("horizons", "Forecast horizons", horizons,
                     selected = horizons, multiple = TRUE),
      selectizeInput("species_sets", "Species sets", species_sets,
                     selected = species_sets, multiple = TRUE),
      # sep = "" because these are newmoon indices, not counts to be read with
      # a thousands separator
      sliderInput("newmoons", "New moon range", min = newmoons[1],
                  max = newmoons[2], value = newmoons, step = 1, sep = ""),
      hr(),
      numericInput("y_min", "Time series y-axis minimum", -5, step = 1),
      checkboxInput("show_regimes", "Shade regime periods", TRUE),
      numericInput("binwidth", "Histogram bin width", 0.1, min = 0.01, step = 0.05),
      numericInput("pool_below", "Pool skill scores below", -1, step = 0.5)
    ),
    navset_card_tab(
      nav_panel(
        "Time series",
        plotOutput("timeseries", height = "auto", click = "timeseries_click"),
        helpText("Click a point to see the forecast panel for that run.")
      ),
      nav_panel("Histogram", plotOutput("histogram", height = "auto"))
    )
  )

  server <- function(input, output) {
    selected <- reactive({
      req(input$newmoons)
      scores |>
        filter(species %in% input$species,
               model %in% input$models,
               eval_horizon %in% as.numeric(input$horizons),
               species_set %in% input$species_sets,
               between(newmoonnumber, input$newmoons[1], input$newmoons[2]))
    })

    plot_height <- function() max(400, 200 * length(input$species))

    output$timeseries <- renderPlot({
      req(input$y_min)
      validate(need(nrow(selected()) > 0, "No scores match the current selection."))
      plot_skill_timeseries(selected(), input$y_min, input$show_regimes)
    }, height = plot_height)

    output$histogram <- renderPlot({
      req(input$binwidth, input$pool_below)
      validate(need(nrow(selected()) > 0, "No scores match the current selection."))
      plot_skill_histogram(selected(), input$binwidth, input$pool_below)
    }, height = plot_height)

    # Points the line was interpolated over have no run behind them, so they
    # are dropped rather than hit tested against
    timeseries_points <- reactive({
      prepare_timeseries(selected(), input$y_min) |>
        filter(!is.na(skill_score))
    })

    # Every score comes from a run that plot_run_forecasts() drew a panel for,
    # named for the run and the newmoon its test window starts at, so a clicked
    # point identifies a file in the figures directory
    observeEvent(input$timeseries_click, {
      point <- nearPoints(timeseries_points(), input$timeseries_click,
                          maxpoints = 1, threshold = 10)
      req(nrow(point) == 1)

      run_label <- paste(point$model, point$species_set, sep = "_")
      figure <- paste0(run_label, "_", point$test_start_newmoonnumber, ".png")
      showModal(modalDialog(
        title = paste0(run_label, " forecast from newmoon ",
                       point$test_start_newmoonnumber, " — ", point$species,
                       ", horizon ", point$eval_horizon, ", skill score ",
                       round(point$skill_score_actual, 2)),
        size = "xl",
        easyClose = TRUE,
        footer = modalButton("Close"),
        if (file.exists(file.path(figures_path, figure))) {
          # The panel is drawn at 1500x1000 and has to shrink a long way to fit
          # the dialog, so it links out to itself at full size
          tags$a(href = file.path("forecast_figures", figure), target = "_blank",
                 tags$img(src = file.path("forecast_figures", figure),
                          width = "100%"))
        } else {
          paste0("No forecast panel at ", file.path(figures_path, figure), ".")
        }
      ))
    })
  }

  shinyApp(ui, server)
}

skill_score_app()
