library(purrr)
library(dplyr)
library(loo)
library(ggplot2)

loos_invndvi <- readRDS("loos.rds")
loos_ndvi <- readRDS("loos_ndvi.rds")

newmoons <- sapply(loos_invndvi, function(x) x$test_start_newmoonnumber)
models <- rep(c("Baseline","AR","GAM_AR","GAM_VAR","Simple"),6)

pairwise_comparisons <- map2(loos_ndvi, loos_invndvi, function(mod_a, mod_b) {
  loo_compare(list(NDVI = mod_a, Inv_NDVI = mod_b))
})

loo_ts <- imap_dfr(pairwise_comparisons, function(comp_matrix, element_name) {
  
  comp_df <- as.data.frame(comp_matrix)
  r_names <- rownames(comp_df)
  
  raw_elpd_ndvi    <- comp_df[r_names == "NDVI", "elpd_loo"][1]
  raw_elpd_invndvi <- comp_df[r_names == "Inv_NDVI", "elpd_loo"][1]
  
  se_difference    <- comp_df[2, "se_diff"][1]

  elpd_diff_fixed  <- raw_elpd_ndvi - raw_elpd_invndvi
  
  data.frame(
    elpd_diff  = as.numeric(elpd_diff_fixed),
    se_diff    = as.numeric(se_difference),
    lower      = as.numeric(elpd_diff_fixed - (1.96 * se_difference)),
    upper      = as.numeric(elpd_diff_fixed + (1.96 * se_difference)),
    stringsAsFactors = FALSE
  )
})

loo_ts$newmoon = newmoons
loo_ts$model = models

ggplot(loo_ts, aes(x = newmoon, y = elpd_diff)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey10", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_line(aes(color = model), linewidth = 1) +
  geom_point(aes(color = model), size = 2) +
  facet_wrap(~model, scales = "free_y") +
  labs(
    subtitle = "Difference in ELPD (NDVI - Inv_NDVI).",
    x = "Test Start Newmoonnumber",
    y = "Δ ELPD (Higher = NDVI Wins)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  ) + scale_y_continuous(trans = "pseudo_log")

loos_invndvi <- Map(c, loos_invndvi, model = models)
loos_ndvi <- Map(c, loos_ndvi, model = models)

k_invndvi <- do.call(rbind, lapply(seq_along(loos_invndvi), function(i) {
  data.frame(
    Observation = seq_along(pareto_k_values(loos_invndvi[[i]])),
    Pareto_k = pareto_k_values(loos_invndvi[[i]]),
    Model = loos_invndvi[[i]]$model,
    newmoonnumber = loos_invndvi[[i]]$test_start_newmoonnumber
  )
}))

k_ndvi <- do.call(rbind, lapply(seq_along(loos_ndvi), function(i) {
  data.frame(
    Observation = seq_along(pareto_k_values(loos_ndvi[[i]])),
    Pareto_k = pareto_k_values(loos_ndvi[[i]]),
    Model = loos_ndvi[[i]]$model,
    newmoonnumber = loos_ndvi[[i]]$test_start_newmoonnumber
  )
}))

ggplot(k_invndvi, aes(x = newmoonnumber, y = Pareto_k, color = Pareto_k > 0.5)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  facet_wrap(~ Model) +
  theme_minimal() +
  labs(x = "Observation Index", y = "Pareto k estimate", color = "Warning (k > 0.5)")

ggplot(k_ndvi, aes(x = newmoonnumber, y = Pareto_k, color = Pareto_k > 0.5)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  facet_wrap(~ Model) +
  theme_minimal() +
  labs(x = "Observation Index", y = "Pareto k estimate", color = "Warning (k > 0.5)")
