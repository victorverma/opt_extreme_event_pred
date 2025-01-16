here::i_am("data_analysis/fit_farima_mods.R")

# Load packages -----------------------------------------------------------

suppressPackageStartupMessages({
  library(doParallel)
  library(extRemes)
  library(here)
  library(logger)
  library(Rcpp)
  library(tidyverse)
})

# Set up the log file -----------------------------------------------------

invisible(file.remove("fit_farima_mods.log"))
log_appender(appender_tee("fit_farima_mods.log"))

# Source code -------------------------------------------------------------

log_info("Sourcing utility functions")
source(here("data_analysis/utils.R"))
source(here("data_analysis/estimate_quantiles.R"))

# Read command-line arguments ---------------------------------------------

if (interactive()) {
  # Some of the choices below are only appropriate when analyzing hourly data
  dataset_file <- "../make_dataset/dataset.RData"
  num_cores <- min(detectCores(), 8)
  alpha_hat <- 1.25
  d_hat <- NULL
  num_coefs <- 1e5
  flux_thresholds <- 1e-5
  quantile_levels <- 0.9
  ords <- c(1, 5, 10, 50, 100)
  leads <- c(1, 2)
  plot_obj_fun <- TRUE
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
  dataset_file <- cmd_args[1]
  num_cores <- as.integer(cmd_args[2])
  alpha_hat <- if (cmd_args[3] != "NULL") as.double(cmd_args[3])
  d_hat <- if (cmd_args[4] != "NULL") as.double(cmd_args[4])
  num_coefs <- as.integer(cmd_args[5])
  flux_thresholds <- cmd_args[6] %>% str_split(" ") %>% unlist() %>% as.double()
  if (identical(flux_thresholds, NA_real_)) flux_thresholds <- double()
  quantile_levels <- cmd_args[7] %>% str_split(" ") %>% unlist() %>% as.double()
  if (identical(quantile_levels, NA_real_)) quantile_levels <- double()
  ords <- cmd_args[8] %>% str_split(" ") %>% unlist() %>% as.integer()
  leads <- cmd_args[9] %>% str_split(" ") %>% unlist() %>% as.integer()
  plot_obj_fun <- as.logical(cmd_args[10])
}

# Fit models --------------------------------------------------------------

log_info("Loading dataset")
load(dataset_file)
registerDoParallel(cores = num_cores)
log_info("Started fitting models")
results <- foreach (i = 1:nrow(training_windows)) %dopar% {
  time_elapsed <- system.time(
    result <- with(
      training_windows,
      fit_farima_mod(
        dataset,
        first_train_time[i], last_train_time[i],
        alpha_hat, d_hat,
        num_coefs,
        flux_thresholds, quantile_levels, ords, leads,
        plot_obj_fun, label[i]
      )
    )
  )["elapsed"]
  time_elapsed <- signif(time_elapsed, 2)
  log_info(
    str_glue("Fit model for {training_windows$label[i]} in {time_elapsed}s")
  )
  result
} %>%
  add_column(training_windows, result = .) %>%
  unnest(result)

# Summarize and save the results ------------------------------------------

log_info("Summarizing and saving model-fitting results")
farima_results <- results %>% select(!obj_fun_plot) %>% unnest(preds)
farima_summary <- farima_results %>%
  select(ord:quantile_level, obs, pred) %>%
  mutate(
    flux_threshold = if_else(type == "threshold", flux_threshold, NA_real_),
    quantile_level = if_else(type == "level", quantile_level, NA_real_)
  ) %>%
  nest(.by = ord:quantile_level, .key = "obs_preds") %>%
  mutate(
    confusion_mat = map(obs_preds, make_confusion_mat),
    precision = map_dbl(confusion_mat, calc_precision),
    event_rate = map_dbl(confusion_mat, calc_event_rate),
    alarm_rate = map_dbl(confusion_mat, calc_alarm_rate),
    tpr = map_dbl(confusion_mat, calc_tpr),
    fpr = map_dbl(confusion_mat, calc_fpr),
    tss = map_dbl(confusion_mat, calc_tss)
  )
save(farima_results, farima_summary, file = "farima_results.RData")

if (plot_obj_fun) {
  log_info("Saving objective function plots")
  pdf("obj_fun_plots.pdf")
  walk(results$obj_fun_plot, print)
  invisible(dev.off())
}
log_info("Done")