here::i_am("data_analysis/fit_ar_mods.R")

# Load packages -----------------------------------------------------------

suppressPackageStartupMessages({
  library(doParallel)
  library(here)
  library(logger)
  library(quantreg)
  library(Rcpp)
  library(tidyverse)
})

# Set up the log file -----------------------------------------------------

invisible(file.remove("fit_ar_mods.log"))
log_appender(appender_tee("fit_ar_mods.log"))

# Source code -------------------------------------------------------------

log_info("Sourcing utility functions")
source(here("data_analysis/utils.R"))
source(here("data_analysis/estimate_quantiles.R"))

# Read command-line arguments ---------------------------------------------

if (interactive()) {
  # Some of the choices below are only appropriate when analyzing hourly data
  dataset_file <- "../make_dataset/dataset.RData"
  num_cores <- min(detectCores(), 8)
  loss_type <- "ols"
  flux_thresholds <- 1e-5
  quantile_levels <- 0.9
  ords <- c(1, 5, 10, 50, 100)
  leads <- c(1, 2)
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
  dataset_file <- cmd_args[1]
  num_cores <- as.integer(cmd_args[2])
  loss_type <- cmd_args[3] # Should be "ols" or "lad"
  flux_thresholds <- cmd_args[4] %>% str_split(" ") %>% unlist() %>% as.double()
  if (identical(flux_thresholds, NA_real_)) flux_thresholds <- double()
  quantile_levels <- cmd_args[5] %>% str_split(" ") %>% unlist() %>% as.double()
  if (identical(quantile_levels, NA_real_)) quantile_levels <- double()
  ords <- cmd_args[6] %>% str_split(" ") %>% unlist() %>% as.integer()
  leads <- cmd_args[7] %>% str_split(" ") %>% unlist() %>% as.integer()
}

# Fit models --------------------------------------------------------------

log_info("Loading dataset")
load(dataset_file)
registerDoParallel(cores = num_cores)
log_info("Started fitting models")
ar_results <- foreach (i = 1:nrow(training_windows)) %dopar% {
  time_elapsed <- system.time(
    result <- with(
      training_windows,
      fit_ar_mod(
        dataset,
        first_train_time[i], last_train_time[i],
        loss_type,
        flux_thresholds, quantile_levels, ords, leads
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
ar_summary <- ar_results %>%
  select(ord:quantile_level, obs, pred) %>%
  mutate(
    flux_threshold = if_else(type == "threshold", flux_threshold, NA_real_),
    quantile_level = if_else(type == "level", quantile_level, NA_real_)
  ) %>%
  nest(.by = ord:quantile_level, .key = "obs_preds") %>%
  mutate(
    confusion_mat = map(obs_preds, make_confusion_mat),
    event_rate = map_dbl(confusion_mat, calc_event_rate),
    alarm_rate = map_dbl(confusion_mat, calc_alarm_rate),
    precision = map_dbl(confusion_mat, calc_precision),
    tpr = map_dbl(confusion_mat, calc_tpr),
    fpr = map_dbl(confusion_mat, calc_fpr),
    tss = map_dbl(confusion_mat, calc_tss)
  )

assign(str_glue("{loss_type}_ar_results"), ar_results)
assign(str_glue("{loss_type}_ar_summary"), ar_summary)
save(
  list = str_c(loss_type, c("ar_results", "ar_summary"), sep = "_"),
  file = str_glue("{loss_type}_ar_results.RData")
)
log_info("Done")