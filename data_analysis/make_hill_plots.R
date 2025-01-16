here::i_am("data_analysis/make_hill_plots.R")

suppressPackageStartupMessages({
  library(doParallel)
  library(here)
  library(logger)
  library(Rcpp)
  library(tidyverse)
})

invisible(file.remove("make_hill_plots.log"))
log_appender(appender_tee("make_hill_plots.log"))

log_info("Sourcing utility functions")
source(here("data_analysis/utils.R"))

if (interactive()) {
  dataset_file <- "../make_dataset/dataset.RData"
  num_cores <- min(detectCores(), 8)
  max_num_top_ord_stats <- 100
  alpha_hat <- 1.25
} else {
  cmd_args <- commandArgs(trailingOnly = TRUE)
  dataset_file <- cmd_args[1]
  num_cores <- as.integer(cmd_args[2])
  max_num_top_ord_stats <- as.integer(cmd_args[3])
  alpha_hat <- as.double(cmd_args[4])
}

log_info("Loading dataset")
load(dataset_file)
registerDoParallel(cores = num_cores)
log_info("Started making plots")
plots <- foreach (i = 1:nrow(training_windows)) %dopar% {
  time_elapsed <- system.time(
    plot <- with(
      training_windows,
      make_hill_plot(
        dataset,
        first_train_time[i], last_train_time[i],
        label[i],
        max_num_top_ord_stats,
        alpha_hat
      )
    )
  )["elapsed"]
  time_elapsed <- signif(time_elapsed, 2)
  log_info(
    str_glue("Made plot for {training_windows$label[i]} in {time_elapsed}s")
  )
  plot
}

log_info("Saving plots")
pdf("hill_plots.pdf")
walk(plots, print)
invisible(dev.off())
log_info("Done")