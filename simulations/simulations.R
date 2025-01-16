here::i_am("simulations/simulations.R")

# Load packages -----------------------------------------------------------

library(doParallel)
library(doRNG)
library(extRemes)
library(here)
library(quantreg)
library(rlang)
library(tidyverse)

# Extract command-line arguments ------------------------------------------

# The command-line arguments should be
#   1. the path of the file with the simulation settings
#     This file should create a tibble called sim_settings that has one row for
#     each combination of simulation settings. See files named sim_settings.R in
#     this directory for examples. The columns should be
#     - dataset_type: the type of the data-generating model (DGM), AR or FARIMA
#     - ar_ord: if dataset_type is AR, the order of the AR model that will be
#       used as the DGM. For AR models fit to the data, this will also be their
#       order. If dataset_type is FARIMA, this should be set to one
#     - ar_roots: the roots of the AR polynomial of the AR model. These should
#       lie outside the closed unit disk to ensure that the model is stationary
#       and causal
#     - min_ar_root_abs_val: the minimum absolute value of an AR root
#     - ar_coefs: the AR coefficients of the AR model; these should be created
#       using make_ar_coefs()
#     - farima_ord: if dataset_type is FARIMA, the differencing order of the
#       FARIMA model that will be used as the DGM. It should satisfy
#       alpha * (farima_ord - 1) < -1; see "Fractional ARIMA with stable
#       innovations" by Kokoszka and Taqqu. If dataset_type is AR, this should
#       be NA_real_
#     - sampler: a function that draws a random sample from the innovation
#       distribution, given alpha; should be either rt() or rstable()
#     - alpha: the parameter of the innovation distribution
#     - lead: how far ahead to predict
#     - train_set_size: the size of the training set
#     - test_set_size: the size of the test set
#     - ps: a vector with the levels of the quantiles that will be the
#       thresholds for extremeness
#     - mod_types: a named list with components mod_type and quantile_type.
#       Valid values of mod_type are "naive", "ols", "lad", and "oracle", and
#       valid values of quantile_type are "empirical" and "extreme". If mod_type
#       is "naive", then take_abs_val should also be a component, and it should
#       be TRUE or FALSE
#     - run_num: the number of the simulation run
#   2. the number of cores to use
#     If the number of cores is zero, then all cores will be used.
sim_settings_file <- commandArgs(trailingOnly = TRUE)[1]
num_cores <- as.double(commandArgs(trailingOnly = TRUE)[2])
num_cores <- if_else(num_cores == 0, detectCores(), num_cores)
cat(str_glue("{num_cores} cores are being used"), "\n")

# Source functions --------------------------------------------------------

source(here("simulations/utils.R"))
source(here("code/calc_tail_dep_coef.R"))
source(here("code/mod_fitters.R"))
source(here("code/summarizers.R"))
source(here("code/tail_index_estimators.R"))

# Run the simulations -----------------------------------------------------

source(sim_settings_file)
# Store these separately to save space by avoiding duplicating them
causal_coefs <- sim_settings %>%
  distinct(ar_ord, ar_coefs, farima_ord) %>% # This step speeds up the next step
  mutate(
    ar_causal_coefs = map2(
      ar_ord, ar_coefs, make_ar_causal_coefs, num_terms = 1e6
    ),
    farima_coefs = map(farima_ord, make_farima_coefs, max_lag = 1e6)
  )
oracle_asymp_precisions <- sim_settings %>%
  distinct(ar_ord, ar_coefs, farima_ord, alpha, lead) %>%
  inner_join(causal_coefs) %>%
  mutate(
    oracle_asymp_precision = pmap_dbl(
      list(ar_causal_coefs, farima_coefs, alpha, lead),
      calc_oracle_asymp_precision
    )
  ) %>%
  select(!c(ar_causal_coefs, farima_coefs))
sim_settings <- sim_settings %>%
  inner_join(oracle_asymp_precisions) %>%
  relocate(oracle_asymp_precision, .after = lead)

set.seed(1, "L'Ecuyer-CMRG")

oracle_quantiles <- sim_settings %>%
  distinct(
    dataset_type, ar_ord, ar_coefs, farima_ord, sampler, alpha, lead, ps
  ) %>%
  inner_join(causal_coefs) %>%
  mutate(
    dataset = pmap(
      list(dataset_type, ar_ord, ar_coefs, farima_coefs, sampler, alpha, lead),
      make_dataset,
      train_set_size = 1e6, test_set_size = 0
    ),
    oracle_quantiles = pmap(
      list(dataset_type, ar_coefs, farima_coefs, ps, dataset),
      calc_oracle_quantiles
    )
  ) %>%
  select(!c(ar_causal_coefs, farima_coefs, dataset)) %>%
  unnest(oracle_quantiles)
sim_settings <- sim_settings %>%
  inner_join(oracle_quantiles) %>%
  relocate(oracle_lin_pred_quantiles, oracle_y_quantiles, .after = ps)

registerDoParallel(cores = num_cores)
results <- foreach (i = 1:nrow(sim_settings)) %dorng% {
  cat(str_glue("Iteration {i}"), "\n")
  with(
    sim_settings,
    {
      farima_coefs <- if (!is.na(farima_ord[i])) {
        # Without this, filtering could cause an error
        farima_ord_ <- farima_ord[i]
        causal_coefs %>%
          filter(near(farima_ord, farima_ord_)) %>%
          pull(farima_coefs) %>%
          `[[`(1)
      }
      dataset <- make_dataset(
        dataset_type[i], ar_ord[i], ar_coefs[[i]], farima_coefs,
        sampler[[i]], alpha[i],
        lead[i], train_set_size[i], test_set_size[i]
      )
      foreach (mod_type = mod_types[[i]], .combine = bind_rows) %:%
        foreach (
          p = ps[[i]],
          oracle_lin_pred_quantile = oracle_lin_pred_quantiles[[i]],
          oracle_y_quantile = oracle_y_quantiles[[i]],
          .combine = bind_rows
        ) %do% {
          validate(
            ar_coefs[[i]], dataset,
            p, oracle_lin_pred_quantile, oracle_y_quantile,
            mod_type
          )
        }
    }
  )
} %>%
  add_column(sim_settings, results = .) %>%
  unnest(results) %>%
  select(!c(ps, oracle_lin_pred_quantiles, oracle_y_quantiles, mod_types))

# Save the results --------------------------------------------------------

save(causal_coefs, results, file = "results.RData")
