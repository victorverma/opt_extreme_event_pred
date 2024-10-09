here::i_am("simulations/farima_asymp_precisions/farima_asymp_precisions.R")

# Load packages -----------------------------------------------------------

library(doParallel)
library(here)
library(tidyverse)

# Source functions --------------------------------------------------------

source(here("simulations/utils.R"))

# Calculate the precisions ------------------------------------------------

# The command-line arguments should be
# 1. The number of alpha values to use
# 2. The number of d values to use
# 3. The number of lead times to use
# 4. The number of cores to use
#   If the number of cores is zero, then all cores will be used.
# For each alpha, d, and lead time, the asymptotic precision at the lead time of
# the optimal predictor for a FARIMA(0, d, 0) model with symmetric alpha-stable,
# unit scale parameter innovations will be calculated
cmd_args <- commandArgs(trailingOnly = TRUE)
num_alphas <- as.integer(cmd_args[1])
num_ds <- as.integer(cmd_args[2])
num_leads <- as.integer(cmd_args[3])
num_cores <- as.integer(cmd_args[4])
num_cores <- if_else(num_cores == 0, detectCores(), num_cores)
cat(str_glue("Number of cores in use: {num_cores}"), "\n")
registerDoParallel(cores = num_cores)

settings <- expand_grid(
  alpha = seq(1, 2, length.out = num_alphas + 2)[-c(1, num_alphas + 2)],
  d = seq(0, 0.5, length.out = num_ds + 2)[-c(1, num_ds + 2)],
  leads = list(seq_len(num_leads))
) %>%
  filter(alpha * (d - 1) < -1)

farima_coefs <- settings %>%
  distinct(d) %>%
  mutate(
    farima_coefs = mclapply(
      d, make_farima_coefs, max_lag = 1e6, mc.cores = num_cores
    ),
    last_to_first_coef_ratio = map_dbl(farima_coefs, ~ last(.x) / first(.x)),
    are_coefs_strictly_decreasing = map_lgl(farima_coefs, ~ all(diff(.x) < 0))
  )

results <- foreach (i = 1:nrow(settings)) %dopar% {
  with(
    settings,
    {
      d_ <- d[i]
      farima_coefs_row <- filter(farima_coefs, near(d_, d))
      tibble(
        last_to_first_coef_ratio = farima_coefs_row$last_to_first_coef_ratio,
        are_coefs_strictly_decreasing =
          farima_coefs_row$are_coefs_strictly_decreasing,
        lead = leads[[i]],
        oracle_asymp_precision = map_dbl(
          lead,
          calc_oracle_asymp_precision,
          ar_causal_coefs = NULL,
          farima_coefs = farima_coefs_row$farima_coefs[[1]],
          alpha = alpha[i]
        )
      )
    }
  )
} %>%
  add_column(settings, results = .) %>%
  select(!leads) %>%
  unnest(results)

# Save the results --------------------------------------------------------

save(results, file = here("simulations/farima_asymp_precisions/results.RData"))