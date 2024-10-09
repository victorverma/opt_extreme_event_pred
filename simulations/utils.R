make_ar_roots <- function(ar_ord) {
  num_conj_pairs <- sample(0:floor(ar_ord / 2), 1)
  # Set the moduli and arguments of the complex roots
  r <- rep(rexp(num_conj_pairs, rate = 0.5) + 1, times = 2)
  theta <- runif(num_conj_pairs, max = 2 * pi)
  theta <- c(theta, -theta)
  # Set the moduli and arguments of the real roots
  r <- c(r, rexp(ar_ord - 2 * num_conj_pairs, rate = 0.5) + 1)
  theta <- c(
    theta, sample(c(0, pi), ar_ord - 2 * num_conj_pairs, replace = TRUE)
  )
  complex(modulus = r, argument = theta)
}

check_make_ar_roots <- function() {
  ar_ord <- rpois(1, lambda = 5)
  ar_roots <- make_ar_roots(ar_ord)
  browser(
    expr = any(Mod(ar_roots) <= 1) ||
      !isTRUE(all.equal(sort(ar_roots), sort(Conj(ar_roots))))
  )
}

plot_ar_roots <- function(ar_roots) {
  ar_roots %>%
    {tibble(re = Re(.), im = Im(.))} %>%
    ggplot(aes(re, im)) +
    geom_point() +
    geom_function(fun = ~ sqrt(1 - .x^2), xlim = c(-1, 1), n = 10001) +
    geom_function(fun = ~ -sqrt(1 - .x^2), xlim = c(-1, 1), n = 10001) +
    coord_fixed() +
    labs(x = "Real Part", y = "Imaginary Part") +
    theme_bw()
}

# Make a real vector (coef[1], ... coef[ar_ord]) such that every root of the
# polynomial 1 - coef[1] * z - ... - coef[ar_ord] * z^ar_ord lies outside the
# closed unit disk. Then the corresponding AR(ar_ord) model has a unique 
# stationary causal solution
make_ar_coefs <- function(ar_roots) {
  prod_coefs <- c(-ar_roots[1], 1)
  for (ar_root in ar_roots[-1]) {
    prod_coefs <- c(0, prod_coefs) - c(ar_root * prod_coefs, 0)
  }
  prod_coefs <- suppressWarnings(as.double(prod_coefs))
  lead_coef <- -1 / prod_coefs[1]
  -lead_coef * prod_coefs
}

make_ar_causal_coefs <- function(ar_ord, ar_coefs, num_terms) {
  if (is.null(ar_coefs)) return(NULL)
  ar_coefs <- -ar_coefs[-1]
  ar_causal_coefs <- double(num_terms + ar_ord - 1)
  ar_causal_coefs[ar_ord] <- 1
  for (i in seq_len(num_terms - 1)) {
    ar_causal_coefs[i + ar_ord] <- sum(
      ar_coefs * ar_causal_coefs[(i - 1 + ar_ord):i]
    )
  }
  tail(ar_causal_coefs, num_terms)
}

calc_oracle_asymp_precision <- function(ar_causal_coefs, farima_coefs,
                                        alpha,
                                        lead) {
  if (!is.null(ar_causal_coefs)) {
    coefs <- ar_causal_coefs
  } else {
    coefs <- farima_coefs
  }
  powers <- abs(coefs)^alpha
  # Note that this assumes that the extremal skewness of the innovations
  # distribution equals 0.5
  sum(powers[-(1:lead)]) / sum(powers)
}

make_farima_coefs <- function(d, max_lag) {
  if (is.na(d)) return(NULL)
  # If $\{X_t\}$ is a FARIMA(0, d, 0) process, then
  # $X_t = \sum_{j \ge 0} a_j\epsilon_{t - j}$, where the $\epsilon_t$'s are
  # iid, a_0 = 1, and a_j = $\Gamma(j + d) / \Gamma(d)\Gamma(j + 1)$ for
  # $j \ge 1$. Simplifying that ratio yields the expression below
  c(1, cumprod(1 + (d - 1) / seq_len(max_lag)))
}

check_farima_coefs <- function(d, max_lag) {
  lags <- seq_len(max_lag)
  coefs1 <- c(1, gamma(lags + d) / (gamma(d) * gamma(lags + 1))) %>%
    discard(~ is.nan(.x) || near(.x, 0))
  coefs2 <- make_farima_coefs(d, length(coefs1) - 1)
  all(near(coefs1, coefs2))
}

rt <- function(n, alpha) {
  stats::rt(n, alpha)
}

# Generates independent symmetric alpha-stable random variables with scale
# parameter one, using the Chambers-Mallows-Stuck method as given in Theorem
# 1.3 of Univariate Stable Distributions by John Nolan
rstable <- function(n, alpha) {
  theta <- runif(n, -pi / 2, pi / 2)
  w <- rexp(n, 1)
  if (alpha != 1) {
    z <- (sin(alpha * theta) / (cos(theta)^(1 / alpha))) *
      (cos((alpha - 1) * theta) / w)^((1 - alpha) / alpha)
  } else {
    z <- tan(theta)
  }
  z
}

make_dataset <- function(dataset_type = c("AR", "FARIMA"),
                         ar_ord, ar_coefs, farima_coefs, sampler, alpha,
                         lead, train_set_size, test_set_size) {
  match.arg(dataset_type)
  if (is.null(ar_coefs) && is.null(farima_coefs)) return(NULL)

  if (dataset_type == "AR") {
    vals <- arima.sim(
      model = list(ar = -ar_coefs[-1]),
      n = ar_ord + lead + train_set_size + test_set_size - 1,
      rand.gen = sampler,
      alpha = alpha
    )    
  } else {
    num_farima_coefs <- length(farima_coefs)
    vals <- sampler(
      n = ar_ord + lead + train_set_size + test_set_size + num_farima_coefs - 2,
      alpha
    ) %>%
      # Equivalent to
      # stats::filter(farima_coefs, sides = 1) %>%
      #   as.vector() %>%
      #   tail(length(.) - num_farima_coefs + 1),
      # but faster since it uses the FFT
      convolve(rev(farima_coefs), type = "filter")
  }
  
  dataset <- vals %>%
    embed(ar_ord) %>%
    # Column names need to start with "flux_lag_" so that sourced code will work
    `colnames<-`(str_c("flux_lag", 0:(ar_ord - 1), sep = "_")) %>%
    as_tibble() %>%
    # The value at a lead time of one will be needed to fit the AR model
    mutate(flux_lead_1 = lead(flux_lag_0, 1), y = lead(flux_lag_0, lead)) %>%
    slice_head(n = train_set_size + test_set_size) %>%
    mutate(
      in_train_set = row_number() <= train_set_size,
      # T_REC and was_y_imputed need to be added to make sourced code work
      T_REC = NA_real_,
      was_y_imputed = FALSE
    )
  # max_lag and lead are needed for sourced code to work
  attributes(dataset) <- c(
    attributes(dataset),
    list(attrs_for_mods = c(max_lag = ar_ord - 1, lead = lead))
  )
  dataset
}

plot_dataset <- function(dataset_type = c("AR", "FARIMA"),
                         ar_ord, farima_ord, alpha,
                         dataset) {
  match.arg(dataset_type)
  if (dataset_type == "AR") {
    title_prefix <- str_glue("AR({ar_ord}) Model")
  } else {
    title_prefix <- str_glue("FARIMA(0, {farima_ord}, 0) Model")
  }
  
  dataset %>%
    mutate(t = row_number()) %>%
    ggplot(aes(t, flux_lag_0, color = in_train_set)) +
    geom_line() +
    scale_color_manual(values = c("orange", "blue")) +
    labs(
      y = expression(Y[t]), color = "In Training Set?",
      title = bquote(list(.(title_prefix), alpha == .(alpha)))
    ) +
    theme_bw() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
}

help_calc_oracle_quantiles <- function(dataset_type, ar_coefs, p, dataset) {
  if (dataset_type == "AR") {
    mod_fitter <- new_ar_mod_fitter(
      p = p,
      estim_type = NULL,
      quantile_type = "empirical",
      tail_dep_coef_type = "empirical"
    )
  } else {
    mod_fitter <- new_naive_mod_fitter(flare_class = NULL, p = p)
  }

  train_set <- dataset %>% filter(in_train_set) %>% select(!in_train_set)
  attr(train_set, "attrs_for_mods") <- attr(dataset, "attrs_for_mods")
  mod <- mod_fitter(train_set, phi = -ar_coefs[-1])
  tibble(
    oracle_y_quantile = mod$flux_threshold,
    oracle_lin_pred_quantile = mod$q_lin_pred
  )
}

calc_oracle_quantiles <- function(dataset_type,
                                  ar_coefs, farima_coefs,
                                  ps,
                                  dataset) {
  foreach (p = ps, .combine = bind_rows) %do% {
    help_calc_oracle_quantiles(dataset_type, ar_coefs, p, dataset)
  } %>%
    summarize(
      oracle_lin_pred_quantiles = list(oracle_lin_pred_quantile),
      oracle_y_quantiles = list(oracle_y_quantile)
    )
}

validate <- function(ar_coefs, dataset,
                     p, oracle_lin_pred_quantile, oracle_y_quantile,
                     mod_type) {
  if (mod_type$mod_type == "naive") {
    mod_fitter <- new_naive_mod_fitter(
      flare_class = NULL, p = p, quantile_type = mod_type$quantile_type
    )
    oracle_lin_pred_quantile <- oracle_y_quantile
  } else {
    mod_fitter <- new_ar_mod_fitter(
      p = p,
      estim_type = if (!mod_type$is_oracle) mod_type$mod_type,
      quantile_type = mod_type$quantile_type,
      tail_dep_coef_type = "empirical"
    )
  }

  train_set <- dataset %>% filter(in_train_set) %>% select(!in_train_set)
  attr(train_set, "attrs_for_mods") <- attr(dataset, "attrs_for_mods")
  test_set <- dataset %>% filter(!in_train_set) %>% select(!in_train_set)
  mod <- mod_fitter(
    train_set,
    phi = if (mod_type$is_oracle) -ar_coefs[-1],
    q_lin_pred = if (mod_type$is_oracle) oracle_lin_pred_quantile,
    q_y = oracle_y_quantile
  )
  confusion_mat <- mod %>% predict(test_set) %>% make_confusion_mat()

  tibble(
    p,
    oracle_lin_pred_quantile, oracle_y_quantile,
    mod_type = list(mod_type),
    flux_threshold = mod$flux_threshold,
    q_lin_pred = mod$q_lin_pred,
    confusion_mat = list(confusion_mat),
    train_precision = mod$tail_dep_coef,
    test_precision = calc_precision(confusion_mat[[1]]),
    tpr = calc_tpr(confusion_mat[[1]]), fpr = calc_fpr(confusion_mat[[1]])
  )
}
