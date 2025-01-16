# naive_model -------------------------------------------------------------

fit_naive_mod <- function(train_set,
                          q_lin_pred,
                          flux_threshold, p, quantile_type,
                          tail_dep_coef_type, block_size) {
  x <- train_set$flux_lag_0
  y <- train_set$y

  if (is.null(flux_threshold)) {
    if (quantile_type == "empirical") {
      flux_threshold <- quantile(y, p)
    } else {
      flux_threshold <- calc_extreme_quantile(y, p)
    }
  } else {
    # If an oracle model is to be fit, then both flux_threshold and p will be
    # non-null
    if (is.null(p)) p <- mean(y <= flux_threshold)
  }
  if (is.null(q_lin_pred)) {
    if (quantile_type == "empirical") {
      q_lin_pred <- quantile(x, p)
    } else {
      q_lin_pred <- calc_extreme_quantile(x, p)
    }
  }
  tail_dep_coef <- calc_tail_dep_coef(
    x, y, tail_dep_coef_type, q_lin_pred, flux_threshold, block_size
  )

  structure(
    list(
      flux_threshold = flux_threshold, p = p, quantile_type = quantile_type,
      tail_dep_coef_type = tail_dep_coef_type, block_size = block_size,
      # So the test set TDC can be computed; change names here and in
      # predict.naive_mod later
      q_lin_pred = q_lin_pred,
      tail_dep_coef = tail_dep_coef
    ),
    class = "naive_mod"
  )
}

predict.naive_mod <- function(mod, test_set) {
  tibble(
    T_REC = test_set$T_REC,
    y = test_set$y,
    was_y_imputed = test_set$was_y_imputed,
    lin_pred = test_set$flux_lag_0,
    obs = y >= mod$flux_threshold,
    pred = lin_pred >= mod$q_lin_pred
  )
}

new_naive_mod_fitter <- function(flare_class = c("A", "B", "C", "M", "X"),
                                 p = NULL,
                                 quantile_type = c("empirical", "extreme"),
                                 tail_dep_coef_type = c(
                                   "empirical", "modified CFG"
                                 ),
                                 block_size = 1L) {
  stopifnot(xor(is.null(flare_class), is.null(p)))
  match.arg(flare_class)
  if (!is.null(p)) stopifnot(between(p, 0, 1))
  quantile_type <- match.arg(quantile_type)
  tail_dep_coef_type <- match.arg(tail_dep_coef_type)
  stopifnot(is.integer(block_size) && (block_size >= 1))

  if (is.null(p)) {
    flux_threshold <- case_when(
      flare_class == "A" ~ 1e-8,
      flare_class == "B" ~ 1e-7,
      flare_class == "C" ~ 1e-6,
      flare_class == "M" ~ 1e-5,
      flare_class == "X" ~ 1e-4
    )
  } else {
    flux_threshold <- NULL
  }
  mod_fitter <- function(train_set, q_lin_pred = NULL, q_y = NULL, ...) {
    if (!is.null(q_y)) flux_threshold <- q_y
    fit_naive_mod(
      train_set,
      q_lin_pred,
      flux_threshold, p, quantile_type,
      tail_dep_coef_type, block_size
    )
  }
  class(mod_fitter) <- c(class(mod_fitter), "naive_mod_fitter")
  mod_fitter
}

# ar_model ----------------------------------------------------------------

calc_mat_pow <- function(A, pow) {
  if (pow == 0) {
    return(diag(nrow(A)))
  }
  else if (pow == 1) {
    return(A)
  } else {
    B <- calc_mat_pow(A, floor(pow / 2))
    B_squared <- B %*% B
    if (pow %% 2 == 0) {
      return(B_squared)
    } else {
      return(B_squared %*% A)
    }
  }
}

calc_opt_pred_coefs <- function(phi_hat, max_lag, lead) {
  phi_hat %>%
    cbind(diag(nrow = max_lag + 1)[, -(max_lag + 1)]) %>%
    calc_mat_pow(lead) %>%
    `[`(, 1)
}

fit_ar_mod <- function(train_set,
                       phi, q_lin_pred,
                       flux_threshold,
                       p,
                       estim_type, quantile_type,
                       tail_dep_coef_type, block_size) {
  max_lag <- unname(attr(train_set, "attrs_for_mods")["max_lag"])
  lead <- unname(attr(train_set, "attrs_for_mods")["lead"])
  if (is.null(phi)) {
    fit_fun <- if (estim_type == "ols") lm else rq
    phi_hat <- train_set %>%
      select(starts_with("flux")) %>%
      fit_fun(flux_lead_1 ~ . - 1, data = .) %>%
      coef()
  } else {
    phi_hat <- phi
    names(phi_hat) <- train_set %>% names() %>% str_subset("^flux_lag")
  }
  beta_hat <- phi_hat %>%
    calc_opt_pred_coefs(max_lag, lead) %>%
    set_names(names(phi_hat))

  X <- train_set %>% select(starts_with("flux_lag")) %>% as.matrix()
  y <- train_set$y
  if (is.null(flux_threshold)) {
    if (quantile_type == "empirical") {
      flux_threshold <- quantile(y, p)
    } else {
      flux_threshold <- calc_extreme_quantile(y, p)
    }
  } else {
    if (is.null(p)) p <- mean(y <= flux_threshold)
  }
  lin_pred <- drop(X %*% beta_hat)
  if (is.null(q_lin_pred)) {
    if (quantile_type == "empirical") {
      q_lin_pred <- quantile(lin_pred, p)
    } else {
      q_lin_pred <- calc_extreme_quantile(lin_pred, p)
    }
  }
  tail_dep_coef <- calc_tail_dep_coef(
    lin_pred, y, tail_dep_coef_type, q_lin_pred, flux_threshold, block_size
  )

  structure(
    list(
      max_lag = max_lag,
      ord = max_lag + 1,
      lead = lead,
      flux_threshold = flux_threshold,
      p = p,
      estim_type = estim_type,
      tail_dep_coef_type = tail_dep_coef_type, block_size = block_size,
      phi_hat = phi_hat, beta_hat = beta_hat,
      q_lin_pred = q_lin_pred,
      tail_dep_coef = tail_dep_coef
    ),
    class = "ar_mod"
  )
}

predict.ar_mod <- function(mod, test_set) {
  if (is.null(mod$beta_hat)) return(NULL)
  X <- test_set %>% select(starts_with("flux_lag")) %>% as.matrix()
  tibble(
    T_REC = test_set$T_REC,
    y = test_set$y,
    was_y_imputed = test_set$was_y_imputed,
    lin_pred = drop(X %*% mod$beta_hat),
    obs = y >= mod$flux_threshold,
    pred = lin_pred >= mod$q_lin_pred
  )
}

new_ar_mod_fitter <- function(flare_class = c("A", "B", "C", "M", "X"),
                              p = NULL,
                              estim_type = c("ols", "lad"),
                              quantile_type = c("empirical", "extreme"),
                              tail_dep_coef_type = c(
                                "empirical", "modified CFG"
                              ),
                              block_size = 1L) {
  if (!is.null(flare_class)) match.arg(flare_class)
  if (!is.null(p)) stopifnot(between(p, 0, 1))
  match.arg(estim_type)
  match.arg(quantile_type)
  match.arg(tail_dep_coef_type)
  stopifnot(is.integer(block_size) && (block_size >= 1))
  if (!is.null(p)) {
    flux_threshold <- NULL
  } else {
    flux_threshold <- case_when(
      flare_class == "A" ~ 1e-8,
      flare_class == "B" ~ 1e-7,
      flare_class == "C" ~ 1e-6,
      flare_class == "M" ~ 1e-5,
      flare_class == "X" ~ 1e-4
    )
  }

  mod_fitter <- function(train_set, phi = NULL, q_lin_pred = NULL, q_y = NULL) {
    if (is.null(estim_type)) stopifnot(!is.null(phi))
    if (!is.null(estim_type)) stopifnot(is.null(phi))
    if (!is.null(q_y)) flux_threshold <- q_y
    fit_ar_mod(
      train_set, phi, q_lin_pred,
      flux_threshold, p, estim_type, quantile_type,
      tail_dep_coef_type, block_size
    )
  }
  class(mod_fitter) <- c(class(mod_fitter), "ar_mod_fitter")
  mod_fitter
}