# See data_analysis/comments_on_utils.html for comments on the code below

# Calculate coefficients --------------------------------------------------

cppFunction('
  NumericVector calc_causal_coefs(double d, int num_coefs) {
    NumericVector causal_coefs(num_coefs);
    causal_coefs[0] = 1;
    for (int j = 1; j < num_coefs; j++) {
      causal_coefs[j] = causal_coefs[j - 1] * (1 + (d - 1) / j);
    }
    return causal_coefs;
  }
')

cppFunction('
  NumericVector calc_invert_coefs(NumericVector causal_coefs,
                                  int num_coefs) {
    NumericVector invert_coefs(num_coefs);
    invert_coefs[0] = 1;
    for (int j = 1; j < num_coefs; j++) {
      for (int k = 1; k <= j; k++) {
        invert_coefs[j] -= causal_coefs[k] * invert_coefs[j - k];
      }
    }
    return invert_coefs;
  }
')

# Try implementing this in C++?
calc_opt_pred_coefs <- function(causal_coefs, invert_coefs, num_coefs,
                                num_obs, lead) {
  stopifnot(num_coefs >= num_obs + lead)
  coefs <- double(num_obs)
  for (r in 0:(num_obs - 1)) {
    coefs[r + 1] <- sum(
      causal_coefs[(lead + 1):(r + lead + 1)] * invert_coefs[(r + 1):1]
    )
  }
  coefs
}

# Work with FARIMA processes ----------------------------------------------

sim_from_farima <- function(d, num_coefs, innov_sampler, n, num_burn_in_obs) {
  causal_coefs <- calc_causal_coefs(d, num_coefs)
  innov_sampler(n + num_burn_in_obs + num_coefs - 1) %>%
    convolve(rev(causal_coefs), type = "filter") %>%
    tail(n)
}

estim_innovs <- function(x, invert_coefs) {
  convolve(x, rev(invert_coefs), type = "filter")
}

suggest_num_coefs <- function(d, alpha, epsilon = 1e-3) {
  ceiling(epsilon^(1 / (alpha * (d - 1) + 1)))
}

# Work with skewed t distributions ----------------------------------------

dt <- function(x, alpha, mu, sigma, xi, log = FALSE) {
  stopifnot(alpha > 0, sigma > 0, xi > 0)
  norm_const <- 2 * xi / (sigma * (xi^2 + 1))
  multipliers <- if_else(x < mu, xi, 1 / xi)
  vals <- stats::dt(multipliers * (x - mu) / sigma, df = alpha, log = log)
  if (log) {
    vals <- log(norm_const) + vals
  } else {
    vals <- norm_const * vals
  }
  vals
}

pt <- function(q, alpha, mu, sigma, xi) {
  stopifnot(alpha > 0, sigma > 0, xi > 0)
  multipliers <- if_else(q < mu, xi, 1 / xi)
  base_probs <- stats::pt(multipliers * (q - mu) / sigma, df = alpha)
  if_else(
    q < mu,
    2 * base_probs / (xi^2 + 1),
    (1 - xi^2 + 2 * xi^2 * base_probs) / (xi^2 + 1)
  )
}

rt <- function(n, alpha, mu, sigma, xi) {
  b <- rbinom(n, size = 1, prob = xi^2 / (xi^2 + 1))
  w <- stats::rt(n, df = alpha)
  z <- ((xi + 1 / xi) * b - 1 / xi) * abs(w)
  sigma * z + mu
}

calc_t_mean <- function(alpha, mu, sigma, xi) {
  stopifnot(alpha > 1, sigma > 0, xi > 0)
  alpha_expr <- 2 * gamma((alpha + 1) / 2)
  alpha_expr <- alpha_expr / ((alpha - 1) * gamma(alpha / 2))
  alpha_expr <- alpha_expr * sqrt(alpha / pi)
  xi_expr <- (xi^4 - 1) / (xi * (xi^2 + 1))
  alpha_expr * sigma * xi_expr + mu
}

calc_t_mu <- function(alpha, sigma, xi) {
  stopifnot(alpha > 1, sigma > 0, xi > 0)
  -calc_t_mean(alpha, mu = 0, sigma, xi)
}

calc_t_log_lik <- function(theta, x, num_coefs) {
  rho <- theta[1]
  alpha <- theta[2]
  sigma <- theta[3]
  xi <- theta[4]
  stopifnot(rho > 0, rho < 1, alpha > 1, sigma > 0, xi > 0)
  mu <- calc_t_mu(alpha, sigma, xi)
  d <- rho * (1 - 1 / alpha) # Denominator should really be min(alpha, 2)?
  causal_coefs <- calc_causal_coefs(d, num_coefs)
  invert_coefs <- calc_invert_coefs(causal_coefs, num_coefs)
  innov_hats <- estim_innovs(x, invert_coefs)
  sum(dt(innov_hats, alpha, mu, sigma, xi, log = TRUE))
}

calc_t_penalty <- function(theta, hyper_theta) {
  sigma <- theta[3]
  xi <- theta[4]
  zeta <- hyper_theta[1]
  eta <- hyper_theta[2]
  log_norm_const <- zeta * log(eta) - log(gamma(zeta))
  sigma_penalty <- log_norm_const - (zeta + 1) * log(sigma) - eta / sigma
  xi_penalty <- log_norm_const - (zeta + 1) * log(xi) - eta / xi
  sigma_penalty + xi_penalty
}

calc_t_obj_fun <- function(x, num_coefs, hyper_theta) {
  force(x)
  force(num_coefs)
  if (is.null(hyper_theta)) {
    obj_fun <- partial(calc_t_log_lik, x = x, num_coefs = num_coefs)
  } else {
    obj_fun <- function(theta) {
      calc_t_log_lik(theta, x, num_coefs) + calc_t_penalty(theta, hyper_theta)
    }
  }
  obj_fun
}

calc_t_mle <- function(x, num_coefs, hyper_theta = NULL, ci_level = 0.95,
                       optimizer, par, ...) {
  obj_fun <- calc_t_obj_fun(x, num_coefs, hyper_theta)
  names(par) <- c("rho", "alpha", "sigma", "xi")
  results <- optimizer(obj_fun, par, ...)
  results$cis <- calc_cis(results$par, results$hessian, ci_level)
  results
}

# Work with Pareto distributions ------------------------------------------

cppFunction('
  NumericVector dpareto(NumericVector x,
                        double alpha, double mu, double sigma,
                        bool log_ = false) {
    if (alpha <= 0)
      stop("alpha must be positive");
    if (sigma <= 0)
      stop("sigma must be positive");
    int n = x.size();
    double support_min = mu + sigma;
    double norm_const = alpha / sigma;
    double log_norm_const = log(norm_const);
    double power = -alpha - 1;
    NumericVector vals(n);
    if (log_) {
      for (int i = 0; i < n; i++) {
        if (x[i] >= support_min)
          vals[i] = log_norm_const + power * log((x[i] - mu) / sigma);
        else
          vals[i] = R_NegInf;
      }
    } else {
      for (int i = 0; i < n; i++) {
        if (x[i] >= support_min)
          vals[i] = norm_const * pow((x[i] - mu) / sigma, power);
        else
          vals[i] = 0;
      }
    }
    return vals;
  }
')

ppareto <- function(q, alpha, mu, sigma) {
  stopifnot(alpha > 0, sigma > 0)
  if_else(q >= mu + sigma, 1 - ((q - mu) / sigma)^(-alpha), 0)
}

rpareto <- function(n, alpha, mu, sigma) {
  u <- runif(n)
  sigma * u^(-1 / alpha) + mu
}

calc_pareto_mean <- function(alpha, mu, sigma) {
  stopifnot(alpha > 1, sigma > 0)
  (alpha * sigma) / (alpha - 1) + mu
}

calc_pareto_mu <- function(alpha, sigma) {
  stopifnot(alpha > 1, sigma > 0)
  -calc_pareto_mean(alpha, mu = 0, sigma)
}

calc_pareto_log_lik <- function(theta, x, num_coefs) {
  rho <- theta[1]
  alpha <- theta[2]
  sigma <- theta[3]
  stopifnot(rho > 0, rho < 1, alpha > 1, sigma > 0)
  mu <- calc_pareto_mu(alpha, sigma)
  d <- rho * (1 - 1 / alpha)
  causal_coefs <- calc_causal_coefs(d, num_coefs)
  invert_coefs <- calc_invert_coefs(causal_coefs, num_coefs)
  innov_hats <- estim_innovs(x, invert_coefs)
  sum(dpareto(innov_hats, alpha, mu, sigma, log = TRUE))
}

calc_pareto_penalty <- function(theta, hyper_theta) {
  sigma <- theta[3]
  zeta <- hyper_theta[1]
  eta <- hyper_theta[2]
  log_norm_const <- zeta * log(eta) - log(gamma(zeta))
  log_norm_const - (zeta + 1) * log(sigma) - eta / sigma
}

calc_pareto_obj_fun <- function(x, num_coefs, hyper_theta) {
  force(x)
  force(num_coefs)
  if (is.null(hyper_theta)) {
    obj_fun <- partial(calc_pareto_log_lik, x = x, num_coefs = num_coefs)
  } else {
    obj_fun <- function(theta) {
      calc_pareto_log_lik(theta, x, num_coefs) +
        calc_pareto_penalty(theta, hyper_theta)
    }
  }
  obj_fun
}

calc_pareto_mle <- function(x, num_coefs, hyper_theta = NULL, ci_level = 0.95,
                            optimizer, par, ...) {
  obj_fun <- calc_pareto_obj_fun(x, num_coefs, hyper_theta)
  names(par) <- c("rho", "alpha", "sigma")
  results <- optimizer(obj_fun, par, ...)
  results$cis <- calc_cis(results$par, results$hessian, ci_level)
  results
}

# Run gradient ascent -----------------------------------------------------

run_line_search <- function(obj_fun,
                            par, val, grad,
                            step_sizes, step_size_multipliers,
                            lower, upper) { # Limit number of iterations
  next_par <- par + step_sizes * grad
  is_out_of_bounds <- !between(next_par, lower, upper)
  while (any(is_out_of_bounds)) {
    step_sizes <- step_size_multipliers^is_out_of_bounds * step_sizes
    next_par <- par + step_sizes * grad
    is_out_of_bounds <- !between(next_par, lower, upper)
  }
  next_val <- coalesce(obj_fun(next_par), -Inf)
  while (next_val < val) {
    step_sizes <- step_size_multipliers * step_sizes
    next_par <- par + step_sizes * grad
    next_val <- coalesce(obj_fun(next_par), -Inf)
  }
  list(next_par = next_par, next_val = next_val)
}

run_grad_ascent <- function(obj_fun,
                            par,
                            lower, upper,
                            step_sizes, step_size_multipliers,
                            max_num_iters,
                            trace = FALSE,
                            tol = .Machine$double.eps^0.5,
                            hessian = FALSE) {
  if (!all(between(par, lower, upper)))
    stop("initial value is outside bounds")
  val <- obj_fun(par)
  if (!is.finite(val)) stop("objective function is not finite at initial value")

  iter_num <- 0
  print(str_glue("iter {iter_num}"))
  start_time <- proc.time()
  if (is.null(names(par))) {
    param_names <- str_c("par", 1:length(par))
  } else {
    param_names <- names(par)
  }
  if (trace) {
    trace_mat <- matrix(nrow = max_num_iters + 1, ncol = length(par) + 1)
    colnames(trace_mat) <- c(param_names, "obj_fun_val")
    trace_mat[1, ] <- c(par, val)
  }
  print(str_glue("par = ({toString(par)})", "val = {val}", .sep = "\n"))
  print(str_glue("time (s) = {(proc.time() - start_time)['elapsed']}"))
  prev_par <- par
  prev_val <- val

  while (iter_num < max_num_iters) {
    iter_num <- iter_num + 1
    print(str_glue("iter {iter_num}"))
    grad <- tryCatch(
      grad(obj_fun, par, method = "simple"), error = function(e) NULL
    )
    if (is.null(grad) || !all(is.finite(grad))) {
      warning("returning as grad is NULL or not all finite", immediate. = TRUE)
      break
    }
    with(
      run_line_search(
        obj_fun, par, val, grad, step_sizes, step_size_multipliers, lower, upper
      ),
      {par <<- next_par; val <<- next_val}
    )
    if (trace) trace_mat[iter_num + 1, ] <- c(par, val)
    print(str_glue("par = ({toString(par)})", "val = {val}", .sep = "\n"))
    print(str_glue("time (s) = {(proc.time() - start_time)['elapsed']}"))
    if (abs(val / prev_val - 1) < tol) break
    prev_par <- par
    prev_val <- val
  }

  results <- list(par = set_names(par, param_names), obj_fun_val = val)
  if (trace) results$trace_mat <- head(trace_mat, iter_num + 1)
  if (hessian) {
    cat("Calculating final Hessian", "\n")
    results$hessian <- tryCatch(
      structure(hessian(obj_fun, par), dimnames = rep(list(param_names), 2)),
      error = function(e) NULL
    )
    if (is.null(results$hessian) || !all(is.finite(results$hessian))) {
      warning(
        "Hessian can't be computed or isn't all finite", immediate. = TRUE
      )
    }
    print(str_glue("time (s) = {(proc.time() - start_time)['elapsed']}"))
  }
  results
}

run_lbfgsb <- function(obj_fun,
                       par,
                       lower, upper,
                       control = list(),
                       hessian = FALSE) {
  if (!all(between(par, lower, upper)))
    stop("initial value is outside bounds")
  val <- obj_fun(par)
  if (!is.finite(val)) stop("objective function is not finite at initial value")

  if (is.null(control$fnscale)) control$fnscale <- -1
  if (control$fnscale >= 0) {
    control$fnscale <- -1
    warning("control$fnscale corrected to -1", immediate. = TRUE)
  }

  cat("Running L-BFGS-B", "\n")
  start_time <- proc.time()
  if (is.null(names(par))) {
    param_names <- str_c("par", 1:length(par))
  } else {
    param_names <- names(par)
  }
  results_ <- tryCatch(
    optim(
      par = par,
      fn = obj_fun,
      method = "L-BFGS-B",
      lower = lower, upper = upper,
      control = control,
      hessian = hessian
    ),
    error = function(e) toString(e)
  )
  print(str_glue("time (s) = {(proc.time() - start_time)['elapsed']}"))

  if (is.character(results_)) {
    warning(
      str_glue("L-BFGS-B failed with this error:\n{results_}"),
      immediate. = TRUE
    )
    results <- NULL
  } else {
    results <- with(
      results_,
      list(
        par = set_names(par, param_names),
        obj_fun_val = value,
        diagnostics = list(
          counts = counts, convergence = convergence, message = message
        ),
        hessian = tryCatch(
          structure(hessian, dimnames = rep(list(param_names), 2)),
          error = function(e) NULL
        )
      )
    )
    if (is.null(results$hessian) || !all(is.finite(results$hessian))) {
      warning(
        "Hessian can't be computed or isn't all finite", immediate. = TRUE
      )
    }
  }
  results
}

# Calculate log-likelihoods -----------------------------------------------

calc_cis <- function(par, hessian, ci_level = 0.95) {
  if (is.null(hessian) || !all(is.finite(hessian))) {
    warning(
      "couldn't calculate CIs as Hessian is NULL or not all finite",
      immediate. = TRUE
    )
    return(NULL)
  }
  hessian_inv <- solve(hessian)
  std_errs <- sqrt(diag(-hessian_inv))
  q <- qnorm((1 + ci_level) / 2)

  cis <- par %>%
    enframe(name = "param", value = "est") %>%
    mutate(ci_left = est - q * std_errs, ci_right = est + q * std_errs) %>%
    relocate(est, .after = ci_left)

  rho <- par["rho"]
  alpha <- par["alpha"]
  d <- rho * (1 - 1 / alpha)
  d_grad <- c(1 - 1 / alpha, rho / alpha^2, 0, 0)
  d_std_err <- sqrt(drop(t(d_grad) %*% -hessian_inv %*% d_grad))
  d_ci <- d + c(-1, 1) * q * d_std_err
  d_ci %>%
    set_names(c("ci_left", "ci_right")) %>%
    as_tibble_row() %>%
    add_column(param = "d", .before = 1) %>%
    add_column(est = d, .after = "ci_left") %>%
    bind_rows(cis, .)
}

# Do EDA ------------------------------------------------------------------

make_run_tbl <- function(dataset) {
  dataset %>%
    pull(flux) %>%
    is.na() %>%
    rle() %>%
    with(
      tibble(
        is_na_run = values,
        run_length = lengths,
        start_index = cumsum(c(1, head(run_length, -1))),
        end_index = start_index + run_length - 1,
        start_time = dataset$time[start_index],
        end_time = dataset$time[end_index]
      )
    )
}

make_minutely_window_plot <- function(time_, radius) {
  time_ <- as_datetime(time_)
  time_str <- format(time_, format = "%Y-%m-%d %H:%M:%S", usetz = TRUE)
  radius_str <- case_when(
    radius@year > 0 ~ str_glue("{radius@year} year(s)"),
    radius@month > 0 ~ str_glue("{radius@month} month(s)"),
    radius@day > 0 ~ str_glue("{radius@day} day(s)"),
    radius@hour > 0 ~ str_glue("{radius@hour} hour(s)"),
    radius@minute > 0 ~ str_glue("{radius@minute} minute(s)")
  )
  flux_tbl %>%
    filter(between(time, time_ - radius, time_ + radius)) %>%
    ggplot(aes(time, flux)) +
    geom_point(alpha = 0.5, shape = "bullet") +
    scale_y_log10() +
    labs(
      x = "Time", y = expression("X-Ray Flux" ~~ (W / m^2)),
      title = bquote(.(time_str) %+-% .(radius_str))
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}

make_window_plot <- function(dataset, time_, radius) {
  units <- units(diff(dataset$time[1:2]))
  y_desc <- case_when(
    units == "days" ~ "Daily Maximum X-Ray Flux",
    units == "hours" ~ "Hourly Maximum X-Ray Flux",
    .default = "Maximum X-Ray Flux"
  )
  time_str <- format(time_,  usetz = TRUE)
  radius_str <- case_when(
    radius@year > 0 ~ str_glue("{radius@year} year(s)"),
    radius@month > 0 ~ str_glue("{radius@month} month(s)"),
    radius@day > 0 ~ str_glue("{radius@day} day(s)"),
    radius@hour > 0 ~ str_glue("{radius@hour} hour(s)"),
    radius@minute > 0 ~ str_glue("{radius@minute} minute(s)")
  )
  dataset %>%
    filter(between(time, time_ - radius, time_ + radius)) %>%
    ggplot(aes(time, flux, color = was_interpolated, size = was_interpolated)) +
    geom_point() +
    geom_line(aes(time, flux), inherit.aes = FALSE) +
    scale_y_log10() +
    scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black")) +
    labs(
      x = "Time", y = bquote(.(y_desc) ~~ (W / m^2)),
      color = "Was Interpolated?",
      size = "Was Interpolated?",
      title = bquote(.(time_str) %+-% .(radius_str))
    ) +
    theme_bw() +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5))
}

# Use the Hill estimator --------------------------------------------------

calc_hill_ests <- function(x) {
  x <- sort(x, decreasing = TRUE)
  log_x <- log(x)
  leaded_log_x <- lead(log_x, 1)
  tibble(
    num_top_ord_stats = seq_along(x),
    hill_est = 1 / (cummean(log_x) - leaded_log_x)
  ) %>%
    slice_head(n = -1)
}

make_hill_plot <- function(dataset,
                           first_train_time, last_train_time, label,
                           max_num_top_ord_stats = 300,
                           alpha_hat = NULL) {
  plot <- dataset %>%
    filter(between(time, first_train_time, last_train_time)) %>%
    pull(flux) %>%
    calc_hill_ests() %>%
    filter(num_top_ord_stats <= max_num_top_ord_stats) %>%
    ggplot(aes(num_top_ord_stats, hill_est)) +
    geom_line() +
    labs(
      x = "Number of Top Order Statistics",
      y = expression("Hill Estimate of " * alpha)
    ) +
    ggtitle(label) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  if (!is.null(alpha_hat)) {
    plot <- plot +
      geom_hline(linetype = "dashed", yintercept = alpha_hat) +
      labs(caption = bquote("Dashed line at " * hat(alpha) == .(alpha_hat)))
  }
  plot
}

# Fit FARIMA models -------------------------------------------------------

cppFunction('
  NumericVector calc_integrand(NumericVector x,
                               double d,
                               NumericVector lambda) {
    int n = x.size();
    int num_lambdas = lambda.size();
    NumericVector integrand_vals(num_lambdas);
    double cos_term, sin_term;
    for (int j = 0; j < num_lambdas; j++) {
      cos_term = 0.0;
      sin_term = 0.0;
      for (int t = 0; t < n; t++) {
        cos_term += x[t] * cos(lambda[j] * (t + 1));
        sin_term += x[t] * sin(lambda[j] * (t + 1));
      }
      integrand_vals[j] = pow(2 - 2 * cos(lambda[j]), d) *
        (pow(cos_term, 2) + pow(sin_term, 2));
    }
    return integrand_vals;
  }
')

calc_integral <- function(x, d, subdivisions = 10000) {
  tryCatch(
    integrate(
      f = calc_integrand,
      lower = 1 / length(x), upper = pi,
      x = x, d = d,
      subdivisions = subdivisions
    )$value,
    error = function(e) NA_real_
  )
}

make_obj_fun_plot <- function(obj_fun, x, from, to, d_hat, label) {
  tibble(
    d = seq(from, to, length.out = 100),
    obj_fun_val = map_dbl(d, obj_fun, x = x)
  ) %>%
    ggplot(aes(d, obj_fun_val)) +
    geom_line() +
    geom_vline(linetype = "dashed", xintercept = d_hat) +
    scale_y_log10() +
    labs(
      y = "Objective Function Value",
      title = label,
      caption = bquote("Dashed line at " * hat(d) == .(d_hat))
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}

calc_alpha <- function(x) {
  tryCatch(
    x %>% fevd(type = "GEV") %>% distill() %>% `[`("shape") %>% `/`(1, .),
    error = function(e) NA_real_
  )
}

calc_d <- function(x, alpha_hat, plot_obj_fun = FALSE, label = NULL) {
  tryCatch(
    {
      d_hat <- optimize(
        f = calc_integral, interval = c(-0.5, 1 - 1 / alpha_hat), x = x
      )$minimum
      if (plot_obj_fun) {
        obj_fun_plot <- make_obj_fun_plot(
          calc_integral, x, from = -0.5, to = 1 - 1 / alpha_hat, d_hat, label
        )
      } else {
        obj_fun_plot <- NULL
      }
      tibble(d_hat, obj_fun_plot = list(obj_fun_plot))
    },
    error = function(e) tibble(d_hat = NA_real_, obj_fun_plot = list(NULL))
  )

}

make_data_plot <- function(dataset, first_train_time, last_train_time, label) {
  dataset %>%
    filter(between(time, first_train_time, last_train_time)) %>%
    ggplot(aes(time, flux)) +
    geom_line() +
    scale_y_log10() +
    labs(x = "Time", y = expression("X-Ray Flux" ~~ (W / m^2)), title = label) +
    theme_bw()
}

calc_asymp_opt_precision <- function(alpha, causal_coefs, lead) {
  powers <- abs(causal_coefs)^alpha
  # Note that this assumes that the extremal skewness of the innovations is 0.5
  sum(powers[-(1:lead)]) / sum(powers)
}

help_calc_farima_preds <- function(thresholds_levels, ords, leads) {
  expand_grid(ord = ords, lead = leads, thresholds_levels) %>%
    add_column(
      asymp_opt_precision = NA_real_,
      opt_pred_coefs = list(NULL),
      lin_pred_threshold = NA_real_,
      obs = NA,
      pred = NA
    )
}

calc_farima_preds <- function(x, y,
                              alpha_hat, d_hat,
                              num_coefs,
                              thresholds_levels, ords, leads) {
  if (is.na(alpha_hat) || is.na(d_hat)) {
    return(help_calc_farima_preds(thresholds_levels, ords, leads))
  }

  causal_coefs <- calc_causal_coefs(d_hat, num_coefs)
  invert_coefs <- calc_invert_coefs(causal_coefs, num_coefs)

  foreach (ord = ords, .combine = bind_rows) %do% {
    X <- embed(x, ord)
    foreach (lead = leads, .combine = bind_rows) %do% {
      asymp_opt_precision <- calc_asymp_opt_precision(
        alpha_hat, causal_coefs, lead
      )
      opt_pred_coefs <- calc_opt_pred_coefs(
        causal_coefs, invert_coefs, num_coefs, ord, lead
      )
      lin_preds <- drop(X %*% opt_pred_coefs)
      foreach (i = 1:nrow(thresholds_levels), .combine = bind_rows) %do% {
        lin_pred_threshold <- estim_quantile(
          lin_preds, thresholds_levels$quantile_level[i]
        )
        tibble(
          ord, lead,
          type = thresholds_levels$type[i],
          flux_threshold = thresholds_levels$flux_threshold[i],
          quantile_level = thresholds_levels$quantile_level[i],
          asymp_opt_precision,
          opt_pred_coefs = list(opt_pred_coefs),
          lin_pred_threshold,
          obs = y[leads == lead] >= flux_threshold,
          pred = last(lin_preds) >= lin_pred_threshold
        )
      }
    }
  }
}

fit_farima_mod <- function(dataset,
                           first_train_time, last_train_time,
                           alpha_hat, d_hat,
                           num_coefs,
                           flux_thresholds, quantile_levels, ords, leads,
                           plot_obj_fun = FALSE, label = NULL) {
  x <- dataset %>%
    filter(between(time, first_train_time, last_train_time)) %>%
    pull(flux)
  x_mean <- mean(x)
  if (is.null(alpha_hat)) alpha_hat <- calc_alpha(x - x_mean)
  if (is.null(d_hat)) {
    result <- calc_d(x - x_mean, alpha_hat, plot_obj_fun, label)
  } else {
    result <- tibble(d_hat, obj_fun_plot = list(NULL))
  }
  i <- which(dataset$time == last_train_time)
  y <- dataset$flux[i + leads]
  thresholds_levels <- make_thresholds_levels(
    x, flux_thresholds, quantile_levels
  )
  x <- x - x_mean
  preds <- calc_farima_preds(
    x, y,
    alpha_hat, result$d_hat,
    num_coefs,
    thresholds_levels, ords, leads
  )
  tibble(alpha_hat, result, preds = list(preds))
}

# Fit AR models -----------------------------------------------------------

calc_ar_preds <- function(x, y, loss_type, thresholds_levels, ords, leads) {
  n <- length(x)
  match.arg(loss_type, c("ols", "lad"))
  mod_fitter <- if (loss_type == "ols") lm else rq

  foreach (ord = ords, .combine = bind_rows) %do% {
    X <- embed(x, ord)
    foreach (lead = leads, .combine = bind_rows) %do% {
      z <- x[(ord + lead):n]
      opt_pred_coefs <- unname(
        coef(mod_fitter(z ~ head(X, n - ord - lead + 1) - 1))
      )
      lin_preds <- drop(X %*% opt_pred_coefs)
      foreach (i = 1:nrow(thresholds_levels), .combine = bind_rows) %do% {
        lin_pred_threshold <- estim_quantile(
          lin_preds, thresholds_levels$quantile_level[i]
        )
        tibble(
          ord, lead,
          type = thresholds_levels$type[i],
          flux_threshold = thresholds_levels$flux_threshold[i],
          quantile_level = thresholds_levels$quantile_level[i],
          opt_pred_coefs = list(opt_pred_coefs),
          lin_pred_threshold,
          obs = y[leads == lead] >= flux_threshold,
          pred = last(lin_preds) >= lin_pred_threshold
        )
      }
    }
  }
}

fit_ar_mod <- function(dataset,
                       first_train_time, last_train_time,
                       loss_type,
                       flux_thresholds, quantile_levels, ords, leads) {
  x <- dataset %>%
    filter(between(time, first_train_time, last_train_time)) %>%
    pull(flux)
  i <- which(dataset$time == last_train_time)
  y <- dataset$flux[i + leads]
  thresholds_levels <- make_thresholds_levels(
    x, flux_thresholds, quantile_levels
  )
  x <- x - mean(x)
  calc_ar_preds(x, y, loss_type, thresholds_levels, ords, leads)
}

# Fit baseline models -----------------------------------------------------

calc_baseline_preds <- function(x, y, thresholds_levels, leads) {
  foreach (lead = leads, .combine = bind_rows) %do% {
    foreach (i = 1:nrow(thresholds_levels), .combine = bind_rows) %do% {
      tibble(
        ord = 1, lead,
        type = thresholds_levels$type[i],
        flux_threshold = thresholds_levels$flux_threshold[i],
        quantile_level = thresholds_levels$quantile_level[i],
        opt_pred_coefs = list(NULL),
        lin_pred_threshold = NA_real_,
        obs = y[leads == lead] >= flux_threshold,
        pred = last(x) >= flux_threshold
      )
    }
  }
}

fit_baseline_mod <- function(dataset,
                             first_train_time, last_train_time,
                             flux_thresholds, quantile_levels, leads) {
  x <- dataset %>%
    filter(between(time, first_train_time, last_train_time)) %>%
    pull(flux)
  i <- which(dataset$time == last_train_time)
  y <- dataset$flux[i + leads]
  thresholds_levels <- make_thresholds_levels(
    x, flux_thresholds, quantile_levels
  )
  calc_baseline_preds(x, y, thresholds_levels, leads)
}

# Summarize results -------------------------------------------------------

make_confusion_mat <- function(obs_preds) {
  # Use factor() to prevent zero rows and columns from being dropped
  with(
    obs_preds,
    table( # Specify the levels so that they'll be sorted properly
      pred = factor(pred, levels = c("TRUE", "FALSE")),
      obs = factor(obs, levels = c("TRUE", "FALSE"))
    )
  )
}

calc_event_rate <- function(confusion_mat) {
  if (is.null(confusion_mat)) return(NA_real_)
  sum(confusion_mat[, "TRUE"]) / sum(confusion_mat)
}

calc_alarm_rate <- function(confusion_mat) {
  if (is.null(confusion_mat)) return(NA_real_)
  sum(confusion_mat["TRUE", ]) / sum(confusion_mat)
}

# The true positive rate (TPR) is the probability of an alarm given an event
calc_tpr <- function(confusion_mat) {
  if (is.null(confusion_mat)) return(NA_real_)
  # If there are no events, return 1, the best possible TPR
  coalesce(confusion_mat["TRUE", "TRUE"] / sum(confusion_mat[, "TRUE"]), 1)
}

# The false positive rate (FPR) is the probability of an alarm given a non-event
calc_fpr <- function(confusion_mat) {
  if (is.null(confusion_mat)) return(NA_real_)
  # If there are no non-events, return 0, the best possible FPR
  coalesce(confusion_mat["TRUE", "FALSE"] / sum(confusion_mat[, "FALSE"]), 0)
}

# The precision is the probability of an event given an alarm
calc_precision <- function(confusion_mat) {
  if (is.null(confusion_mat)) return(NA_real_)
  confusion_mat["TRUE", "TRUE"] / sum(confusion_mat["TRUE", ])
}

# The True Skill Statistic (TSS) is defined as the difference between the TPR
# and the FPR
calc_tss <- function(confusion_mat) {
  if (is.null(confusion_mat)) return(NA_real_)
  calc_tpr(confusion_mat) - calc_fpr(confusion_mat)
}

# Miscellaneous -----------------------------------------------------------

make_thresholds_levels <- function(x, flux_thresholds, quantile_levels) {
  thresholds_levels1 <- tibble(
    type = "threshold",
    flux_threshold = flux_thresholds,
    quantile_level = map_dbl(
      flux_threshold, function(flux_threshold) mean(x <= flux_threshold)
    )
  )
  thresholds_levels2 <- tibble(
    type = "level",
    quantile_level = quantile_levels,
    flux_threshold = map_dbl(quantile_level, estim_quantile, x = x)
  )
  bind_rows(thresholds_levels1, thresholds_levels2) %>% arrange(flux_threshold)
}