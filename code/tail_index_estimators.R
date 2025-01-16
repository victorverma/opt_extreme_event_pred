fit_gp_mod <- function(x, threshold) {
  tryCatch(
    decluster(x, threshold, method = "intervals") %>%
      fevd(threshold = threshold, type = "GP"),
    error = function(e) NULL
  )
}

calc_clust_maxes <- function(x, threshold) {
  x %>%
    decluster(threshold, method = "intervals") %>%
    attr("clusters") %>%
    tapply(x[x > threshold], ., max) %>%
    as.vector()
}

calc_prob_int_trans <- function(x, threshold, coefs) {
  x %>%
    calc_clust_maxes(threshold) %>%
    `-`(threshold) %>%
    pevd(scale = coefs["scale"], shape = coefs["shape"], type = "GP")
}

do_ks_test <- function(x, threshold, coefs) {
  if (is.null(coefs)) {
    return(tibble(ks_test_stat = NA_real_, ks_p_val = NA_real_))
  }
  x %>%
    calc_prob_int_trans(threshold, coefs) %>%
    ks.test(punif) %>%
    with(tibble(ks_test_stat = statistic, ks_p_val = p.value))
}

calc_extreme_quantile <- function(x, p) {
  if (length(x) < 10) stop("length(x) < 10")
  ord <- order(x, decreasing = TRUE)
  x10 <- x[ord[10]]
  if (quantile(x, 0.75) > x10) return(NA_real_)
  results <- tibble(
    threshold = seq(quantile(x, 0.75), x10, length.out = 10),
    mod = map(threshold, fit_gp_mod, x = x),
    coefs = map(mod, ~ tryCatch(distill(.x), error = function(e) NULL)),
    aic = map_dbl(
      mod, ~ tryCatch(summary(.x)$AIC, error = function(e) NA_real_)
    ),
    ks_test_results = map2(threshold, coefs, do_ks_test, x = x)
  ) %>%
    unnest(ks_test_results)
  if (all(is.na(results$ks_test_stat))) return(NA_real_)
  results <- slice_min(results, ks_test_stat)
  coefs <- results$coefs[[1]]
  threshold <- results$threshold[1]
  p0 <- mean(x <= threshold)
  if (!between((p - p0) / (1 - p0), 0, 1)) return(NA_real_)
  q <- qevd(
    (p - p0) / (1 - p0),
    scale = coefs["scale"], shape = coefs["shape"], threshold = threshold,
    type = "GP"
  )
  threshold + unname(q)
}