sim_settings <- tibble(
  dataset_type = "AR",
  ar_ord = 5,
  ar_roots = list(c(-3 + 1i, -3 - 1i, 4 + 2i, 4 - 2i, 2)),
  min_ar_root_abs_val = 2,
  ar_coefs = list(make_ar_coefs(ar_roots[[1]])),
  farima_ord = NA_real_
) %>%
  expand_grid(
    sampler = list(rt),
    alpha = 1,
    lead = 1,
    train_set_size = 1e4,
    test_set_size = 1e6,
    ps = list(c(0.9, 0.95, 0.99, 0.999)),
    mod_types = list(
      list(
        list(
          mod_type = "oracle", is_oracle = TRUE, quantile_type = "empirical"
        ),
        list(mod_type = "lad", is_oracle = FALSE, quantile_type = "empirical"),
        list(mod_type = "lad", is_oracle = FALSE, quantile_type = "extreme")
      )
    ),
    run_num = 1:100
  )
