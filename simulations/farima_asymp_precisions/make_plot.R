here::i_am("simulations/farima_asymp_precisions/make_plot.R")

library(foreach)
library(grid)
library(gridExtra)
library(here)
library(tidyverse)

source(here("simulations/utils.R"))

load(here("simulations/farima_asymp_precisions/results.RData"))

plot <- results %>%
  filter(lead == 1) %>%
  ggplot(aes(alpha, d, z = oracle_asymp_precision)) +
  geom_contour_filled() +
  geom_point(x = 1.4, y = 0.19, color = "red") +
  scale_fill_brewer(palette = "Blues") +
  labs(x = expression(alpha), fill = expression(lambda[h]**(opt))) +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    legend.position = "inside",
    legend.position.inside = c(0.12, 0.7),
    legend.background = element_rect(fill = NA)
  )

settings <- expand_grid(
  alpha = 1.4,
  d = 0.19,
  leads = list(1:24)
)
farima_coefs <- settings %>%
  distinct(d) %>%
  mutate(
    farima_coefs = lapply(
      d, make_farima_coefs, max_lag = 1e6
    ),
    last_to_first_coef_ratio = map_dbl(farima_coefs, ~ last(.x) / first(.x)),
    are_coefs_strictly_decreasing = map_lgl(farima_coefs, ~ all(diff(.x) < 0))
  )
results2 <- foreach (i = 1:nrow(settings)) %do% {
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
plot2 <- results2 %>%
  ggplot(aes(lead, oracle_asymp_precision)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 20, by = 5)) +
  labs(x = "h", y = expression(lambda[h]**(opt))) +
  theme_bw() +
  theme(text = element_text(size = 20))

plot3 <- arrangeGrob(plot, plot2, nrow = 1, widths = c(2, 1))

ggsave(
  filename = here("figures/farima_asymp_precisions.pdf"),
  plot3,
  width = 8.75, height = 5, units = "in", dpi = 600
)