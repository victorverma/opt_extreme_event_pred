here::i_am("simulations/quantile_type_impact/make_plot.R")

library(here)
library(tidyverse)

load(here("simulations/quantile_type_impact/results.RData"))

plot <- results %>%
  mutate(
    mod_type2 = map_chr(mod_type, ~ str_c(.x[-2], collapse = "/")),
    mod_type2 = case_when(
      mod_type2 == "lad/empirical" ~ "Non-Oracle/Empirical",
      mod_type2 == "lad/extreme" ~ "Non-Oracle/Extreme",
      mod_type2 == "oracle/empirical" ~ "Oracle"
    ),
    p = as_factor(p)
  ) %>%
  ggplot(aes(p, test_precision, color = mod_type2)) +
  geom_boxplot() +
  geom_hline(
    aes(yintercept = oracle_asymp_precision),
    color = "orange", linetype = "dashed"
  ) +
  scale_color_manual(values = c("blue", "black", "orange")) +
  ylim(c(0, 1)) +
  labs(y = "Test Set Precision", color = "Predictor") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
ggsave(
  filename = here("figures/quantile_type_impact.pdf"),
  plot,
  width = 6, height = 3
)