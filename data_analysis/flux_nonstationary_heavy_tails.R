here::i_am("data_analysis/flux_nonstationary_heavy_tails.R")

library(grid)
library(gridExtra)
library(here)
library(tidyverse)

load("~/research/flux_forecasting/data/processed/flux/1m_tbls.RData")

panel1 <- flux_tbl %>%
  mutate(row_num = row_number()) %>%
  filter(between(year(time), 2010, 2022), row_num %% 10 == 0) %>%
  ggplot(aes(time, flux)) +
  geom_point(shape = ".") +
  scale_y_log10() +
  labs(x = NULL, y = NULL) +
  theme_classic()
panel2 <- flux_tbl %>%
  filter(
    between(
      time,
      as_datetime("2011-09-01 00:00:00"),
      as_datetime("2014-05-31 23:59:59")
    )
  ) %>%
  ggplot(aes(time, flux)) +
  geom_line() +
  labs(x = NULL, y = NULL) +
  theme_classic()
pdf(
  here("oberwolfach/figures/flux_nonstationary_heavy_tails.pdf"),
  width = 6, height = 3
)
grid.arrange(
  panel1, panel2,
  nrow = 1,
  bottom = "Time",
  left = textGrob(expression("X-Ray Flux" ~~ (W / m^2)), rot = 90)
)
dev.off()