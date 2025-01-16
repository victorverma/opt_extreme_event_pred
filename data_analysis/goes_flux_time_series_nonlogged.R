here::i_am("data_analysis/goes_flux_time_series_nonlogged.R")

library(grid)
library(gridExtra)
library(here)
library(tidyverse)

load("~/research/flux_forecasting/data/processed/flux/1m_tbls.RData")

panel1 <- flux_tbl %>%
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
  theme_bw()
panel2 <- flux_tbl %>%
  filter(
    between(
      time,
      as_datetime("2012-03-06 23:30:00"),
      as_datetime("2012-03-07 01:00:00")
    )
  ) %>%
  ggplot(aes(time, flux)) +
  geom_line() +
  labs(x = NULL, y = NULL) +
  theme_bw()
pdf(here("figures/goes_flux_time_series_nonlogged.pdf"), width = 6, height = 3)
grid.arrange(
  panel1, panel2,
  nrow = 1,
  bottom = "Time",
  left = textGrob(expression("X-Ray Flux" ~~ (W / m^2)), rot = 90)
)
dev.off()

# read_csv("~/Dropbox (University of Michigan)/FlareData/GOES_dataset.csv") %>%
#   filter(between(event_date, as_date("2011-09-01"), as_date("2014-05-31"))) %>%
#   rename(peak_intensity_str = class) %>%
#   mutate(
#     class = case_when(
#       str_sub(peak_intensity_str, 1, 1) == "A" ~ 1e-8,
#       str_sub(peak_intensity_str, 1, 1) == "B" ~ 1e-7,
#       str_sub(peak_intensity_str, 1, 1) == "C" ~ 1e-6,
#       str_sub(peak_intensity_str, 1, 1) == "M" ~ 1e-5,
#       str_sub(peak_intensity_str, 1, 1) == "X" ~ 1e-4,
#     ),
#     multiplier = as.double(map_chr(peak_intensity_str, str_sub, start = 2)),
#     peak_intensity = multiplier * class
#   ) %>%
#   select(!c(class, multiplier)) %>%
#   slice_max(peak_intensity) %>%
#   print(width = Inf)