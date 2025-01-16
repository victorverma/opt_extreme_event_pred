#!/bin/bash

dataset_file='../make_dataset/dataset.RData'
num_cores=8
if [ ! -z "$1" ]; then
  num_cores=$1
fi
max_num_top_ord_stats=100
alpha_hat=NULL
d_hat=NULL
num_coefs=1e5
flux_thresholds=""
quantile_levels="0.9 0.95 0.99"
ords="1 168"
leads="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24"
plot_obj_fun=TRUE

cd ../make_hill_plots/
Rscript --vanilla ../../make_hill_plots.R \
  $dataset_file $num_cores $max_num_top_ord_stats $alpha_hat &
cd ../fit_baseline_mods/
Rscript --vanilla ../../fit_baseline_mods.R \
  $dataset_file \
  $num_cores \
  "$flux_thresholds" "$quantile_levels" "$leads" &
cd ../fit_farima_mods/
Rscript --vanilla ../../fit_farima_mods.R \
  $dataset_file \
  $num_cores \
  $alpha_hat $d_hat \
  $num_coefs \
  "$flux_thresholds" "$quantile_levels" "$ords" "$leads" \
  $plot_obj_fun &
cd ../fit_ols_ar_mods/
Rscript --vanilla ../../fit_ar_mods.R \
  $dataset_file \
  $num_cores \
  ols \
  "$flux_thresholds" "$quantile_levels" "$ords" "$leads" &
cd ../fit_lad_ar_mods/
Rscript --vanilla ../../fit_ar_mods.R \
  $dataset_file \
  $num_cores \
  lad \
  "$flux_thresholds" "$quantile_levels" "$ords" "$leads" &
wait
cd ../review_results/
quarto render review_results.qmd --to html
cd ..