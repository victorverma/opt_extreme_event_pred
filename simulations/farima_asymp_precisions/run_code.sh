#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <num_cores>"
  exit 1
fi
num_cores=$1

Rscript --vanilla farima_asymp_precisions.R 400 400 5 $num_cores
Rscript --vanilla make_plot.R
