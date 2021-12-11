#!/bin/sh

Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 1 -r 1 -s 4423 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 1 -o 10
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 2 -r 1 -s 6212 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 1 -o 10
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 4 -r 1 -s 6100 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 1 -o 10
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 5 -r 1 -s 5589 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 1 -o 10