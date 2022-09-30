#!/bin/sh

Rscript ../src_MCMC/MCMC_source.R -p enron -v 1 -K 1 -r 1 -s 247 -M 5 -i 100000 -b 50000 -c 1 -l 50 -m 5 -n 1 -o 10
Rscript ../src_MCMC/MCMC_source.R -p enron -v 1 -K 3 -r 1 -s 7993 -M 5 -i 100000 -b 50000 -c 1 -l 50 -m 5 -n 1 -o 10
Rscript ../src_MCMC/MCMC_source.R -p enron -v 1 -K 4 -r 1 -s 5640 -M 5 -i 100000 -b 50000 -c 1 -l 50 -m 5 -n 1 -o 10
Rscript ../src_MCMC/MCMC_source.R -p enron -v 1 -K 5 -r 1 -s 8952 -M 5 -i 100000 -b 50000 -c 1 -l 50 -m 5 -n 1 -o 10