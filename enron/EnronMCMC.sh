#!/bin/sh

Rscript src/MCMC_source.R -p enron -v 1 -K 2 -r 1 -s 9822 - M 5 -i 100000 -b 50000 -c 1 -l 50 -m 5 -n 1 -o 10

Rscript src/MCMC_revised.R -p enron -v 2 -K 2 -r 1 -s 3508 -M 5 -i 500000 -b 250000 -c 1 -l 50 -m 1 -n 1 -o 10 -w 10

Rscript src/MCMC_revised.R -p enron -v 3 -K 2 -r 1 -s 4228 -M 5 -i 100000 -b 50000 -c 1 -l 50 -m 5 -n 1 -o 10 -d 1 -w 10

