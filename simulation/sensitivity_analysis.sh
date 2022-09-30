#!/bin/sh

# For different values of the hyperparameter for piecewise constant hazards lambda, please run

Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 2 -s 5994 -i 50000 -b 25000 -c 1 -l 30 -m 5 -n 1 -o 10
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 3 -s 5994 -i 50000 -b 25000 -c 1 -l 40 -m 5 -n 1 -o 10
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 4 -s 5994 -i 50000 -b 25000 -c 1 -l 60 -m 5 -n 1 -o 10
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 5 -s 5994 -i 50000 -b 25000 -c 1 -l 70 -m 5 -n 1 -o 10

# For different values of the hyperparameter for coefficients beta, please run

Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 6 -s 5994 -i 50000 -b 25000 -c 1 -l 50 -m 1 -n 1 -o 10
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 7 -s 5994 -i 50000 -b 25000 -c 1 -l 50 -m 3 -n 1 -o 10
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 8 -s 5994 -i 50000 -b 25000 -c 1 -l 50 -m 7 -n 1 -o 10
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 9 -s 5994 -i 50000 -b 25000 -c 1 -l 50 -m 9 -n 1 -o 10

# For different values of the hyperparameter for Dirichlet parameters xi, please run

Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 10 -s 5994 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 0.2 -o 10
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 11 -s 5994 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 0.5 -o 10
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 12 -s 5994 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 2 -o 10
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 13 -s 5994 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 5 -o 10

# For testing the selection of knots, please run

Rscript ../src_MCMC/MCMC_source.R -p simknots -v 1 -K 3 -r 1 -s 4146 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 1 -o 10 -M 3
Rscript ../src_MCMC/MCMC_source.R -p simknots -v 1 -K 3 -r 2 -s 2941 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 1 -o 10 -M 4
Rscript ../src_MCMC/MCMC_source.R -p simknots -v 1 -K 3 -r 3 -s 2981 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 1 -o 10 -M 5
Rscript ../src_MCMC/MCMC_source.R -p simknots -v 1 -K 3 -r 4 -s 2882 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 1 -o 10 -M 6
Rscript ../src_MCMC/MCMC_source.R -p simknots -v 1 -K 3 -r 5 -s 7161 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 1 -o 10 -M 7
Rscript ../src_MCMC/MCMC_source.R -p simknots -v 1 -K 3 -r 6 -s 3364 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 1 -o 10 -M 8