# SurNet: Survival Mixed Membership Blockmodel

We propose a Bayesian hierarchical model---Survival Mixed Membership Blockmodel (SMMB) to predict and infer the response behaviors for any given pair of actors. You can reproduce our results step by step as follows.

## Compile C++ source code 

First, please compile C++ source code for the MCMC algorithm by running

```
cd src_MCMC
chmod 755 compile.sh
./compile.sh
```

## Simulation study 

1. Please first enter the simulation directory

```
cd ../simulation/
```

2. To generate synthetic data, please run 

```
R --vanilla --slave < simulate_SMMB.R
```

It will generate a 150-actor network with time-to-event data between pairs of users. The observed data, including survival times, censoring indicators and covariates, will be stored as the "Dat" variable in the workspace "ObsData/Obs_SMMSB_v1.RData", while the whole workspace with the true values of parameters will be saved in the workspace `ObsData/workspace_SMMSB_v1.RData` to verify the inference results. The adjacent matrix of the synthetic network is drawn as `Images/Adjacent_matrix_grey_v1.jpg` as Figure 1(b) in the main text.

3. To reproduce the statistical inference in our manuscript, please run

```
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 1 -s 5994 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 1 -o 10
```
where 

*-p: the project name, 
*-v: the version,
*-K: the number of roles,
*-r: the replication,
*-s: the seed for MCMC algorithm,
*-i: the iteration number of MCMC algorithm,
*-b: the number of iterations as burnins,
*-c: the number of cores for parallel computing
*-l: the hyperparameter kappa for the piecewise constant hazards,
*-m: the hyperparameter sigma^2 for the coefficients of covariates,
*-n: the hyperparameter a for the Dirichlet parameters,
*-o: the number of extra iterations to sample for calculating the DIC.

The posterior inference of parameters can be found in the folder Inference.  The file is named `Inference_simulation_SMMB_K3_v1_r1.RData` to present its corresponding project, its number of roles, its version number and its replication number. 

4. To select the number of roles, we vary the number of roles from 1 to 5 to compare their corresponding DIC values by running

```
chmod 755 select_role_number.sh
./select_role_number.sh
```

We suggest to run the command for each role number in parallel. 

5. Then, please reproduce the figures and tables by running

```
R --vanilla --slave < simulation_analysis.R
```

First, DIC can correctly select the optimal role number as the true value K = 3 in the manuscript. Then, we compare the estimated parameter values with the true values. To compare the estimation of the role pair for each actor pair, we generate the contingency table between the estimated role pairs and their true values as Table 2. Next, we compare the estimated baseline survival probabilities at each knot with the true values to show that our proposed model can correctly recover the piecewise constant hazards and the intercept terms as Figure 2. Besides, we list the true values, the posterior mean and variance, and the 95% credible intervals of all coefficients in Table 1.


## Sensitivity analysis 

We conduct sensitivity analysis to check whether our inference is robust to the selection of hyparameters or knots. 
1. To conduct the posterior inference for each case, please run

```
chmod 755 sensitivity_analysis.sh
./sensitivity_analysis.sh
```

2. Then, we calculate the ARIs between the estimated role pairs by the original setting and by the alternative. Please run

```
R --vanilla --slave < sensitivity_analysis.R
```

It will list the ARIs under different settings shown in Tables 3 and 4.

## Model comparison 

We would like to compare the performance of the SMMB with that of two alternative strategies. The first strategy is to fit a unique semiparametric cure rate model (SCRM) for all the survival times, called an overall SCRM. If the number of roles K is set as one, then the SMMB degenerate to an overall SCRM. Thus, we can directly infer the piesewise constant hazards and the coefficients for the covariates by setting K=1 in our proposed MCMC algorithm. The second strategy is to fit a separate SCRM for the survival times between each pair of users, which is termed as pairwise SCRMs. To infer the pair-specific hazards and coefficients for each pair, we apply an MCMC algorithm for each pair separately. Finally, to calibrate the performance of different models, we take the L measure of the posterior predictive distribution as the criterion.

1. To obtain the posterior sampling of pairwise SCRMs, we run the MCMC algorithm by

```
R --vanilla --slave < MCMC_pairwiseSCRM.R
```

The posterior sampling will be stored in the folder Inference as `Inference_simulation_pairwiseSCRM_pairX.RData` for `X`-th user pair.

2. To generate data replicate from the posterior predictive distribution and calculate the L measure, we run

```
R --vanilla --slave < PPC_Lmeasure.R
```

We collect the L measure of three models as "Lmeasure_Collect" varaible, which is saved in `PPC_results/Lmeasure_PPC.RData`.


## Enron email corpus 


The Enron corpus is the largest publicly available email dataset to date and was released by the Federal Energy Regulatory Commission during its investigation of Enron's bankruptcy. It contains the user information and the timestamps of each email. Following [the enrondata GitHub repository](https://github.com/enrondata/enrondata/blob/master/data/misc/edo\_enron-custodians.txt), we focus on the email folders of 148 Enron users whose positions in the company are available. 

1. Please first enter the "enron" directory

```
cd ../enron/
```

2. To extract the response or censoring time of email interactions among 148 users, please conduct the preprocessing by running

```
R --vanilla --slave < preprocessing.R
```

As a result, the observed data are stored as a matrix called `sur_collect` in the workspace `Reponse_time_collection.RData`, including censoring indicators nu, survival outcomes y and covariate matrix X. At the same time, we also generate the binary adjacent matrix stored in `adj_enron_communications.txt` for applying MMSB. 

3. To summarize the Enron dataset, please run

```
R --vanilla --slave < summarize_enron.R
```

We first provide the summary statistics shown in Table 5. Then, we draw the heatmap of the adjacent matrix as Figure 1(a).

4. To conduct statistical inference, please run

```
Rscript ../src_MCMC/MCMC_source.R -p enron -v 1 -K 2 -r 1 -s 9822 -k 5 -i 100000 -b 50000 -c 1 -l 50 -m 5 -n 1 -o 10
```

where 
	
*-k: the number of knots.

Here, we run totally 100,000 iterations and regard the first 50,000 as burnins in the MCMC algorithm.

5. To select the number of roles, we also vary the number of roles from 1 to 5 to compare their corresponding DIC values by running

```
chmod 755 select_role_number.sh
./select_role_number.sh
```

6. To reproduce all figures and tables in the real data analysis, please run

```
R --vanilla --slave < real_data_analysis.R
```

The code draws the scatter plots of the actor-specific role proportions as Figure 3 and the line charts of the baseline survival functions for all role pairs as Figure 4. It then gives the estimated role proportions of four CEOs in Table 3 by Mixed Membership Stochastics Blockmodel (MMSB) and SMMB, respectively. Finally, it outputs the posterior mean, the posterior SD and the 95% credible interval of coefficients as Table 4. 
