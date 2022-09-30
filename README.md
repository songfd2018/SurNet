# SurNet: Survival Mixed Membership Blockmodel

We propose a Bayesian hierarchical model---Survival Mixed Membership Blockmodel (SMMB) to predict and infer the response behaviors for any given pair of actors. You can reproduce our results step by step as follows.

## Compile C++ source code

First, please compile C++ source code for the MCMC algorithm by running

```
cd src_MCMC
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

It will generate a 150-actor network with time-to-event data between pairs of users. The observed data, including survival times, censoring indicators and covariates, will be stored as the `Dat` variable in the workspace `ObsData/Obs_SMMSB_v1.RData`, while the whole workspace with the true values of parameters will be saved in the workspace `ObsData/workspace_SMMSB_v1.RData` to verify the inference results. The adjacent matrix of the synthetic network is drawn as `Images/Adjacent_matrix_grey_v1.jpg` as Figure 1(b) in the main text.

3. To reproduce the statistical inference in our manuscript, please run

```
Rscript ../src_MCMC/MCMC_source.R -p simulation -v 1 -K 3 -r 1 -s 5994 -i 50000 -b 25000 -c 1 -l 50 -m 5 -n 1 -o 10
```

where 

* -p: the project name, 
* -v: the version,
* -K: the number of roles,
* -r: the replication,
* -s: the seed for MCMC algorithm,
* -i: the iteration number of MCMC algorithm,
* -b: the number of iterations as burnins,
* -c: the number of cores for parallel computing
* -l: the hyperparameter kappa for the piecewise constant hazards,
* -m: the hyperparameter sigma^2 for the coefficients of covariates,
* -n: the hyperparameter a for the Dirichlet parameters,
* -o: the number of extra iterations to sample for calculating the DIC.

The posterior inference of parameters can be found in the folder Inference.  The file is named `Inference_simulation_SMMB_K3_v1_r1.RData` to present its corresponding project, its number of roles, its version number and its replication number. 

4. To select the number of roles, we vary the number of roles from 1 to 5 to compare their corresponding DIC values by running

```
./select_role_number.sh
```

We suggest to run the command for each role number in parallel. 

5. Then, please reproduce the figures and tables by running

```
R --vanilla --slave < simulation_analysis.R
```

First, DIC can correctly select the optimal role number as the true value `K = 3` in the manuscript. Then, we compare the estimated parameter values with the true values. To compare the estimation of the role pair for each actor pair, we generate the contingency table between the estimated role pairs and their true values. Next, we compare the estimated baseline survival probabilities at each knot with the true values to show that our proposed model can correctly recover the piecewise constant hazards and the intercept terms. Besides, we list the true values, the posterior mean and variance, and the 95% credible intervals of all coefficients.

### Sensitivity analysis

We conduct sensitivity analysis to check whether our inference is robust to the selection of hyparameters or knots. 

1. To conduct the posterior inference for each case, please run

```
./sensitivity_analysis.sh
```

2. Then, we calculate the ARIs between the estimated role pairs by the original setting and by the alternative. Please run

```
R --vanilla --slave < sensitivity_analysis.R
```

It will list the ARIs under different settings shown in Tables 3 and 4.

### Model comparison

We would like to compare the performance of the SMMB with that of two alternative strategies. The first strategy is to fit a unique semiparametric cure rate model (SCRM) for all the survival times, called an overall SCRM. If the number of roles `K` is set as one, then the SMMB degenerate to an overall SCRM. Thus, we can directly infer the piesewise constant hazards and the coefficients for the covariates by setting `K=1` in our proposed MCMC algorithm. The second strategy is to fit a separate SCRM for the survival times between each pair of users, which is termed as pairwise SCRMs. To infer the pair-specific hazards and coefficients for each pair, we apply an MCMC algorithm for each pair separately. Finally, to calibrate the performance of different models, we take the L measure of the posterior predictive distribution as the criterion.

1. To obtain the posterior sampling of pairwise SCRMs, we run the MCMC algorithm by

```shell
R --vanilla --slave < MCMC_pairwiseSCRM.R
```

The posterior sampling will be stored in the folder Inference as `Inference_simulation_pairwiseSCRM_pairX.RData` for `X`-th user pair.

2. To generate data replicate from the posterior predictive distribution and calculate the L measure, we run

```shell
R --vanilla --slave < PPC_Lmeasure.R
```

We collect the L measure of three models as "Lmeasure_Collect" varaible, which is saved in `PPC_results/Lmeasure_PPC.RData`.

## Simulation with different heterogeneity and sparsity

In section 5.2, we evaluate the performance of our proposed MCMC algorithm with different heterogeneity levels, different numbers of active edges and different numbers of observed survival times per active edge. We consider 12 different simulation settings and generate 100 replicated dataset for each setting. Here, we take the setting with `K=3` roles, the connectivility probability `p_c=0.3` and the average number of response time per edge `n_ij` follows `Pois(25)+5` as an example to introduce the workflow.

1.  Please first enter the simulation directory
   
   ```shell
   cd more_simulation
   ```

2. Generate 100 replicated dataset with the same parameter values in the `K3p20r25` folder by running
   
   ```shell
   R --vanilla --slave < Dataset_generation_K3p30r25.R
   ```

3. Enter each folder and conduct the MCMC algorithm by running
   
   ```shell
   cd K3p20r25/sim_v1
   Rscript ../../..//src_MCMC/MCMC_revised.R -p sim -v 1 -K 3 -r 1 -s 7158 -i 50000 -b 25000 -t 10 -l 50 -m 5 -n 1 -o 10 -w 10 (-b 1)
   ```
   
   Compared with the previous version `MCMC_source.R`, `MCMC_revised.R` incorporate a global swapping Metropolis step to speed up convergence. Besides the arguments of `MCMC_source.R`, the argument `-w` controls the number of iteration to conduct a global swapping metropolis step (by default, 10).

### Output of the MCMC algorithm

The posterior inference of parameters can be found in the folder "Inference".  The file is named as `Inference_sim_SMMB_K2_v1_r1.RData` to show its corresponding project `sim`, its number of roles `K2`,  its version number `v1` and its replication number `r1`. 

#### Posterior samples, mean, sd, and credible intervals

The workspace stores the posterior samples, mean, sd and 95% credible intervals of parameters. For example, "lambda_record", "lambda.est", "lambda.sd", and "lambda.ci" give the posterior samples, mean, standard deviation and 95% credible interval of lambda. Based on the workspace, we can compute the biases, SD, RMSE and Coverage probabilities across 100 replicates. In general, we run multiple chains from different seeds to check the convergence of the MCMC algorithm.

#### Likelihood and information criteria

The function will return BIC, WAIC and log-likelihood values by variables "BIC", "WAIC" and "log_like". Only when the arugment -d is 1, then the conditional DIC is computed and stored in the variable "DIC".

#### Running time

"end_time" and "start_time" will store the starting and end time of the MCMC algorithm, so the running time can be obtain by running the following command in R

```r
difftime(end_time, start_time, units = "mins")
```

### Reproduce Figures 1 and 2

In Figure 1, we draw the coverage probablility of 95% credible intervals for each parameter in 100 replicated datasets. We provide the cached workspace `RecoveryResults.Rdata` containing the converage probabilities of parameters for 12 different settings, so Figure 1 can be reproduced by running

```shell
R --vanilla --slave < RecoveryPerformance.R
```

Similarly, we store the running time of the first replication under different settings in `RunningTime.Rdata`. To reproduce Figure 2, please run

```shell
R --vanilla --slave < Scalability.R
```

## Enron email corpus

The Enron corpus is the largest publicly available email dataset to date and was released by the Federal Energy Regulatory Commission during its investigation of Enron's bankruptcy. It contains the user information and the time stamps of each email. Following [the enrondata GitHub repository](https://github.com/enrondata/enrondata/blob/master/data/misc/edo\_enron-custodians.txt), we focus on the email folders of 148 Enron users whose positions in the company are available. You can download the pulicly available email corpus from [the website](https://www.cs.cmu.edu/~./enron/). You can reproduce the figures in the manuscript by the following steps.

1. (Optional) To extract the response or censoring time of email interactions among 148 users, please conduct the preprocessing by running

```shell
R --vanilla --slave < EnronPreprocessing.R
```

As a result, the observed data are stored as `Data/Reponse_time_collection.RData` for the analysis without confidential information, including censoring indicators `nu`, survival outcomes `y` and covariate matrix `X`. At the same time, the data frame for the analysis with confidential information and the data frame considering sending behaviors are stored as `Data/Respons_confidential.RData` and `Data/Response_time.RData`, respectively. This step is time-consuming, so we also attach these workspaces in the supplementary materials.

Morevoer, we also generate the binary adjacent matrix stored in `adj_enron_communications.txt` for applying MMSB.

2. Run the MCMC algorithm for the above three circumstances by running 

```shell
./EnronMCMC.sh
```

3. To determine the number of roles, we also vary the number of roles from 1 to 5 to compare their corresponding DIC values by running

```shell
./EnronSelection.sh
```

The posterior inference of parameters can be found in the `Inference` folder. `v1` stands for analyzing without confidential information, `v2` represents the analysis with confidiential information, and `v3` indicates to consider sending behaviors. 

4. Finally, to reproduce all figures and tables in the real data analysis, including the scatter plot and boxplot of role proportions, the baseline survival functions, and output the coefficient effects, please run

```shell
R --vanilla --slave < EnronAnalysis.R
```

For the analysis without confidiential information, the R code draws the scatter plots of the actor-specific role proportions and the line charts of the baseline survival functions for all role pairs in `Figure 3`. It then gives the estimated role proportions of four CEOs in `Table 2` by Mixed Membership Stochastics Blockmodel (MMSB) and SMMB, respectively. Finally, it outputs the posterior mean, the posterior SD and the 95% credible interval of coefficients in `Table 3`.

Meanwhile, the R code will reproduce the corresponding scatter plots, line charts and coefficient estimation for another two circumstances, that is, analysis with confidential information and considering sending behaviors.