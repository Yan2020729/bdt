### bdt
The R package **bdt** implements adapted bootstrap diagnostic tool (BDT), originally described in Bahamyirou A, Blais L, Forget A, Schnitzer ME. *Understanding and diagnosing the potential for bias when using machine learning methods with doubly robust causal estimators*. Statistical Methods in Medical Research. 2019 Jun;28(6):1637-50. Data-adaptive methods can produce a separation of the two exposure groups in terms of propensity score densities which can lead to biased finite-sample estimates of the treatment effect. Bahamyirou et al. presented a diagnostic tool to identify the instability in the estimation of marginal causal effects in the presence of near practical positivity violations.  To expedite the application of this diagnostic tool, we develop **bdt** package which is based on a bootstrap resampling of the subjects and simulation of the outcome data conditional on observed treatment and covariates.  You can install the current version of with: 

```{r, eval=FALSE}
install.packages("devtools")
devtools::install_github("Yan2020729/bdt")
```

This repository hosts the (commented) source code. For more theoretical details, please refer to the original paper: https://doi.org/10.1177%2F0962280218772065.

## stages
- Stage 1: Propensity scores estimation. For a given observed dataset $(Y, A, W)$ with $n$ i.i.d. observations where the outcome $Y$ is either a continuous or a binary variable, the treatment $A$ is a binary variable and $W$ could be a vector, matrix or dataframe containging baseline covariates, the propensity scores are estimated by logistic regression and/or SL methods. 
- Stage 2: Fitting outcome models. For the two treatment subgroups of subjects with $A=1$ and $A=0$, fit a linear regression of the outcome $Y$ on baseline covariates $W$ respectively where the model only contains main terms of $W$. We denote the estimated coefficients and their corresponding variance as $\hat{\beta}_{0,a}$, $\hat{\beta}_{W,a}$ and $\hat{\sigma}^2_a$ for the realization of treatment $a\in\{0,1\}$.
- Stage 3: Computing the "true" effect for the bootstrap data generating distribution $\hat{P}_0$. It is obtained by computing the difference of the two potential outcomes for all subjects using the estimated coefficients from 2nd step.
- Stage 4: Simulation of the bootstrap datasets. Sample the $n$ subjects with replacement then replace their observed outcomes with new outcome values which are estimated based on the coefficients obtained in the 2nd step. Specifically, if the outcome $Y$ is continuous,  the two potential outcomes $Y^1$ and $Y^0$ are generated from a $\mathcal{N}(\mu_a, \sigma^2_a)$ distribution with mean $\mu_a=\hat{\beta}_{0,a}+W\hat{\beta}_{W,a}, a\in\{0,1\}$. Similarly, if $Y$ is binary, $Y^1$ and $Y^0$ are generated from a $Binomial(n_a, p_a)$ distribution where $n_a$ represents the number of subjects with $A=a$ in the resampled data and $logit(p_a)=\hat{\beta}_{0,a}+W\hat{\beta}_{W,a}, a\in\{0,1\}$. 
- Stage 5: Estimation of target parameters. For the resampled data with simulated outcomes, three candidate estimators TMLE, A-IPTW and IPTW are applied to estimate the targeted parameter using the``true'' specification of the outcome model (same with the linear regression model in step 2) and the propensity score is modeled by logistic regression and/or SL. User could define specific bounds for the values of $g_n$ where the default value is 0.025. 
- Stage 6: Repeat steps (3) and (5) $M$ times. Then bias could be calculated by comparing the mean of the estimators $\hat{\psi}$ across $M$ bootstrap samples with the true value of the target parameter under the bootstrap data generating distribution. Also, the coverage rates of three estimators are obtained based on the $95\%$ confidence intervals.



