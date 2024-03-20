# hdmax2

<img src="https://raw.githubusercontent.com/bcm-uga/hdmax2/package/hdmax2_hex.png" width="130" align="right">

R package `hdmax2` performs high dimension mediation analysis. This package allows users to estimate the indirect effects of mediators, calculate the overall indirect effect of mediators, and facilitates the execution of high-dimensional mediation analysis. 

The HDMAX2 method includes unobserved confounding factors through a latent factor mixed model approach. It supports the use of exposures of various types and the consideration of both continuous and binary outcomes.


## Installation 

```
conda create --name hdmax2  

conda activate hdmax2
conda install -c conda-forge r-base
conda install conda-forge::r-devtools
conda install conda-forge::r-rcppeigen
conda install conda-forge::r-lme4

R
> install.packages("ggplot2")
> install.packages("prettydoc")
> install.packages("fdrtool")
> install.packages("mediation")


> devtools::install_github("bcm-uga/hdmax2")

```

## Usage
Load our simulated data and choose K from pca analysis (for example) on potential mediators.
```
# load data
simu_data = hdmax2::simu_data

# exposure
X = simu_data$X_continuous

# outcome
Y = simu_data$Y_continuous

# potential mediators
M = simu_data$M

# Observed confounding factors
covar = simu_data$age

# latent factor number estimated with pca 
K = 5
```

### Run step 1

Step 1 involves estimating latent factors using the `LFMM` algorithm. Subsequently, association studies are conducted between exposure and potential mediators, and between potential mediators and the outcome. These studies include both latent factors and observed confounding factors (when applicable). The mediation test using the *max-squared* test is performed using p-values obtained from the previous association studies.

```
hdmax2_step1 = hdmax2::run_AS(X = X,
                              Y = Y,
                              M = M,
                              K = K)
```

### Select mediators

Mediators are selected using *max-squared* test p-values. Here's we choose to select top 10 p-values as dection method.

```
mediators_top10 = M[,names(sort(hdmax2_step1$max2_pvalues)[1:10])]

```
### Run step 2

With selected mediators step 2 proceed to different effect estimations.

```
hdmax2_step2 = hdmax2::estimate_effect(object =hdmax2_step1,
                                    m = mediators_top10)
```

### Plot results

Those effects can be analyze with proposed plots.

```
library(ggplot2)
hdmax2::plot_hdmax2(hdmax2_step2, N_med = 10)
```

## Bug report / Help

Please open an issue if you find a bug.

## Jumentier et al. 2023 release access

Previous code version associated with cited publication is always available. To access this code version choose v1.0.0.0 in tags tab.

To install this version:
```
devtools::install_github("bcm-uga/hdmax2", ref="v1.0.0.0")
```

## References 

- Jumentier B, CC Barrot, M Estavoyer, J Tost, B Heude, O Francois, and J Lepeule (2023). High-Dimensional Medi- 220
ation Analysis: A New Method Applied to Maternal Smoking, Placental DNA Methylation, and Birth Out- 221
comes. Environmental Health Perspectives 131 (). Publisher: Environmental Health Perspectives, 047011. 222
https://doi.org/10.1289/EHP11559.

