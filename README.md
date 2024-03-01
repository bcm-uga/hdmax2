# hdmax2

<img src="https://raw.githubusercontent.com/bcm-uga/hdmax2/package/hdmax2_hex.png" width="130" align="right">

R package {hdmax2} performs high dimension mediation analysis. This package allows users to estimate the indirect effects of mediators, calculate the overall indirect effect of mediators, and facilitates the execution of high-dimensional mediation analysis. 

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
> install.packages("prettydoc")
> install.packages("FactoMineR")
> install.packages("factoextra")
> install.packages("fdrtool")
> install.packages("mediation")
> devtools::install_github("bcm-uga/LEA")

> devtools::install_github("bcm-uga/hdmax2")

```

## Usage

```
library(hdmax2)
load("data/simu_data.RData")

X_binary = as.matrix(simu_data$X_binary)
Y_continuous = as.matrix(simu_data$Y_continuous)
M = as.matrix(simu_data$M)

res_step1 = hdmax2::run_AS(X_matrix = X_binary ,
                 Y_matrix = Y_continuous,
                 M_matrix = M, 
                 K = K,
                 X_type = "binary",
                 Y_type = "continuous",
                 M_type = "methylation",
                 multivariate = FALSE,
                 covar = cbind(age, gender),
                 diagnostic.plot = F)


mediators_top10 = M[,names(sort(res_step1$max2)[1:10])]

res_step2 = hdmax2::estimate_effect(object = res_step1,
                                    m = mediators_top10,
                                     boots = 100,
                                    sims = 100)

hdmax2::plot(res_step2, N_med = 10)

```

## Bug report / Help

Please open an issue if you find a bug.

## References 

- Jumentier B, CC Barrot, M Estavoyer, J Tost, B Heude, çO Fran, and J Lepeule (2023). High-Dimensional Medi- 220
ation Analysis: A New Method Applied to Maternal Smoking, Placental DNA Methylation, and Birth Out- 221
comes. Environmental Health Perspectives 131 (). Publisher: Environmental Health Perspectives, 047011. 222
https://doi.org/10.1289/EHP11559.

