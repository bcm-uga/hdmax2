# hdmax2

<img src="https://raw.githubusercontent.com/bcm-uga/hdmax2/package/hdmax2_hex.png" width="130" align="right">

The `R` package `hdmax2` offers powerful tools for conducting high-dimensional mediation analysis. This method investigates the causal pathways linking exposure variables to outcome variables through intermediary factors known as mediators. These mediators often stem from biological metrics like transcriptomes or methylomes, which present data in high dimensions.

The package is capable of detecting individual mediators, accurately estimating the indirect effects of exposure variables associated with each mediator, and determining an overall indirect effect encompassing all detected mediators. Utilizing a latent factor mixed model methodology, the method effectively mitigates unobserved confounding factors. It accommodates exposures of diverse types and allows for the analysis of both continuous and binary outcomes.




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

# Installing R packages from GitHub may require that users remove older versions and restart their R session
# The package might have been installed in your computer (even though it does not work). Remove it using remove.packages()
#  rs.restartR() if in RStudio
```

## Usage


```
library(hdmax2)
library(ggplot2)

# Loading some simulated data
attach(simu_data)

# Exposure variables
  X = X_continuous

# Outcome variable
  Y = Y_continuous

# Intermediate variables including mediators
  M = simu_data$M
  
# Choose K (latent factors number) from pca analysis (for example) on potential mediators  
  K = 5
  
detach(simu_data)

# Computing significance values for intermediate variables
# This step uses LFMMs and max-squared tests 
  hdmax2_step1 = run_AS(X = X,
                        Y = Y,
                        M = M,
                        K = K)

# Selecting mediators (ten variables having the lowest p-values)
  mediators_top10 = order(hdmax2_step1$max2_pvalues)[1:10]
  M_10 = M[,mediators_top10]
  
# Ids of selected mediators  
  colnames(M_10)

# Estimating indirect and direct effects of exposure on outcome
  hdmax2_step2 = estimate_effect(object = hdmax2_step1,
                                 m = M_10)
  
# Showing some results
  plot_hdmax2(hdmax2_step2, plot_type = "all_plot")
```


## Bug report / Help

If you encounter a problem, please open a GitHub issue or contact the program developers.


## References 

- Jumentier B, CC Barrot, M Estavoyer, J Tost, B Heude, O Francois, and J Lepeule (2023). High-Dimensional Mediation Analysis: A New Method Applied to Maternal Smoking, Placental DNA Methylation, and Birth Outcomes. Environmental Health Perspectives 131. Publisher: Environmental Health Perspectives, 047011. 
https://doi.org/10.1289/EHP11559.

