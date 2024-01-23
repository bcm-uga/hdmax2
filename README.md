# hdmax2

## Conda Installation for Mac OS arm64


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


