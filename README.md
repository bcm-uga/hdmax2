# hdmax2

## Installation 


```
install.packages (c("survival", "FDRtools", "mediation", "FactoMineR"), depandencies = TRUE)
devtools::install_github("bcm-uga/LEA")
```

***
#### Mac OS arm64 

You need to compile from source a couple of packages `svglite`, `gifski`, `jpeg`, `hdf5r`, `nloptr` and `showtext`.  

Download binaries on CRAN (r-release (arm64)) and locally install them using R CMD INSTALL

- https://cran.r-project.org/web/packages/svglite/index.html
- https://cran.r-project.org/web/packages/gifski/index.html
- https://cran.r-project.org/web/packages/jpeg/index.html
- https://cran.r-project.org/web/packages/hdf5r/index.html
- https://cran.r-project.org/web/packages/nloptr/index.html
- https://cran.r-project.org/web/packages/showtext/index.html
