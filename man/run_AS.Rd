% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step1_run_AS.R
\name{run_AS}
\alias{run_AS}
\title{The function hdmax2::run_AS() evaluates the association
between exposure variables, intermediary variables and the outcome variable
using a latent factor mixed model (LFMM Caye et al. 2019) to
estimate K unobserved latent factors  U.}
\usage{
run_AS(
  exposure,
  outcome,
  M,
  K,
  covar = NULL,
  suppl_covar = NULL,
  each_var_pval = FALSE
)
}
\arguments{
\item{exposure}{An explanatory variable data frame with n rows and d columns.
Each column corresponds to a distinct explanatory variable (exposure).
Continuous and binary variables must be encoded in numeric format. categorical variables are factor objects. The user can use the as.factor function to encode categorical variables, and  levels() and ordered() functions to define the modal order of categorical variables.}

\item{outcome}{An explanatory variable matrix with n rows and 1 columns, corresponds to a vector, which supports both continuous and binary formats.}

\item{M}{Continuous intermediary variables matrix  encompassing potential mediators with n rows and p columns.
Must be encoded as numeric. No NAs allowed.}

\item{K}{an integer for the number of latent factors in the regression model.}

\item{covar}{set of adjustment factors, must be numeric. No NAs allowed}

\item{suppl_covar}{possible supplementary adjustment factors for the second association study (must be nested within the first set of adjustment factors )}

\item{each_var_pval}{A logical to indicate if p-values must be estimated for each exposure variables (each_var_pval = TRUE) in addition to the pvalue of the global model (each_var_pval = FALSE, by default)
Useful to visually check the fit of the estimated proportion of null p-values.}
}
\value{
an object with the following attributes

for first association study (mod1):
\itemize{
\item pValue, estimation of the effects of exposure X and outcome on the matrix M.
\item U, scores matrix for the K latent factors computed from the for first regression
\item zscore, a score matrix for the exposure X and the outcome Y.
\item fscore, a score matrix for the exposure X and the outcome Y.
\item adj_rsquared
\item gif, Genomic inflation factor for X and Y, expressing the deviation of the distribution of the observed test statistic compared to the distribution of the expected test statistic
}

for second association study (mod2):
\itemize{
\item pValue, zscore, fscore,  adj_rsquared, gif
}

results of max2 test:
\itemize{
\item pval, results of max2 test
}

input element:
exposition , outcome and covariates
}
\description{
The function run_AS() evaluates the association between exposure variables X ,
intermediary variables M and the outcome variable Y, using a latent factor mixed model
(LFMM Caye et al. 2019) to estimate K unobserved latent factors  U.
First this function tests the significance of association between
the exposure variables and the potential mediator variables.
Then it tests association between the potential mediator variables and the outcome variable.
Finally it evaluates the significance of the indirect effects by computing
the squared maximum of two series of P-values with max2 test. This rejects
the null-hypothesis that either the effect of X on M, or the effect of M on Y is null.
Optional covariates Z, can be included as observed adjustment factors in the model.
}
\details{
For each argument, missing values must be imputed: no NA allowed. K (number of latent factors) can be estimated
with the eigenvalues of a PCA.
Max2 test The P-value is computed for each markers following this formula
\deqn{pV = max(pVal1, pVal2)^2}
}
\examples{
# Load example dataset
attach(simu_data)
 K = 5
# Run {hdmax2} step 1
hdmax2_step1 = run_AS(exposure = simu_data$X_continuous ,
                      outcome = simu_data$Y_continuous,
                      M = simu_data$M1,
                      K = K)

head(hdmax2_step1$max2_pvalues)
}
\author{
Florence Pittion, Magali Richard, Olivier Francois, Basile Jumentier
}
