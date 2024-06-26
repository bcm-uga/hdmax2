% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step2_estimate_effect.R
\name{estimate_effect}
\alias{estimate_effect}
\title{The function hdmax2::estimate_effect() takes as input an object
hdmax2_step1 and a list of potential mediators MS to be analyzed in subsequent steps.}
\usage{
estimate_effect(object, m, boots = 1000, ...)
}
\arguments{
\item{object}{results from hdmax2 step 2}

\item{m}{a response variable matrix with n rows and p columns corresponding to mediators selected at step1.
Response variables must be encoded as numeric. No NAs allowed.}

\item{boots}{number of bootstrap}

\item{...}{arguments for inner functions}
}
\value{
hdmax2_step2 object
\itemize{
\item ACME, estimation of the average causal mediation effect (the indirect effect)
\item ADE, estimation average direct effect
\item PM, estimation of the proportion mediated
\item TE, estimation of the total effect
}

Regressions:
\itemize{
\item xm, regression X on M
\item my, regression M on Y
}

Overall effect
\itemize{
\item oie, overall indirect effect
\item oie_med , oie median
\item oie_sd , oie standard deviation
\item ote, overall total effect
\item ode, overall direct effect
}
}
\description{
For each univariate exposure variable and the subset of mediators MS,
the hdmax2::estimate_effect() function computes several estimates
to evaluate the indirect effects in the path between exposure variables
and the outcome variable. Initially, this function assesses each mediator
variable MS_j individually and computes causal measures of interest
such as (i) the Average Causal Mediated Effect (ACME, corresponding to
the indirect effect) and (ii) the Proportion Mediated (PM). The ACME
differs from the Average Direct Effect (ADE), which represents the
unmediated effect, and from the Total Effect (TE) which is equal to the sum
of direct and indirect effect. PM corresponds to the proportion of
the total effect that is mediated by the mediator (ratio of the indirect
effect to the total effect). ACME and PM are computed by the
mediation::mediate() function of the package mediation, that
automatically detects the type of statistical model used in the mediation
analysis (Tingley et al. 2014). The function mediation::mediate()
calculates uncertainty estimates by a quasi-Bayesian Monte Carlo approach
described (Imai et al. 2010). In addition, it estimates the intermediary
effect sizes  a_j and b_j and their standard deviations. Eventually,
hdmax2 calculates an Overall Indirect Effect (OIE) from a single model
that includes all mediators MS simultaneously. The OIE corresponds
to the sum of the indirect effect associated with all mediators.
The confidence interval (CI) of the OIE is estimated by a bootstrap
approach. Along with the OIE, hdmax2 estimates the Overall Total
Effect (OTE) corresponding to the effect of exposure variables on
the outcome variable, as well as the Overall Direct Effect (ODE)
corresponding to the effect of exposure variables on the outcome variable
when the mediators MS are included in the model.
}
\details{
We use the mediate function of the mediation package on the set of markers having Q-value lower
than the FDR threshold. It estimates their indirect effects and
tests their significance.
}
\examples{
# Load example dataset
attach(simu_data)
K = 5
# Run {hdmax2} step 1
hdmax2_step1 = hdmax2::run_AS(
  exposure = simu_data$X_binary,
  outcome =  simu_data$Y_continuous,
  M =  simu_data$M1,
  K = K
)
# Select mediators
mediators_subset = names(sort(hdmax2_step1$max2_pvalues)[1:10])
mediators_top10 = simu_data$M1[, mediators_subset]
# Run {hdmax2} step 2
hdmax2_step2 = hdmax2::estimate_effect(object = hdmax2_step1, 
                                       m = mediators_top10)
}
\author{
Florence Pittion, Magali Richard, Olivier Francois, Basile Jumentier
}
