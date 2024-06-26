% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdmax2_lfmm2_functions.R
\name{lfmm2_med}
\alias{lfmm2_med}
\title{Fitting Latent Factor Mixed Models (Least squares algorithm)}
\usage{
lfmm2_med(input, env, K, lambda = 1e-05, effect.sizes = FALSE)
}
\arguments{
\item{input}{Continuous intermediary variables matrix  encompassing potential mediators with n rows and p columns.}

\item{env}{An explanatory variable matrix with n rows and d columns.}

\item{K}{latent factor number}

\item{lambda}{ridge penalization parameter}

\item{effect.sizes}{true or false to obtain effect sizes}
}
\value{
an object with the following attributes
U
V
}
\description{
Latent Factor Mixed Models (LFMMs) are factor regression models.
The lfmm2 function estimates latent factors based on an exact least-squares
approach.
}
\examples{
data(simu_data)
K = 5
mod.lfmm1 = lfmm2_med(input = simu_data$M1, 
env = simu_data$X_binary, 
K = K,
effect.sizes = FALSE)
}
\author{
Florence Pittion, Magali Richard, Olivier Francois
}
