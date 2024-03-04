% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdmax2_plot.R
\name{plot_hdmax2}
\alias{plot_hdmax2}
\title{Summary plot for HDMAX2
This function draw a summary plot of hdmax2}
\usage{
plot_hdmax2(object, N_med = 10)
}
\arguments{
\item{object}{results from hdmax2 step 2}

\item{N_med}{number of selected mediators}
}
\value{
Summary plot for ACME, PM, OIE, effects interpretation
}
\description{
Summary plot for HDMAX2
This function draw a summary plot of hdmax2
}
\examples{
# Load example dataset
simu_data = hdmax2::simu_data
# Run {hdmax2} step 1
hdmax2_step1 = hdmax2::run_AS(
  X_matrix = as.matrix(simu_data$X_continuous) ,
  Y_matrix =  as.matrix(simu_data$Y_continuous),
  M_matrix =  as.matrix(simu_data$M)
)
# Select mediators
mediators_subset = names(sort(hdmax2_step1$max2_pvalues)[1:10])
mediators_top10 = simu_data$M[, mediators_subset]
# Run {hdmax2} step 2
hdmax2_step2 = hdmax2::estimate_effect(
   object = hdmax2_step1, 
   m = mediators_top10)
hdmax2::plot_hdmax2(hdmax2_step2, N_med = 10)
}
\author{
Florence Pittion, Magali Richard
}