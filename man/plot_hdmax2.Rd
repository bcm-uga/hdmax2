% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdmax2_plot.R
\name{plot_hdmax2}
\alias{plot_hdmax2}
\title{Summary plot for HDMAX2}
\usage{
plot_hdmax2(object, plot_type = "all_plot")
}
\arguments{
\item{object}{results from hdmax2 step 2}

\item{plot_type}{"all_plot" by default generate all proposed plots , elsewhere choose one individual plot among : "plot_ACME", "plot_Prop_Med", "plot_overall_effects", "plot_effect_size"}
}
\value{
Summary plot for ACME, PM, OIE, effects interpretation
}
\description{
This function draw a summary plot of hdmax2 method results
}
\examples{
# Load example dataset
simu_data = hdmax2::simu_data
K = 5
# Run {hdmax2} step 1
hdmax2_step1 = hdmax2::run_AS(
  exposure = simu_data$X_continuous ,
  outcome =  simu_data$Y_continuous,
  M =  simu_data$M,
  K = K
)
# Select mediators
mediators_subset = names(sort(hdmax2_step1$max2_pvalues)[1:10])
mediators_top10 = simu_data$M[, mediators_subset]
# Run {hdmax2} step 2
hdmax2_step2 = hdmax2::estimate_effect(object = hdmax2_step1, 
                                       m = mediators_top10)
# Generate plot
hdmax2::plot_hdmax2(hdmax2_step2, plot_type = "all_plot")

}
\author{
Florence Pittion, Magali Richard
}
