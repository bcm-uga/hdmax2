% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdmax2_helper_functions.R
\name{combp2}
\alias{combp2}
\title{Identification of differentially methylated regions (DMR)}
\usage{
combp2(data, dist.cutoff = 1000, bin.size = 310, seed = 0.01, nCores = 10)
}
\arguments{
\item{data}{A data frame from bed format file with colname name
"V1","V2", "V3","V4","V5",V1 indicate chromosome (1,2,3,...,X,Y),
V2 is chromosome position, V4 is for P value and V5 for name of CpGs.}

\item{dist.cutoff}{Maximum distance in base pair to combine adjacent DMRs.}

\item{bin.size}{bin size for autocorrelation calculation.}

\item{seed}{FDR significance threshold for initial selection of DMR region.}

\item{nCores}{Number of computer cores used in calculation}
}
\value{
Results of the DMRs analysis.
\itemize{
\item result.fdr, table of selected AMRs. For each AMR include chromosomic position, P-value, and FDR
}
}
\description{
Function adapt from the combp function() of the ENmix package
}
\details{
The input should be a data frame with column name V1-V5, indicating chromosome, start position,end position,
pValues and probe names. The function will use a modified comb-p method to identify
differentially methylated regions.
}
\author{
Basile Jumentier
}
