% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdmax2_helper_functions.R
\name{AMR_search}
\alias{AMR_search}
\title{Identifying aggregated mediator regions (AMR)}
\usage{
AMR_search(chr, start, end, pval, cpg, ...)
}
\arguments{
\item{chr}{chromosomes}

\item{start}{chromosomal position of markers (start)}

\item{end}{chromosomal position of markers (end)}

\item{pval}{P-values for each markers, from the max2 function}

\item{cpg}{name of each markers}

\item{...}{see help of combp of ENmix package}
}
\value{
\itemize{
\item res, table of selected AMRs. For each AMR include chromosomic position, P-value, and FDR
\item data, matrix of all cpg, with annotation and provided P-values
}
}
\description{
Identify aggregated methylated regions (AMR) from the P-values from function max2 using a modified comb-p method.
Compute the P-value and the FDR for each AMR detected.
}
\details{
The function uses a modified comb-p method to identify
aggregated methylated regions (AMRs).
}
\examples{
data = hdmax2::helper_ex
K=5
## run hdmax2 step1
hdmax2_step1 = hdmax2::run_AS(exposure = data$exposure,
                             outcome = data$phenotype,
                             M = data$methylation,
                             K = K)

##Detecting AMR
chr = data$annotation$chr
start = data$annotation$start
end = data$annotation$end
pval = hdmax2_step1$max2_pvalues
cpg = data$annotation$cpg

res.amr_search = hdmax2::AMR_search(
chr = data$annotation$chr,
start = data$annotation$start,
end = data$annotation$end,
pval = hdmax2_step1$max2_pvalues,
cpg = data$annotation$cpg,
seed = 0.7, #Careful to change this parameter when working with real data
nCores = 2)
res.amr_search$res

}
\author{
Basile Jumentier
}
