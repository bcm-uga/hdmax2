#' Dataset for example
#'
#' @format A list with  objects
#' \itemize{
#'   \item  M methylation matrix, 100 samples and 200 potential mediators
#'   \item  M1 methylation matrix, 100 samples and 200 potential mediators
#'   \item  M2 methylation matrix, 100 samples and 200 potential mediators
#'   \item  X_continuous Continuous Exposure for 100 samples
#'   \item  X_continuous2 Continuous Exposure for 100 samples
#'   \item  X_binary Binary Exposure for 100 samples
#'   \item  X_categorial Categorial Exposure for 100 samples
#'   \item  Y_continuous Continuous Outcome for 100 samples
#'   \item  Y_binary Binary Outcome for 100 samples
#'   \item  age age covariable
#'   \item  gender gender covariable
#' }
#'
"simu_data"

#-----------------------------------------

#' Dataset for helper function example
#'
#' @format a list of 5 objects
#' \itemize{
#'   \item{methylation}{methylation matrix, 100 individuals and 200 probes}
#'   \item{exposure}{Exposure for 100 individuals}
#'   \item{phenotype}{Ouctome for 100 individuals}
#'   \item{annotation}{Annotation for the 200 probes}
#'   \item{covariables}{Covariables for 100 individuals : sex, bmi, age, treatment}
#' }
#'
"helper_ex"