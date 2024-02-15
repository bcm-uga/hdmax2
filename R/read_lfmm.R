##' read lfmm, 
##' @param input.file input file in lfmm object
##' @return an object with the following attributes
##' @export
##' @author Olivier Francois
##' @examples 


read_lfmm = function(input.file) {
  
  # test arguments
  if(missing(input.file))
    stop("'input.file' argument is missing.")
  else if (!is.character(input.file))
    stop("'input.file' argument has to be of type character.")
  # check extension 
  test_extension(input.file, "lfmm")
  
  return(as.matrix(read.table(input.file)))
}


##' read env 
##' @param input.file input file in lfmm object
##' @return an object with the following attributes
##' @export
##' @author Olivier Francois
##' @examples 

read_env <- function(input.file) {
  
  # test arguments
  if(missing(input.file))
    stop("'input.file' argument is missing.")
  else if (!is.character(input.file))
    stop("'input.file' argument has to be of type character.")
  # check extension 
  test_extension(input.file, "env")
  
  return(as.matrix(read.table(input.file)));
}