##' read lfmm, 
##' @param input.file input file in lfmm object
##' @return an object with the following attributes
##' @importFrom utils read.table
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
##' @importFrom utils read.table
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



##' test_extension 
##' @param name input.file
##' @param extension extension name
##' @return results 
##' @importFrom Smisc getExtension
##' @export
##' @author Olivier Francois
##' @examples 

test_extension <- function(name, extension)
{
  # obtain the extension of name
  ext = getExtension(basename(name))
  
  # if not the correct extension, stop
  if (ext != extension) {
    p = paste("'input_file' format and extension have to be \".", 
              extension, "\" (not \".",ext,"\").", sep="")
    stop(p)
  } 
  
  return(ext);
}