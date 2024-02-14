##' lfmm2_test adapted function for hdmax2 from LEA::lfmm2.test
##' @param object lfmm2Class object
##' @param input a response variable matrix with n rows and p columns
##' @param env An explanatory variable matrix with n rows and d columns.
##' @param covar covariables
##' @param genomic.control correct pvalue with genomic inflation factor
##' @param linear true or false (else is logistic)
##' @param family of logistic reg
##' @param full compute partial regression FALSE/TRUE
##' @return an object with the following attributes 
##' @export
##' @author Florence Pittion
##' @examples 
##' data(simu_data)
##' 
##' mod.lfmm1 = lfmm2_med(input = simu_data$M, 
##' env = simu_data$X_binary, 
##' K = 5,
##' effect.sizes = FALSE)
##' 
##' res_reg1 = lfmm2_med_test(mod.lfmm1, 
##' input = simu_data$M, 
##' env = simu_data$X_binary,
##' covar = cbind(simu_data$age, simu_data$gender),
##' genomic.control = TRUE)
##' 

lfmm2_med_test= function(object, 
                         input,
                         env,
                         covar,
                         full=FALSE,
                         genomic.control=TRUE, 
                         linear=TRUE,
                         family=binomial(link = "logit"))
{
  ## check object
  if (class(object)!="lfmm2"){stop("the object is not lfmm2 type")}
  ## Check input matrix   
  ## LEA  
  if (is.character(input)){
    warning("Reading large input files with 'read.lfmm()' may be slow. See 'data.table::fread()' for fast import.")
    Y <- read.lfmm(input)
    lst.unique <- unique(as.numeric(Y))
    if (9 %in% lst.unique){
      stop("'input' file contains missing data (9's). Use the 'impute()' function to impute them.")
    }
    if (-9 %in% lst.unique){
      stop("'input' file contains missing data (-9's). Use the 'impute()' function to impute them.")
    }
  } else {
    ## Y is an R object       
    if (is.null(input)){
      stop("NULL value for argument 'input'.")
    }
    Y <- as.matrix(input)
    Y[Y == 9] <- NA
    Y[Y == -9] <- NA
    if (anyNA(Y)) {
      stop("The input matrix contains missing values (NA or 9).")
    }
  }
  
  ## Check independent/covariate matrix  
  ## LEA 
  if (is.character(env)){
    X <- read.env(env)
    if (anyNA(X)){
      stop("'env' file contains missing data (NA).")
    }
  } else {
    if (is.null(env)){
      stop("NULL value for argument 'env'.")
    }
    X <- as.matrix(env)
    if (anyNA(X)) {
      stop("The environmental matrix contains NA.")
    }
  }
  
  d <-  ncol(X) #number of environmental variables
  n <-  nrow(X) #number of individuals
  
  if (nrow(Y) != n){
    stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
  }
  
  if (n < d) {
    stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
  }
  
  p <- ncol(Y)
  gif <-  NULL
  p_value <- NULL
  z_score <- NULL
  f_score <- NULL           
  r_squared <- NULL
  
  if (full){
    
    ## a single p-value is returned (f-test)
    ## Check linear models  
    if (linear == FALSE){
      stop("Option full == TRUE is available only for linear models.")
    }
    
    if (is.null(covar)){
      ## partial regression
      mod_Y = lm(Y ~ ., data = data.frame(object$U)) 
      res_Y = mod_Y$residuals
      mod_X = lm(X ~ ., data = data.frame(object$U))
      res_X = mod_X$residuals
      
      mod_lm =  lm(res_Y ~ res_X)
      sm = summary(mod_lm)
      r_squared <- sapply(sm, FUN = function(x) x$adj.r.squared)
      f_score <- sapply(sm, FUN = function(x) x$fstat[1])
      p_value <- sapply(sm, FUN = function(x) pf(x$fstat[1], x$fstat[2], x$fstat[3], lower.tail = F))
    } else {
      mod_Y = lm(Y ~ ., data = data.frame(covar, object$U)) 
      res_Y = mod_Y$residuals
      mod_X = lm(X ~ ., data = data.frame(covar, object$U))
      res_X = mod_X$residuals
      
      mod_lm =  lm(res_Y ~ res_X)
      sm = summary(mod_lm)
      r_squared <- sapply(sm, FUN = function(x) x$adj.r.squared)
      f_score <- sapply(sm, FUN = function(x) x$fstat[1])
      p_value <- sapply(sm, FUN = function(x) pf(x$fstat[1], x$fstat[2], x$fstat[3], lower.tail = F))
    }
  } else {
    
    ## All p-values returned    
    if (linear){
      mod_lm <- lm(Y ~ ., data = data.frame(X, covar,  object$U)) 
      sm <- summary(mod_lm)
      p_value <- sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 4])
      z_score <- as.matrix(sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 3]))
    } else {
      for (j in 1:p) {
        mod_glm <- glm(Y[, j] ~ ., data = data.frame(X, covar, object$U), family = family)
        sm <- summary(mod_glm)
        p_value <- rbind(p_value, sm$coeff[2:(d + 1), 4])
        z_score <- rbind(z_score, sm$coeff[2:(d + 1), 3])
      }
    }
    
  }
  
  
  if (genomic.control){
    
    if (!full){
      
      if (d == 1){
        gif <- median(z_score^2)/qchisq(0.5, df = 1, lower.tail = FALSE)
      } else {
        gif <- apply(z_score^2, 1, median)/qchisq(0.5, df = 1, lower.tail = FALSE)
      }
      p_value <- pchisq(z_score^2/gif, df = 1, lower.tail = FALSE)
      
    } else {
      
      gif <- median(f_score)/qf(0.5, d, n - d - 1, lower.tail = FALSE)
      p_value <- pf(f_score/gif, d, n - d - 1, lower.tail = FALSE)
    }
    
  }
  
  if (!full & anyNA(z_score)) {
    warning("NA in significance values and z-scores. Check the input matrix for non-variable genomic sites.")
  }
  if (full & anyNA(r_squared)) {
    warning("NA in significance values and R-squared. Check the input matrix for non-variable genomic sites.")
  }
  
  res <- list(pvalues = p_value, zscores = z_score, fscores = f_score, adj.r.squared = r_squared,  gif = gif)
  return(res)
}
