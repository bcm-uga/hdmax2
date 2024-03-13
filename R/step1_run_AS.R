##' Association Study with both exposure and outcome
##'
##' This function uses lfmm (latent factor mixed models) to estimate
##' the effects of exposures and outcomes on a response matrix.
##' applied to estimate the effects of exposure $X$ on a matrix $M$
##' of potential mediators, and the effect of each marker on outcome $Y$.
##' It uses the covariables matrix $conf$ and $K$ latent factors.
##' Test all possible markers to determine potential mediators in the exposure-outcome association.
##' Then compute the squared maximum of two series of P-values with $max2$ test
##' It computes the squared maximum of two series of P-values from the association studie.
##' This rejects the null-hypothesis that either the effect of X on M, or the effect of M on Y is null.
##' @param M a response variable matrix with n rows and p columns.
##' Response variables must be encoded as numeric. No NAs allowed.
##' @param X Exposure. An explanatory variable data frame with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Exposure). 
##' Continuous and binary variables must be encoded in numeric format. categorical variables are factor objects. The user can use the as.factor function to encode categorical variables, and  levels() and ordered() functions to define the modal order of categorical variables.
##' @param Y Outcome. An explanatory variable matrix with n rows and 1 columns, corresponds to a vector, which supports both continuous and binary formats.
##' @param K an integer for the number of latent factors in the regression model.
##' @param covar set of covariable, must be numeric. No NAs allowed
##' @param detailed A logical to indicate if p-values must be estimated for each explanatory variables (detailed = TRUE) in addition to the pvalue of the global model (detailed = FALSE, by default)
##' @param diagnostic.plot if TRUE the histogram of the p-values together
##' with the estimate of eta0 null line is plotted.
##' Useful to visually check the fit of the estimated proportion of null p-values.
##' @param genomic.control correct pvalue with genomic inflation factor
##' @param effect.sizes if effect sizes from lfmm are needed
##' @return an object with the following attributes 
##' 
##' for first association study (mod1):
##'   
##'  - pValue, estimation of the effects of X and Y on the matrix M.
##'
##'  - U, scores matrix for the K latent factors computed from the for first regression
##'  
##'  - zscore, a score matrix for the exposure X and the outcome Y.
##'  
##'  - fscore, a score matrix for the exposure X and the outcome Y.
##'  
##'  - adj_rsquared
##'  
##'  - gif, Genomic inflation factor for X and Y, expressing the deviation of the distribution of the observed test statistic compared to the distribution of the expected test statistic
##'  
##'  
##' for second association study (mod2):
##'  
##'  - pValue, zscore, fscore,  adj_rsquared, gif
##'  
##' results of max2 test:
##'    
##'  - pval, results of max2 test
##'  
##' input element:  
##'  exposition , outcome, matrix  (element and type) and covar
##'  
##'  
##' @details
##' The response variable matrix Y and the explanatory variable are centered.
##' For each argument, missing values must be imputed: no NA allowed. K (number of latent factors) can be estimated
##' with the eigenvalues of a PCA.
##' Possibility of calibrating the scores and pValues by the GIF (Genomic Inflation Factor).
##' See LEA package for more information.
##' Max2 test The P-value is computed for each markers following this formula
##' \deqn{pV = max(pVal1, pVal2)^2}
##' @export
##' @author Florence Pittion, Magali Richard
##' @examples
##' # Load example dataset
##' simu_data = hdmax2::simu_data
##' 
##' # Run {hdmax2} step 1
##' hdmax2_step1 = run_AS(X_matrix = as.matrix(simu_data$X_continuous) ,
##'                       Y_matrix =  as.matrix(simu_data$Y_continuous),
##'                       M_matrix =  as.matrix(simu_data$M),
##'                       K=5)
##' 
##' head(hdmax2_step1$max2_pvalues)

run_AS = function(X,
                  Y,
                  M, 
                  K,
                  covar = NULL,
                  detailed = FALSE,
                  diagnostic.plot = FALSE ,
                  genomic.control = TRUE,
                  effect.sizes = FALSE
) {
  
  ## Check Exposure is a data.frame
  check_argument_exposure(X) 
  
  ## Check Outcome is vector or single column data.frame
  check_argument_outcome(Y)
  
  ## Check Mediator matrix is a numeric matrix
  check_argument_mediators_matrix(M)
  
  ## Check K provided and is integer
  check_K(K)
  
  # Exposure and Outcome before pretreatment
  X_input = X
  Y_input = Y
  
  ## Exposure data frame pretreatment
  # numeric are needed
  if (is.vector(X)){
    expo_var_n = 1
    expo_var_types =  typeof(X)
    expo_var_ids = 1
    if (length(unique(X))<=1){
      stop("Categorial exposome must have at least two levels")
    }
    if (expo_var_types == "character"){
      print("The input exposome is categorial")
      # model matrix transformation
      X = as.factor(X)
      X = model.matrix(~X)
      X = X[,-1]
      new_expo_var_type = typeof(X)
      
    } else if (expo_var_types== "integer"||expo_var_types== "logical"||expo_var_types== "double"){
      print("The input exposome is continuous or binary" )
      X = as.numeric(X)
      new_expo_var_type = typeof(X)
    } 
    
  } else if(is.data.frame(X)){
    expo_var_n = dim(X)[2]
    expo_var_ids = colnames(X)
    expo_var_types = sapply(X, typeof)
    new_expo_var_types = list()
    Xs = NULL
    for(expo_var in 1:expo_var_n) {
      if (expo_var_types[expo_var] == "character"){
        print(paste("The input exposome no ", expo_var," is categorial"))
        # model matrix transformation
        new_X = as.factor(X[,expo_var])
        new_X = model.matrix(~new_X)
        new_X = new_X[,-1]
        new_expo_var_type =  typeof(X)
        
      } else if (expo_var_types[expo_var]== "integer"||expo_var_types[expo_var]== "logical"|| expo_var_types[expo_var]== "double"){
        print(paste("The input exposome no ",expo_var, "is continuous or binary" ))
        new_X = X[,expo_var]
        new_expo_var_type = typeof(X)
      } 
      
      col_name = paste("Var", expo_var, sep="_")
      Xs= cbind(Xs, setNames(new_X,col_name))
      new_expo_var_types[expo_var] = new_expo_var_type
    }
    expo_var_ids = colnames(X)
  } else {
    stop("Unsupported exposure variable type")
  }
  

  ## Outcome pretreatment
  outcome_var_type = NULL
  
  if(is.logical(Y)){
    print("The outcome vector is logical and tranformed in numeric, TRUE become 1 and FALSE become 0.")
    Y = as.matrix(as.numeric(Y))
    outcome_var_type = "binary"
  }
  
  if (all(Y %in% c(0, 1))) {
    if (is.integer(Y)) {
      print("The outcome vector is integer and contains only 0s and 1s.")
      Y = as.matrix(as.double(Y))
      outcome_var_type = "binary"
    } else if (is.double(Y)) {
      print("The outcome vector is numeric and contains only 0s and 1s.")
      Y = as.matrix(Y)
      outcome_var_type = "binary"
    }
  } else if (is.integer(Y)||is.double(Y)) {
    print("The outcome vector is numeric and DON'T contains only 0s and 1s, it is assimilated as continous variable.")
    Y = as.matrix(as.double(Y))
    outcome_var_type = "continuous"
  } else {
    stop("The outcome vector is neither numeric nor logical, therefore it is not supported")
  }
  
  # Exposure and Outcome after pretreatment
  if(expo_var_n == 1){
    X_output = X
  }
  if(expo_var_n > 1){
    X_output = Xs
  }
  Y_output = Y
  
  
  res = list()
  ##################################
  # Run first regression : M ~ X ###
  ##################################
  
   # In univariate situation
  
  if(expo_var_n == 1){
    
    print("Running first regression with univariate explanatory variable.")
    if (expo_var_types == "character"){
      
      mod.lfmm1 = lfmm2_med(input = M, 
                            env = X, 
                            K = K,
                            effect.sizes = effect.sizes)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M, 
                                env = X,
                                covar = covar,
                                full = TRUE, #parameter to compute a single p-value for the global categorial design matrix using partial regressions
                                genomic.control = genomic.control)
    } else if (expo_var_types== "integer"||expo_var_types== "logical"|| expo_var_types== "double"){
      
      mod.lfmm1 = lfmm2_med(input = M, 
                            env = X, 
                            K = K,
                            effect.sizes = effect.sizes)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M, 
                                env = X,
                                covar = covar,
                                genomic.control = genomic.control)
    }
    pval1 = as.double(res_reg1$pvalues)
    names(pval1) = colnames(M)
    U1 = mod.lfmm1$U
    V1 = mod.lfmm1$V
    zscores1 = res_reg1$zscores
    fscores1 = res_reg1$fscores
    adj_rsquared1 = res_reg1$adj.r.squared
    gif1 = res_reg1$gif
    reg1 = list(pval1,
                U1, 
                zscores1,
                fscores1,
                adj_rsquared1,
                gif1)
    names(reg1) = c("pval","U","zscores","fscores", "adj_rsquared", "gif")
  }
  # In multivariate situation
  
  if(expo_var_n > 1){
    X = Xs
    print("Running first regression with multivariate explanatory variables.")
    
    # Computes a global pvalue for regression 1
    
    mod.lfmm1 = lfmm2_med(input = M, 
                          env = X, 
                          K = K,
                          effect.sizes = effect.sizes)
    res_reg1 = lfmm2_med_test(mod.lfmm1, 
                              input = M, 
                              env = X, 
                              full = TRUE, #parameter to compute a single p-value for the global multivariate model using partial regressions 
                              covar = covar,
                              genomic.control = genomic.control)
    pval1 = res_reg1$pvalues
    names(pval1) = colnames(M)
    U1 = mod.lfmm1$U
    V1 = mod.lfmm1$V
    zscores1 = res_reg1$zscores
    fscores1 = res_reg1$fscores
    adj_rsquared1 = res_reg1$adj.r.squared
    gif1 = res_reg1$gif
    
    # Computed a single p-value for each explanatory variable
    
    if (detailed == TRUE & expo_var_n == 1) {
      stop("Cannot perform detailed analysis for univariate exposome. Detailed analysis is only applicable for multivariate exposomes.")
    }
    
    if(detailed == TRUE){
      
      print("Generating detailed pvalues for each explanatory variable.")
      
      mod.lfmm1 = lfmm2_med(input = M, 
                            env = X, 
                            K = K,
                            effect.sizes = effect.sizes)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M, 
                                env = X,
                                full = FALSE,
                                covar = covar,
                                genomic.control = genomic.control)
      pvals_1 = as.matrix(res_reg1$pvalues)
      names(pvals_1) = colnames(M)
    } else {
      pvals_1 = NA
    }
    
    reg1 = list(pval1,
                U1, 
                V1, 
                zscores1,
                fscores1,
                adj_rsquared1, 
                gif1,
                pvals_1 )
    names(reg1) = c("pval","U","V","zscores","fscores", "adj_rsquared","gif","detailed_pval")

  }
  
  res[[1]] = reg1  
  
  #########################################
  ### Run second regression : Y ~ X + M ###
  #########################################
  
  
 
  # The model run is actually M ~ X + Y, i.e. independent of the type of Y (continuous or binary)
  
  print("Running second regression.")
  
  res_reg2 = lfmm2_med_test(mod.lfmm1, #the function will use the latent factors U1 estimated in linear regression 1
                            input = M, 
                            env = cbind(X, Y),
                            covar = covar,
                            genomic.control = genomic.control,
                            full = FALSE)
  pval2 = as.double(res_reg2$pvalues[2,])
  names(pval2) = colnames(M)
  zscores2 = res_reg2$zscores
  fscores2 = res_reg2$fscores
  adj_rsquared2 = res_reg2$adj.r.squared
  gif2 = res_reg2$gif
  reg2 = list(pval2,
              zscores2,
              fscores2,
              adj_rsquared2,
              gif2)
  names(reg2) = c("pval", "zscores", "fscores", "adj_rsquared", "gif")

  
  res[[2]] = reg2
  
  ########################
  ### max-squared test ###
  ########################
  
  print("Running max-squared test.")
  
  # max2 test for global model
  max2_pval <- apply(cbind(pval1, pval2), 1, max)^2
  max2 = max2_pval
  
  res[[3]] = max2
  
  
  if (detailed == TRUE){
    print("Generating max2 pvalues for each explanatory variable.")
    max2_detailed = list()
    for (x in 1:dim(Xs)[2]){
      max2_pval <- apply(cbind(pvals_1[x,], pval2), 1, max)^2
      names(max2_pval) = colnames(M)
      max2_detailed[[colnames(Xs)[x]]] = max2_pval 
    }
    res[[4]] = max2_detailed
  } else {
    print("Not generating max2 pvalues for each explanatory variable.")
    res[[4]] = NA
  }
  
  input = list(
    X_input,
    Y_input,
    X_output,
    Y_output,
    expo_var_types,
    expo_var_ids,
    outcome_var_type,
    covar
  )
  
  names(input) = c("X_input","Y_input","X_output","Y_output", "expo_var_types", "expo_var_ids" , "outcome_var_type", "covar")
  
  res[[5]] = input
  
  names(res) <- c("modele_1", "modele_2", "max2_pvalues",  "max2_pvalues_detailed", "input")
  
  class(res) = "hdmax2_step1"
  return(res)
}


check_argument_exposure = function(argument){
  if(is.data.frame(argument)) {
    print("The exposure argument is a data frame")
    if (ncol(argument) == 1) {
      print("The exposure argument is a data frame with a single column.")
    } else if (ncol(argument) > 1) {
      print("The exposure argument is a data frame with more than one column.")
    }
  } else if (is.vector(argument)) {
    print("The exposure argument is a vector.")
  } else {
    stop("The exposure  is not a data frame,  nor a vector ")
  }
}

check_argument_outcome = function(argument) {
  if (is.vector(argument)) {
    print("The outcome argument is a vector.")
  } else if (is.data.frame(argument)) {
    if (ncol(argument) == 1) {
      print("The outcome argument is a data frame with a single column.")
    } else {
      stop("The outcome data frame must have a single column.")
    }
  } else if (is.matrix(argument)) {
    if (ncol(argument) == 1) {
      print("The outcome matrix has a single column.")
    } else {
      stop("The outcome matrix must have a single column.")
    }
  } else {
    stop("The outcome argument is neither a vector, nor a data frame, nor a matrix with a single column.")
  }
  if (is.numeric(argument)) {
    print("The outcome argument is numeric")
  } else if (is.integer(argument)) {
    print("The outcome argument is integer")
  } else if (is.logical(argument)) {
      print("The outcome argument is logical")
  } else {
    stop("The outcome argument is neither numeric, nor integer, nor logical")
  }
}

check_argument_mediators_matrix = function(argument){
  if(is.matrix(argument)){
    print("Potential mediators matrix is actually a matrix")
  } else {
    stop("Potential mediators matrix must be a matrix")
  }
}

check_K = function(argument){
  if (!is.null(argument)) {
    print(paste("provided K =",argument))
    if(is.integer(argument)) {
      print("K value is integer")
    } else {
      K= as.integer(argument)
      print("K value has been transformed as integer")
    }
  } else {
    stop("K is not provided")
  }
}
