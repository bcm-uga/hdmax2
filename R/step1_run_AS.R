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
##' @param M_matrix a response variable matrix with n rows and p columns.
##' Response variables must be encoded as numeric. No NAs allowed.
##' @param X_matrix Exposure. An explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Exposure). All columns should be in numeric format.
##' @param Y_matrix Outcome. An explanatory variable matrix with n rows and 1 columns.
##' @param X_type Type of exposition, can be either "univariate" (for continuous and binary variables), or "multivariate". Categorial variables must be one-hot encoded as multivariate variables. See helper_functions for example.
##' @param Y_type type of outcome, , can be either "binary" or "continuous"
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
##' @author Florence Pittion
##' @examples
##' # Load example dataset
##' simu_data = hdmax2::simu_data
##' # Run {hdmax2} step 1
##' hdmax2_step1 = run_AS(X_matrix = as.matrix(simu_data$X_continuous) ,
##'                       Y_matrix =  as.matrix(simu_data$Y_continuous),
##'                       M_matrix =  as.matrix(simu_data$M))
##' # Print max2 test results
##' head(hdmax2_step1$max2_pvalues)

run_AS = function(X_matrix,
                  Y_matrix,
                  M_matrix, 
                  K = 5,
                  X_type = "univariate",
                  Y_type = "continuous",
                  covar = NULL,
                  detailed = FALSE,
                  diagnostic.plot = FALSE ,
                  genomic.control = TRUE,
                  effect.sizes = FALSE
) {
  
  res = list()
  
  ##################################
  # Run first regression : M ~ X ###
  ##################################
  
  # Check explanatory variables are in numeric format
  numeric_columns <- apply(X_matrix, 2, is.numeric)
  if (!all(numeric_columns)) {
    non_numeric_column_indices <- which(!numeric_columns)
    stop(paste("Exposome columns", non_numeric_column_indices, "are not in numeric format."))
  } else {
    print("All exposome columns are in numeric format.")
  }
  
  # In univariate situation
  
  if(X_type=="univariate"){
    
    print("Running first regression with univariate explanorty variable.")
    
    mod.lfmm1 = lfmm2_med(input = M_matrix, 
                          env = X_matrix, 
                          K = K,
                          effect.sizes = effect.sizes)
    res_reg1 = lfmm2_med_test(mod.lfmm1, 
                              input = M_matrix, 
                              env = X_matrix,
                              covar = covar,
                              genomic.control = genomic.control)
    pval1 = as.double(res_reg1$pvalues)
    names(pval1) = colnames(M_matrix)
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
  
  if(X_type=="multivariate"){
    
    print("Running first regression with multivariate explatory variables.")
    
    # Computes a global pvalue for regression 1
    
    mod.lfmm1 = lfmm2_med(input = M_matrix, 
                          env = X_matrix, 
                          K = K,
                          effect.sizes = effect.sizes)
    res_reg1 = lfmm2_med_test(mod.lfmm1, 
                              input = M_matrix, 
                              env = X_matrix, 
                              full = TRUE, #parameter to compute a single p-value for the global multivariate model using partial regressions 
                              covar = covar,
                              genomic.control = genomic.control)
    pval1 = res_reg1$pvalues
    names(pval1) = colnames(M_matrix)
    U1 = mod.lfmm1$U
    V1 = mod.lfmm1$V
    zscores1 = res_reg1$zscores
    fscores1 = res_reg1$fscores
    adj_rsquared1 = res_reg1$adj.r.squared
    gif1 = res_reg1$gif
    
    # Computed a single p-value for each explanatory variable
    
    if (detailed == TRUE & X_type=="univariate") {
      stop("Cannot perform detailed analysis for univariate exposome. Detailed analysis is only applicable for multivariate exposomes.")
    }
    
    if(detailed == TRUE){
      
      print("Generating detailed pvalues for each explanatory variable.")
      
      mod.lfmm1 = lfmm2_med(input = M_matrix, 
                            env = X_matrix, 
                            K = K,
                            effect.sizes = effect.sizes)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M_matrix, 
                                env = X_matrix,
                                full = FALSE,
                                covar = covar,
                                genomic.control = genomic.control)
      pvals_1 = as.matrix(res_reg1$pvalues)
      names(pvals_1) = colnames(M_matrix)
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
  
  
  # if(Y_type=="continuous"){
  
  # The model run is actually M ~ X + Y, i.e. independant of the type of Y (continuous or binary)
  
  print("Running second regression.")
  
  Y_matrix = as.numeric(Y_matrix)
  res_reg2 = lfmm2_med_test(mod.lfmm1, #the function will use the latent factors U1 estimated in linear regression 1
                            input = M_matrix, 
                            env = cbind(X_matrix, Y_matrix),
                            covar = covar,
                            genomic.control = genomic.control,
                            full = FALSE)
  pval2 = as.double(res_reg2$pvalues[2,])
  names(pval2) = colnames(M_matrix)
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
  # }
  
  
  
  # if(Y_type=="binary"){
  #   Y_matrix = as.numeric(Y_matrix)
  #   res_reg2 = lfmm2_med_test(mod.lfmm1, 
  #                             input = M_matrix, 
  #                             env = cbind(X_matrix, Y_matrix),
  #                             genomic.control = genomic.control,
  #                             covar = covar,
  #                             full = FALSE,
  #                             linear = FALSE)
  #   pval2 = as.double(res_reg2$pvalues[,2])
  #   names(pval2) = colnames(M_matrix)
  #   zscores2 = res_reg2$zscores
  #   fscores2 = res_reg2$fscores
  #   adj_rsquared2 = res_reg2$adj.r.squared
  #   gif2 = res_reg2$gif
  #   reg2 = list(pval2,
  #               zscores2,
  #               fscores2,
  #               adj_rsquared2, 
  #               gif2)
  #   names(reg2) = c("pval", "zscores", "fscores", "adj_rsquared", "gif")
  # }
  
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
    for (x in 1:dim(X_matrix)[2]){
      max2_pval <- apply(cbind(pvals_1[x,], pval2), 1, max)^2
      names(max2_pval) = colnames(M_matrix)
      max2_detailed[[colnames(X_matrix)[x]]] = max2_pval
    }
    res[[4]] = max2_detailed
  } else {
    print("Not generating max2 pvalues for each explanatory variable.")
    res[[4]] = NA
  }
  
  
  input = list(
    X_matrix,
    Y_matrix,
    X_type,
    Y_type,
    covar
  )
  
  names(input) = c("X_matrix", "Y_matrix", "X_type", "Y_type", "covar")
  
  
  res[[5]] = input
  
  names(res) <- c("modele_1", "modele_2", "max2_pvalues",  "max2_pvalues_detailed", "input")
  
  class(res) = "hdmax2_step1"
  return(res)
  
  
}
