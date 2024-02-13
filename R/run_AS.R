##' Association Study with both exposure and outcome
##'
##' This function uses lfmm (latent factor mixed models) to estimate
##' the effects of exposures and outcomes on a response matrix.
##' applied to estimate the effects of exposure $X$ on a matrix $M$
##' of potential mediators, and the effect of each marker on outcome $Y$.
##' It uses the covariables matrix $conf$ and $K$ latent factors.
##' Compute the squared maximum of two series of P-values
##' Test all possible markers to determine potential mediators in the exposure-outcome association.
##' It computes the squared maximum of two series of P-values from the association studie.
##' This rejects the null-hypothesis that either the effect of X on M, or the effect of M on Y is null.
##'
##'
##' @param M_matrix a response variable matrix with n rows and p columns.
##' Response variables must be encoded as numeric. No NAs allowed.
##' @param X_matrix Exposure. An explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Exposure).
##' Explanatory variables must be encoded as numeric variables.
##' @param Y_matrix Outcome. An explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Outcome).
##' Explanatory variables must be encoded as numeric variables.
##' @param M_type type of potential mediators matrix "methylation", "transcriptome"
##' @param X_type type of exposition "univariate", "multivariate"
##' @param Y_type type of outcome "binary', "continuous"
##' @param K an integer for the number of latent factors in the regression model.
##' @param covar set of covariable, must be numeric. No NAs allowed
##' @param diagnostic.plot if TRUE the histogram of the p-values together
##' with the estimate of eta0 null line is plotted.
##' Useful to visually check the fit of the estimated proportion of null p-values.
##' @param genomic.control correct pvalue with genomic inflation factor
##' @return an object with the following attributes 
##'   for each association study:
##'
##'  - U, scores matrix for the K latent factors.
##'  
##'  - V, latent factors loadings
##'
##'  
##'  - gif, Genomic inflation factor for X and Y, expressing the deviation of the distribution of the observed test statistic compared to the distribution of the expected test statistic
##'  
##'  - pValue, estimation of the effects of X and Y on the matrix M.
##' 
##'  - zscore, a score matrix for the exposure X and the outcome Y.
##'  
##'  - fscore, a score matrix for the exposure X and the outcome Y.
##'
##'    results of max2 test:
##'    
##'  - pval, results of max2 test
##'  
##'  
##' @details
##' The response variable matrix Y and the explanatory variable are centered.
##' For each argument, missing values must be imputed: no NA allowed. K (number of latent factors) can be estimated
##' with the eigenvalues of a PCA.
##' Possibility of calibrating the scores and pValues by the GIF (Genomic Inflation Factor).
##' See LEA package for more information.
##' Max2 test The P-value is computed for each markers following this formula
##'
##' \deqn{pV = max(pVal1, pVal2)^2}
##' 
##' This quantity eta0, i.e. the proportion eta0 of null p-values in a given vector of p-values,
##' is an important parameter when controlling the false discovery rate (FDR).
##' A conservative choice is eta0 = 1 but a choice closer to the true value will
##' increase efficiency and power - see Benjamini and Hochberg (1995, 2000) and Storey (2002) for details.
##' We use the fdrtool package to transform pValues into qValues,
##' which allows us to control the FDR.
##' @export
##' @author Florence Pittion
##' @examples 

run_AS = function(X_matrix,
                  Y_matrix,
                  M_matrix, 
                  K,
                  X_type,
                  Y_type,
                  M_type,
                  covar,
                  multivariate = FALSE,
                  diagnostic.plot = FALSE ,
                  genomic.control = TRUE,
                  effect.sizes = FALSE
) {
  res = list()
  
  if(X_type=="continuous"){
    #regression 1: M ~ X
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
                V1, 
                zscores1,
                fscores1,
                adj_rsquared1,
                gif1)
    names(reg1) = c("pval","U","V","zscores","fscores", "adj_rsquared", "gif")
  }
  
  if(X_type=="binary"){
    #transfo X binary in continuous for lm
    X_matrix = as.numeric(X_matrix)
    #regression 1: M ~ X
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
                V1, 
                zscores1,
                fscores1,
                adj_rsquared1,
                gif1)
    names(reg1) = c("pval","U","V","zscores","fscores", "adj_rsquared", "gif")
  }
  
  if(X_type=="categorial"){
    X_matrix = droplevels(X_matrix)
    X_matrix = mltools::one_hot(data.table::as.data.table(X_matrix))
    
    if(multivariate){
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
      pval1 = as.matrix(res_reg1$pvalues)
      names(pval1) = colnames(M_matrix)
      U1 = mod.lfmm1$U
      V1 = mod.lfmm1$V
      zscores1 = res_reg1$zscores
      fscores1 = res_reg1$fscores
      adj_rsquared1 = res_reg1$adj.r.squared
      gif1 = res_reg1$gif
      reg1 = list(pval1,
                  U1, 
                  V1, 
                  zscores1,
                  fscores1,
                  adj_rsquared1, 
                  gif1)
      names(reg1) = c("pval","U","V","zscores","fscores", "adj_rsquared","gif")
    } 
    if (!multivariate) {
      
      mod.lfmm1 = lfmm2_med(input = M_matrix, 
                            env = X_matrix, 
                            K = K,
                            effect.sizes = effect.sizes)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M_matrix, 
                                env = X_matrix,
                                full = TRUE,
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
      reg1 = list(pval1,
                  U1, 
                  V1, 
                  zscores1,
                  fscores1,
                  adj_rsquared1, 
                  gif1)
      names(reg1) = c("pval","U","V","zscores","fscores", "adj_rsquared","gif")
      
    }
  }
  
  res[[1]] = reg1  
  
  if(Y_type=="continuous"){
    Y_matrix = as.numeric(Y_matrix)
    res_reg2 = lfmm2_med_test(mod.lfmm1, 
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
  }
  
  
  
  if(Y_type=="binary"){
    Y_matrix = as.numeric(Y_matrix)
    res_reg2 = lfmm2_med_test(mod.lfmm1, 
                              input = M_matrix, 
                              env = cbind(X_matrix, Y_matrix),
                              genomic.control = genomic.control,
                              covar = covar,
                              full = FALSE,
                              linear = FALSE)
    pval2 = as.double(res_reg2$pvalues[,2])
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
  }
  
  res[[2]] = reg2
  
  # max2 test
  if (multivariate){
    max2 = c()
    for (x in 1:dim(pval1)[1]){
    max2_pval <- apply(cbind(pval1[x,], pval2), 1, max)^2
    max2 = rbind(max2, max2_pval)
    #rownames(max2) = colnames(X)
    }
    rownames(max2) = colnames(X_matrix)
    colnames(max2) = colnames(M_matrix)
    
  }else {
  max2_pval <- apply(cbind(pval1, pval2), 1, max)^2
  max2 = max2_pval
  }
  res[[3]] = max2
  
  names(res) <- c("mod1", "mod2", "max2")
  return(res)
  
}
