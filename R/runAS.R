##' Association Study with both exposure and outcome
##'
##' This function uses lfmm (latent factor mixed models) to estimate
##' the effects of exposures and outcomes on a response matrix.
##' applied to estimate the effects of exposure $X$ on a matrix $M$
##' of potential mediators, and the effect of each marker on outcome $Y$.
##' It uses the covariables matrix $conf$ and $K$ latent factors.
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
##' @param M_type
##' @param X_type
##' @param Y_type
##' @param K an integer for the number of latent factors in the regression model.
##' @param conf set of covariable, must be numeric. No NAs allowed
##' @param diagnostic.plot
##' @param genomic.control correctef pvalue with genomic inflation factor
##' @return an object with the following attributes 
##'   for each association study:
##'
##'  - U, scores matrix for the K latent factors.
##'  
##'  - V, latent factors loadings
##'
##'  - effect.sizes , the effect size matrix for the exposure X and the outcome Y.
##'  
##'  - lambda , use in ridge lfmm
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
##'  - eta0, for each test, the local false discovery rate (FDR) parameter
##'  
##'  - qval, results of max2 test
##'  
##' @details
##' The response variable matrix Y and the explanatory variable are centered.
##' For each argument, missing values must be imputed: no NA allowed. K (number of latent factors) can be estimated
##' with the eigenvalues of a PCA.
##' Possibility of calibrating the scores and pValues by the GIF (Genomic Inflation Factor).
##' See LEA package for more information.
##' @export
##' @author Florence Pittion
##' @examples
##'
##' library(hdmax2)
##'
##' # Exemple 1
##' res <- runAS(X_matrix = example$X, Y_matrix = example$Y, M_matrix = example$M, X_type = "binary", Y_type = "continuous", K = 5)
##'
##' 
runAS = function(X_matrix,
                 Y_matrix,
                 M_matrix, 
                 K,
                 X_type,
                 Y_type,
                 M_type,
                 conf = NULL,
                 diagnostic.plot = F ,
                 genomic.control = T
                 ) {
  res = list()
  
    if(X_type=="continuous"){
    #regression 1: M ~ X
    mod.lfmm1 = LEA::lfmm2(input = M_matrix, 
                           env = X_matrix, 
                           K=5)
    res_ewas1 = LEA::lfmm2.test(mod.lfmm1, 
                                input = M_matrix, 
                                env = X_matrix,
                                genomic.control = genomic.control)
    pval1 = as.double(res_ewas1$pvalues)
    names(pval1) = colnames(M)
    length(pval1)
    U1 = mod.lfmm1@U
    V1 = mod.lfmm1@V
    effect.sizes1 = mod.lfmm1@B
    lambda1 = mod.lfmm1@lambda
    
    zscores1 = res_ewas1$zscores
    fscores1 = res_ewas1$fscores
    #adj_rsquared1 = res_ewas1$adj.r.squared
    gif1 = res_ewas1$gif
    reg1 = list(pval1,
                U1, 
                V1, 
                effect.sizes1,
                lambda1, 
                zscores1,
                fscores1,
                #adj_rsquared1 = res_ewas1$adj.r.squared
                gif1)
    names(reg1) = c("pval","U","V","effect.sizes","lambda","zscores","fscores","gif")
  }
  

  if(X_type=="binary"){
    mod.lfmm1 = LEA::lfmm2(input = M_matrix, 
                           env = X_matrix, 
                           K=5)
    res_ewas1 = LEA::lfmm2.test(mod.lfmm1, 
                                input = M_matrix, 
                                env = X_matrix,
                                genomic.control = genomic.control)
    pval1 = as.double(res_ewas1$pvalues)
    names(pval1) = colnames(M)
    length(pval1)
    U1 = mod.lfmm1@U
    V1 = mod.lfmm1@V
    effect.sizes1 = mod.lfmm1@B
    lambda1 = mod.lfmm1@lambda
    
    zscores1 = res_ewas1$zscores
    fscores1 = res_ewas1$fscores
    #adj_rsquared1 = res_ewas1$adj.r.squared
    gif1 = res_ewas1$gif
    reg1 = list(pval1,
                U1, 
                V1, 
                effect.sizes1, 
                lambda1, 
                zscores1,
                fscores1,
                #adj_rsquared1 = res_ewas1$adj.r.squared
                gif1)
    names(reg1) = c("pval","U","V","effect.sizes","lambda","zscores","fscores","gif")
  }
  
  res[[1]] = reg1
  
  ## TODO  
  # if(X_type=="categorial")
  #     #else ...
  
  if(Y_type=="continuous"){
    mod.lfmm2 = LEA::lfmm2(input = M_matrix, 
                           env = cbind(X_matrix, Y_matrix, conf), 
                           K=5)
    res_ewas2 = LEA::lfmm2.test(mod.lfmm2, 
                                input = M_matrix, 
                                env = cbind(X_matrix, Y_matrix, conf),
                                genomic.control = genomic.control,
                                full = T)
    pval2 = as.double(res_ewas2$pvalues)
    names(pval2) = colnames(M)
    length(pval2)
    U2 = mod.lfmm2@U
    V2 = mod.lfmm2@V
    effect.sizes2 = mod.lfmm2@B
    lambda2 = mod.lfmm2@lambda
    
    zscores2 = res_ewas2$zscores
    fscores2 = res_ewas2$fscores
    #adj_rsquared2 = res_ewas2$adj.r.squared
    gif2 = res_ewas2$gif
    reg2 = list(pval2,
                U2, 
                V2, 
                effect.sizes2, 
                lambda2, 
                zscores2,
                fscores2,
                #adj_rsquared2 = res_ewas2$adj.r.squared
                gif2)
    names(reg2) = c("pval", "U", "V", "effect.sizes", "lambda", "zscores", "fscores", "gif")
  }
    
  
  
   if(Y_type=="binary"){
  mod.lfmm2 = LEA::lfmm2(input = M_matrix, 
                         env = cbind(X_matrix, Y_matrix, conf), 
                         K=5)
  res_ewas2 = LEA::lfmm2.test(mod.lfmm2, 
                              input = M_matrix, 
                              env = cbind(X_matrix, Y_matrix, conf),
                              genomic.control = genomic.control,
                              full = F,
                              linear = F)
  pval2 = as.double(res_ewas2$pvalues)
  names(pval2) = colnames(M)
  length(pval2)
  U2 = mod.lfmm2@U
  V2 = mod.lfmm2@V
  effect.sizes2 = mod.lfmm2@B
  lambda2 = mod.lfmm2@lambda
  
  zscores2 = res_ewas2$zscores
  fscores2 = res_ewas2$fscores
  #adj_rsquared2 = res_ewas2$adj.r.squared
  gif2 = res_ewas2$gif
  reg2 = list(pval2,
              U2, 
              V2, 
              effect.sizes2, 
              lambda2, 
              zscores2,
              fscores2,
              #adj_rsquared2 = res_ewas2$adj.r.squared
              gif2)
  names(reg2) = c("pval", "U", "V", "effect.sizes", "lambda", "zscores", "fscores", "gif")
   }
  
  
  if(Y_type=="surv_Cox"){
    #Y = survival::Surv(OT, status)
    #regression 2: Y ~ M +X
    p = ncol(M)
    pval2 = c()
    for (j in 1:p) {
      #cox_model = survival::coxph(Y ~ M[,j] + X)
      cox_model = survival::coxph(Y ~ M[,j] + X + U1)
      pval2 = c(pval2, summary(cox_model)$coefficients[1,5])
      # beta_est = c(beta_est, summary(cox_model)$coefficients)
    }
    names(pval2)=  colnames(M)
    reg2 = list(pval2#,
                # effect.sizes2,
                # zscores2,
                # fscores2,
                # #adj_rsquared2 = res_ewas2$adj.r.squared
                # gif2)
    )
    names(reg2) = c("pval") #, "effect.sizes", "zscores", "fscores", "gif")
    
  }
  # TODO transform pval in calibrated pval  
  
  # TODO
  # if(Y_type=="surv_AG"){
  #}
  
  res[[2]] = reg2
  
  # max2 test
  max2_pval <- apply(cbind(pval1, pval2), 1, max)^2
  eta0 <- fdrtool::pval.estimate.eta0(max2_pval, diagnostic.plot = diagnostic.plot)
  qval <- fdrtool::fdrtool(max2_pval,statistic = "pvalue", plot = F, verbose = F)$qval
  max2 = list(max2_pval, eta0, qval)
  names(max2) = c("pval", "eta0", "qval")
  res[[3]] = max2
  
  names(res) <- c("mod1", "mod2", "max2")
  return(res)
 
}
