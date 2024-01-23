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
  # return(list(calibrated_pval1 = pval1, 
  #             calibrated_pval2 = pval2, 
  #             U1 = U1,
  #             V1 = V1,
  #             effect.sizes1 = effect.sizes1,
  #             lambda1 = lambda1,
  #             zscores1 = zscores1,
  #             fscores1 = fscores1,
  #             #adj_rsquared1 = adj.r.squared1,
  #             gif1 = gif1,
  #             U2 = U2,
  #             V2 = V2,
  #             effect.sizes2 = effect.sizes2,
  #             lambda2 = lambda2,
  #             zscores2 = zscores2,
  #             fscores2 = fscores2,
  #             #adj_rsquared2 = adj.r.squared2,
  #             gif2 = gif2,
  #             max2_pval = max2_pval,
  #             max2_eta0 = eta0,
  #             max2_qval = qval
  # ))
}
