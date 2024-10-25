run_AS_surv_aalen = function(exposure,
                           survival_time,
                           censoring_status,
                           M, 
                           K,
                           covar = NULL,
                           suppl_covar = NULL,
                           each_var_pval = FALSE
) {
  
  ## Check exposure is a data.frame
  check_argument_exposure(exposure) 
  
  # ## Check outcome is vector or single column data.frame
  # check_argument_outcome(outcome)
  
  ## Check survival_time
  
  ## Check censoring_status
  
  ## Check Mediator matrix is a numeric matrix
  check_argument_mediators_matrix(M)
  
  ## Check K provided and is integer
  check_K(K)
  
  if(!is.null(covar)){
    check_covar(covar)
  }
  
  if(!is.null(suppl_covar)){
    check_covar(suppl_covar)
  }
  
  # Exposure and outcome before pretreatment
  exposure_input = exposure
  survival_time_input = survival_time
  censoring_status_input = censoring_status
  
  ## Exposure data frame pretreatment
  # numeric are needed
  if (is.vector(exposure)){
    expo_var_n = 1
    expo_var_types =  typeof(exposure)
    expo_var_ids = "univariate"
    if (length(unique(exposure))<=1){
      stop("Categorial exposome must have at least two levels")
    }
    if (expo_var_types == "character"){
      #message("The input exposome is categorial")
      # model matrix transformation (-1 column to avoid colinearity)
      exposure = as.factor(exposure)
      exposure = model.matrix(~exposure)
      exposure = exposure[,-1]
      new_expo_var_type = typeof(exposure)
      
    } 
    else if (expo_var_types== "integer"||expo_var_types== "logical"||expo_var_types== "double"){
      #message("The input exposome is continuous or binary" )
      exposure = as.numeric(exposure)
      new_expo_var_type = typeof(exposure)
    } 
    
  } else if (is.factor(exposure)){
    expo_var_n = 1
    expo_var_types =  typeof(exposure)
    expo_var_ids = "univariate"
    #message("The input exposome is categorial")
    # model matrix transformation (-1 column to avoid colinearity)
    exposure = model.matrix(~exposure)
    exposure = exposure[,-1]
    new_expo_var_type = typeof(exposure)
    
  } else if(is.data.frame(exposure)){
    expo_var_n = dim(exposure)[2]
    expo_var_ids = colnames(exposure)
    expo_var_types = sapply(exposure, typeof)
    new_expo_var_types = list()
    exposures = c()
    for(expo_var in 1:expo_var_n) {
      if (expo_var_types[expo_var] == "character"){
        #message(paste("The input exposome no ", expo_var," is categorial"))
        # model matrix transformation
        new_exposure = as.factor(exposure[,expo_var])
        new_exposure = model.matrix(~new_exposure)
        new_exposure = new_exposure[,-1]
        new_expo_var_type =  typeof(new_exposure)
        
      } else if (is.factor(exposure)){ 
        #message(paste("The input exposome no ", expo_var," is categorial"))
        # model matrix transformation
        new_exposure = model.matrix(~new_exposure)
        new_exposure = new_exposure[,-1]
        new_expo_var_type =  typeof(new_exposure)
      } else if (expo_var_types[expo_var]== "integer"||expo_var_types[expo_var]== "logical"|| expo_var_types[expo_var]== "double"){
        #message(paste("The input exposome no ",expo_var, "is continuous or binary" ))
        new_exposure = exposure[,expo_var]
        new_expo_var_type = typeof(new_exposure)
      } 
      
      col_name = paste("Var", expo_var, sep="_")
      exposures= cbind(exposures, stats::setNames(new_exposure,col_name))
      new_expo_var_types[expo_var] = new_expo_var_type
    }
    expo_var_ids = colnames(exposure)
  } else {
    stop("Unsupported exposure variable type")
  }
  
  
  ## outcome pretreatment
  
  
  
  # exposure and outcome after pretreatment
  # if(expo_var_n == 1){
  #   if(is.data.frame(X)){
  #     exposure_output = exposures
  #   } else {
  #     exposure_output = exposure
  #   }
  # }
  # if(expo_var_n > 1){
  #   exposure_output = exposures
  # }
  # outcome_output = outcome
  
  
  res = list()
  ##################################
  # Run first regression : M ~ X ###
  ##################################
  
  # In univariate situation
  
  if(expo_var_n == 1){
    
    message("Running first regression with univariate exposure variable.")
    if (expo_var_types == "character"|| is.factor(exposure_input)){
      
      mod.lfmm1 = lfmm2_med(input = M, 
                            env = exposure, 
                            K = K,
                            effect.sizes = FALSE)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M, 
                                env = exposure,
                                covar = covar,
                                full = TRUE, #parameter to compute a single p-value for the global categorial design matrix using partial regressions
                                genomic.control = TRUE)
      if(each_var_pval == TRUE){
        
        #message("Generating detailed pvalues for each explanatory variable.")
        
        mod.lfmm1 = lfmm2_med(input = M, 
                              env = exposure, 
                              K = K,
                              effect.sizes = FALSE)
        res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                  input = M, 
                                  env = exposure,
                                  full = FALSE,
                                  covar = covar,
                                  genomic.control = TRUE)
        pvals_1 = as.matrix(res_reg1$pvalues)
        names(pvals_1) = colnames(M)
      } else {
        pvals_1 = NA
      }
    } else if (is.vector(exposure) && (expo_var_types== "integer"||expo_var_types== "logical"|| expo_var_types== "double")){
      
      mod.lfmm1 = lfmm2_med(input = M, 
                            env = exposure, 
                            K = K,
                            effect.sizes = FALSE)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M, 
                                env = exposure,
                                covar = covar,
                                genomic.control = TRUE)
      
      
      
    }else if (is.data.frame(exposure) && (expo_var_types== "integer"||expo_var_types== "logical"|| expo_var_types== "double")){
      
      mod.lfmm1 = lfmm2_med(input = M, 
                            env = exposures, 
                            K = K,
                            effect.sizes = FALSE)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M, 
                                env = exposures,
                                covar = covar,
                                genomic.control = TRUE)
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
    exposure = exposures
    message("Running first regression with multivariate exposure variables.")
    
    # Computes a global pvalue for regression 1
    
    mod.lfmm1 = lfmm2_med(input = M, 
                          env = exposure, 
                          K = K,
                          effect.sizes = FALSE)
    res_reg1 = lfmm2_med_test(mod.lfmm1, 
                              input = M, 
                              env = exposure, 
                              full = TRUE, #parameter to compute a single p-value for the global multivariate model using partial regressions 
                              covar = covar,
                              genomic.control = TRUE)
    pval1 = res_reg1$pvalues
    names(pval1) = colnames(M)
    U1 = mod.lfmm1$U
    V1 = mod.lfmm1$V
    zscores1 = res_reg1$zscores
    fscores1 = res_reg1$fscores
    adj_rsquared1 = res_reg1$adj.r.squared
    gif1 = res_reg1$gif
    
    # Computed a single p-value for each explanatory variable
    
    if (each_var_pval == TRUE & expo_var_n == 1) {
      stop("Cannot perform detailed analysis for univariate exposome. Detailed analysis is only applicable for multivariate exposomes.")
    }
    
    if(each_var_pval == TRUE){
      
      message("Generating detailed pvalues for each explanatory variable.")
      
      mod.lfmm1 = lfmm2_med(input = M, 
                            env = exposure, 
                            K = K,
                            effect.sizes = FALSE)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M, 
                                env = exposure,
                                full = FALSE,
                                covar = covar,
                                genomic.control = TRUE)
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
    names(reg1) = c("pval","U","V","zscores","fscores", "adj_rsquared","gif","each_var_pval")
    
  }
  
  res[[1]] = reg1  
  
  
  
  #########################################
  ### Run second regression : Y ~ X + M ###
  #########################################
  
  message("Running second regression.")
  if(!is.null(suppl_covar)){
    covars = cbind(covar, suppl_covar)
  } else {
    covars = covar
  }
  
  if(expo_var_n == 1){
    if(is.vector(exposure)||is.factor(exposure)||is.matrix(exposure)){
      res_reg2 = surv_aalen_med_test(mod.lfmm1,
                                   M = M,
                                   survival_time = survival_time,
                                   censoring_status = censoring_status,
                                   exposure = exposure,
                                   covar = covars)
    }else if(is.data.frame(exposure)){
      res_reg2 = surv_aalen_med_test(mod.lfmm1,
                                   M = M,
                                   survival_time = survival_time,
                                   censoring_status = censoring_status,
                                   exposure = exposures,
                                   covar = covars)
    }
  }else if(expo_var_n > 1){
    res_reg2 = surv_aalen_med_test(mod.lfmm1,
                                 M = M,
                                 survival_time = survival_time,
                                 censoring_status = censoring_status,
                                 exposure = exposures,
                                 covar = covars)
  }
  
  pval2 = as.double(res_reg2$pvalues)
  names(pval2) = colnames(M)
  #zscores2 = res_reg2$zscores
  
  reg2 = list(pval2 = pval2)
  # ,
  #             zscores2 = zscores2
  # )
  names(reg2) = c("pval")#, "zscores")
  
  
  res[[2]] = reg2
  
  ########################
  ### max-squared test ###
  ########################
  
  message("Running max-squared test.")
  
  # max2 test for global model
  max2_pval <- apply(cbind(pval1, pval2), 1, max)^2
  max2 = max2_pval
  
  res[[3]] = max2
  
  if (each_var_pval == TRUE){
    #message("Generating max2 pvalues for each explanatory variable.")
    max2_each_var_pval = list()
    if(expo_var_n == 1){
      for (x in 1:dim(exposure)[2]){
        max2_pval <- apply(cbind(pvals_1[x,], pval2), 1, max)^2
        names(max2_pval) = colnames(M)
        max2_each_var_pval[[colnames(exposures)[x]]] = max2_pval 
      }
      res[[4]] = max2_each_var_pval
    } else {
      for (x in 1:dim(exposures)[2]){
        max2_pval <- apply(cbind(pvals_1[x,], pval2), 1, max)^2
        names(max2_pval) = colnames(M)
        max2_each_var_pval[[colnames(exposures)[x]]] = max2_pval 
      }
      res[[4]] = max2_each_var_pval
    } 
  } else {
    #message("Not generating max2 pvalues for each explanatory variable.")
    res[[4]] = NA
  }
  
  input = list(
    exposure_input,
    expo_var_types,
    expo_var_ids,
    covar, 
    suppl_covar,
    survival_time_input, 
    censoring_status_input#,
    #survival_distribution
  )
  
  names(input) = c("exposure_input", "expo_var_types", "expo_var_ids" , "covar", "suppl_covar", "survival_time_input", "censoring_status_input")#, "survival_distribution")
  
  
  res[[5]] = input
  
  names(res) <- c("AS_1", "AS_2", "max2_pvalues",  "max2_each_var_pvalues", "input")
  
  class(res) = "hdmax2_step1"
  return(res)
}






check_argument_exposure = function(argument){
  if(is.data.frame(argument)) {
    #message("The exposure argument is a data frame")
    if (ncol(argument) == 1) {
      #message("The exposure argument is a data frame with a single column.")
    } else if (ncol(argument) > 1) {
      #message("The exposure argument is a data frame with more than one column.")
    }
  } else if (is.vector(argument)) {
    #message("The exposure argument is a vector.")
  } else if (is.factor(argument)) {
    #message("The exposure argument is a factor.")
  } else {
    stop("The exposure  is not a data frame,  nor a vector ")
  }
}

check_argument_outcome = function(argument) {
  if (is.vector(argument)) {
    #message("The outcome argument is a vector.")
  } else if (is.data.frame(argument)) {
    if (ncol(argument) == 1) {
      #message("The outcome argument is a data frame with a single column.")
    } else {
      stop("The outcome data frame must have a single column.")
    }
  } else if (is.matrix(argument)) {
    if (ncol(argument) == 1) {
      #message("The outcome matrix has a single column.")
    } else {
      stop("The outcome matrix must have a single column.")
    }
  } else {
    stop("The outcome argument is neither a vector, nor a data frame, nor a matrix with a single column.")
  }
  if (is.numeric(argument)) {
    #message("The outcome argument is numeric")
  } else if (is.integer(argument)) {
    #message("The outcome argument is integer")
  } else if (is.logical(argument)) {
    #message("The outcome argument is logical")
  } else {
    stop("The outcome argument is neither numeric, nor integer, nor logical")
  }
}

check_argument_mediators_matrix = function(argument){
  if(is.matrix(argument)){
    #message("Potential mediators matrix is actually a matrix")
  } else {
    stop("Potential mediators matrix must be a matrix")
  }
}

check_K = function(argument){
  if (!is.null(argument)) {
    #message(paste("provided K =",argument))
    if(is.integer(argument)) {
      #message("K value is integer")
    } else {
      K= as.integer(argument)
      #message("K value has been transformed as integer")
    }
  }
  else {
    stop("K is not provided")
  }
}

check_covar = function(argument){
  if (is.data.frame(argument)){
    #message("Adjutment factors is data frame")
  }else if(is.matrix(argument)){
    #message("Adjutment factors is matrix")
  } else{
    stop("Adjutment factors must be a data frame or a matrix")
  }
  for (i in 1:dim(argument)[2]){
    if(!is.numeric(argument[,1])){
      stop("adjusment factors must be numeric")
    }
  }
}

