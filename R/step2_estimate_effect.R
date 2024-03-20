##' Estimate effects for a set of mediation markers
##' V2 this version allow to integrate covariates in mediation function 
##'
##' Estimate various quantities for causal mediation analysis for a set of
##' markers, including average causal mediation effects
##' (indirect effect), average direct effects, proportions mediated,
##' and total effect.
##'
##' @param object results from hdmax2 step 2
##' @param m a response variable matrix with n rows and p columns corresponding to mediators selected at step1.
##' Response variables must be encoded as numeric. No NAs allowed.
##' @param boots number of bootstrap
##' @return hdmax2_step2 object
##' 
##'  - ACME, estimation of the average causal mediation effect (the indirect effect)
##'  - ADE, estimation average direct effect
##'  - PM, estimation of the proportion mediated
##'  - TE, estimation of the total effect
##'  
##' Regressions:
##'  - xm, regression X on M
##'  - my, regression M on Y
##' 
##' 
##' Overall effect
##'  - oie, overall indirect effect
##'  - oie_med , oie median
##'  - oie_sd , oie standard deviation
##'  - ote, overall total effect
##'  - ode, overall direct effect
##'  
##' @details
##'
##' We use the mediate function of the mediation package on the set of markers having Q-value lower
##' than the FDR threshold. It estimates their indirect effects and 
##' tests their significance.
##'
##' @export
##' @author Florence Pittion, Magali Richard, Olivier Francois, Basile Jumentier
##' @examples 
##' # Load example dataset
##' simu_data = hdmax2::simu_data
##' K = 5
##' # Run {hdmax2} step 1
##' hdmax2_step1 = hdmax2::run_AS(
##'   X = simu_data$X_continuous,
##'   Y =  simu_data$Y_continuous,
##'   M =  simu_data$M,
##'   K = K
##' )
##' # Select mediators
##' mediators_subset = names(sort(hdmax2_step1$max2_pvalues)[1:10])
##' mediators_top10 = simu_data$M[, mediators_subset]
##' # Run {hdmax2} step 2
##' hdmax2_step2 = hdmax2::estimate_effect(object = hdmax2_step1, 
##'                                        m = mediators_top10)

estimate_effect <- function(object , m, boots = 100) {
  
  if (class(object)!="hdmax2_step1"){
    stop("The object is not of class hdmax2_step1. This function only compute hdmax2 objects generated by hdmax2::run_AS")
  }
  
  X_mat = object$input$X_input
  Y = object$input$Y_input
  
  M = m
  if (is.null(colnames(M))) {
    colnames(M) <- 1:ncol(M)
  }
  
  expo_var_ids =object$input$expo_var_ids
  ncol_var = length(expo_var_ids)
  Y_type = object$input$outcome_var_type
  
  
  effects = list()
  # TODO x = as.dataframe de X
  
  for(expo_var_id in expo_var_ids){
   if (is.vector(X_mat)){
     X = X_mat
   } else if (is.data.frame(X_mat)|| is.matrix(X_mat)){
     X = X_mat[, expo_var_id]
   }
    
    if( ncol_var == 1){
      message("Estimating indirect effect for univariate exposome.")  
      
      if (is.null(object$input$covar)) {
        covars = data.frame(latent_factors = object$modele_1$U)
      } else  {
        covars = data.frame(obs_covar = object$input$covar, latent_factors = object$modele_1$U)
      } 
    } else if( ncol_var > 1){
      message("Estimating indirect effect for multivariate exposome.") 
      extra_expo_vars = expo_var_ids[-which(expo_var_ids %in% expo_var_id)]
      df_extra = X_mat[,which(expo_var_ids %in% extra_expo_vars)]
      if(is.vector(df_extra)){
        df_extra = t(t(df_extra))
      }
      for(col in 1:dim(df_extra)[2]){
        if(typeof(df_extra[,col])=="character"){
          df_extra[,col] = as.factor(df_extra[,col])
          message(paste("Categorial column " , col , " transformed in factors in covariable data frame"))
        } else if (is.factor(df_extra[,col])) {
          message(paste("Categorial column " , col , " is factors in covariable data frame"))
        } else if (typeof(df_extra[,col])=="integer"||typeof(df_extra[,col])== "logical"||typeof(df_extra[,col])== "double"){
          df_extra[,col] = as.numeric(df_extra[,col])
          message(paste("Column " , col , " is Continuous or Binary in covariable data frame"))
        }
      }
      if (is.null(object$input$covar)) {
        covars = data.frame(latent_factors = object$modele_1$U, df_extra = df_extra)
      } else  {
        covars = data.frame(obs_covar = object$input$covar, latent_factors = object$modele_1$U, df_extra = df_extra)
      } 
    }
    
    expo_var_type =  typeof(X)
    
    if (expo_var_type == "character"||is.factor(X)){
      message("The input exposome is categorial")
      if(is.factor(X)==FALSE){
      X_fact = as.factor(X)
      } else {
        X_fact = X
      }
      X_dm = stats::model.matrix(~X_fact)
      message("categorial exposome design matrix transformation")
      X = X_dm[,-1]
      cn = colnames(X)
      
        
      k_effects = list()
      
      for (k in 1:length(cn)) {
        
        # from package mediation
        ACME <- matrix(ncol = 4, nrow = ncol(M))
        ADE <- matrix(ncol = 4, nrow = ncol(M))
        PM <- matrix(ncol = 4, nrow = ncol(M))
        TE <- matrix(ncol = 4, nrow = ncol(M))
        
        # from linear models
        xm <- matrix(ncol = 4, nrow = ncol(M))
        my <- matrix(ncol = 4, nrow = ncol(M))
        
        for (i in 1:ncol(M)) {#numeric
          
          dat.x <- data.frame(Xk = X[,k], Xmk = X[,-k], Mi = M[, i], covars = covars)
          dat.y <- data.frame(Xk = X[,k], Xmk = X[,-k], Mi = M[, i], covars = covars, Y = Y)
          
          mod1 = stats::lm(Mi ~ Xk + ., data = dat.x)
          message(paste0("Generate regression 1 for categorial exposure and mediator ", i))
          
          if(Y_type=="continuous"){
            message(paste0("Generate regression 2 for continuous outcome and mediator ", i))   
            mod2 <- stats::lm(Y ~ Xk + Mi + ., data = dat.y)
            
          } else if(Y_type=="binary"){
            message(paste0("Generate regression 2 for binary outcome and mediator ", i))
            mod2 <- stats::glm(Y ~ Xk + Mi + ., family = "binomial", data = dat.y)
          }
          
          xm[i, ] <- summary(mod1)$coeff[2, ] # effect of X
          my[i, ] <- summary(mod2)$coeff[3, ] # effect of M
          
          med = mediation::mediate(mod1, mod2, treat = "Xk", mediator = "Mi")
          
          ACME[i, ] <- c(med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p)
          ADE[i, ] <- c(med$z0, med$z0.ci[1], med$z0.ci[2], med$z0.p)
          PM[i, ] <- c(med$n0, med$n0.ci[1], med$n0.ci[2], med$n0.p)
          TE[i, ] <- c(med$tau.coef, med$tau.ci[1], med$tau.ci[2], med$tau.p)
        }
        
        ACME <- as.data.frame(ACME)
        ADE <- as.data.frame(ADE)
        PM <- as.data.frame(PM)
        TE <- as.data.frame(TE)
        xm <- as.data.frame(xm)
        my <- as.data.frame(my)
        
        colnames(ACME) <- c("est", "CI_2.5", "CI_97.5", "pval")
        colnames(ADE) <- c("est", "CI_2.5", "CI_97.5", "pval")
        colnames(PM) <- c("est", "CI_2.5", "CI_97.5", "pval")
        colnames(TE) <- c("est", "CI_2.5", "CI_97.5", "pval")
        colnames(xm) <- c("Estimate", "Std.Error", "t.Value", "pValue")
        colnames(my) <- c("Estimate", "Std.Error", "t.Value", "pValue")
        
        ACME$feat <- colnames(M)
        ADE$feat <- colnames(M)
        PM$feat <- colnames(M)
        TE$feat <- colnames(M)
        xm$feat <- colnames(M)
        my$feat <- colnames(M)
        
        # bootstrap
        acme_sum <- matrix(nrow = 1, ncol = boots)
        #covars = as.matrix(covars)

        for (i in 1:ncol(acme_sum)) {
          
          if (is.data.frame(X)||is.matrix(X)){
            samp <- sample(dim(X)[1], replace = T)
          } else if (is.vector(X)) {
            samp <- sample(length(X[,k]), replace = T)
          }
          
          data_samp <- data.frame(Xk = X[samp,k], X = X[samp,-k], m = m[samp, ], covars = covars[samp, ])
          # effect A X -> M
          mod1 <- stats::lm(m ~ Xk + ., data = data_samp)
          A <- t(sapply(summary(mod1), function(x) x$coeff[2, ]))
          A <- data.frame(feat = rownames(A), A)
          
          
          # effect B m -> Y
          if(Y_type=="binary"){
            dat.1 <- data.frame(Y = Y[samp], Xk = X[samp,k], X = X[samp,-k], M=M[samp,], covars = covars[samp,])
            #check that glm converges to include this iteration
            mod2 <- tryCatch(
              stats::glm(Y ~ Xk + .,family = "binomial" , data = dat.1),
              warning = function(w) {
                message(paste0("glm.fit produces warning, iteration ", i, " of bootstrap is not used"))
                return(NULL)
              }
            )
            if (!is.null( mod2)) { 
              B <- as.data.frame(summary(mod2)$coeff[3:(ncol(m) + 2), ])
            }
          } 
          
          if(Y_type=="continuous"){
            dat.1 <- data.frame(Y = Y[samp], Xk = X[samp,k], X = X[samp,-k], M=M[samp,], covars = covars[samp,])
            mod2 <- stats::lm(Y ~ Xk+ ., data = dat.1)
            B <- as.data.frame(summary(mod2)$coeff[3:(ncol(m) + 2), ])
          }
          
          colnames(B) <- c("B", "B_sd", "B_tv", "B_pv")
          colnames(A)[2:5] <- c("A", "A_sd", "A_tv", "A_pv")
          
          ab <- cbind(A, B)
          rownames(ab) <- NULL
          
          # effect A*B
          ab$AB <- ab$A * ab$B
          acme_sum[i] <- sum(ab$AB)
        }
        
        ### Compute ODE and OTE for the given model
        
        if(Y_type == "continuous") {
          message("Computing ODE and OTE for continuous outcome.")
          data_total = data.frame(X = X, Y= Y, covars = covars)
          mod_total_effect = stats::lm(Y ~ . , data =  data_total)
          data_direct = data.frame(X = X, Y= Y, M = M, covars = covars)
          mod_direct_effect = stats::lm(Y ~ ., data =  data_direct)
        }
        
        if(Y_type == "binary") {
          message("Computing ODE and OTE for binary outcome.")
           data_total = data.frame(X = X, Y= Y, covars = covars)
          mod_total_effect = stats::glm(Y ~ . , family = "binomial",  data =  data_total)
           data_direct = data.frame(X = X, Y= Y, M = M ,covars = covars)
          mod_direct_effect = stats::glm(Y ~ . , family = "binomial",  data =  data_direct)
        }
        
        ote = summary(mod_total_effect)$coefficients[2,]
        ode = summary(mod_direct_effect)$coefficients[2,]
        
        oie = as.vector(acme_sum)
        oie_med = median(as.vector(acme_sum))
        oie_sd = sd(as.vector(acme_sum))
        
        tmp = list(ACME=ACME, ADE=ADE, PM=PM, TE=TE, xm=xm, my=my, oie=oie, oie_med=oie_med, oie_sd=oie_sd, ote=ote, ode=ode)
       
        
        k_effects[[paste0("cat_", k)]] = tmp
        
      }
      
     effects[[expo_var_id]] = k_effects
      
      
    } else if (expo_var_type== "integer"||expo_var_type== "logical"||expo_var_type== "double"){
      message("The input exposome is continuous or binary" )
      # boolean transformed as numeric
      X = as.numeric(X)
      
      # To collect data from package mediation
      ACME <- matrix(ncol = 4, nrow = ncol(M))
      ADE <- matrix(ncol = 4, nrow = ncol(M))
      PM <- matrix(ncol = 4, nrow = ncol(M))
      TE <- matrix(ncol = 4, nrow = ncol(M))
      
      # To collect data from linear models
      xm <- matrix(ncol = 4, nrow = ncol(M))
      my <- matrix(ncol = 4, nrow = ncol(M))
      
      for (i in 1:ncol(M)) {#numeric
        dat.x <- data.frame(X = X, Mi = M[, i], covars = covars)
        dat.y <- data.frame(X = X, Mi = M[, i], covars = covars, Y = Y)
        mod1 = stats::lm(Mi ~ X + ., data = dat.x)
        message(paste0("Generate regression 1 for continuous or binary exposure and mediator ", i))
        
        if(Y_type=="continuous"){
          message(paste0("Generate regression 2 for continuous outcome and mediator ", i))   
          mod2 <- stats::lm(Y ~ X + Mi + ., data = dat.y)
        } else if(Y_type=="binary"){
          message(paste0("Generate regression 2 for binary outcome and mediator ", i))
          mod2 <- stats::glm(Y ~ X + Mi + ., family = "binomial", data = dat.y)
        }
        
        xm[i, ] <- summary(mod1)$coeff[2, ] # effect of X
        my[i, ] <- summary(mod2)$coeff[3, ] # effect of M
        
        med = mediation::mediate(mod1, mod2, treat = "X", mediator = "Mi")
        
        ACME[i, ] <- c(med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p)
        ADE[i, ] <- c(med$z0, med$z0.ci[1], med$z0.ci[2], med$z0.p)
        PM[i, ] <- c(med$n0, med$n0.ci[1], med$n0.ci[2], med$n0.p)
        TE[i, ] <- c(med$tau.coef, med$tau.ci[1], med$tau.ci[2], med$tau.p)
      } 
      
      ACME <- as.data.frame(ACME)
      ADE <- as.data.frame(ADE)
      PM <- as.data.frame(PM)
      TE <- as.data.frame(TE)
      xm <- as.data.frame(xm)
      my <- as.data.frame(my)
      
      colnames(ACME) <- c("est", "CI_2.5", "CI_97.5", "pval")
      colnames(ADE) <- c("est", "CI_2.5", "CI_97.5", "pval")
      colnames(PM) <- c("est", "CI_2.5", "CI_97.5", "pval")
      colnames(TE) <- c("est", "CI_2.5", "CI_97.5", "pval")
      colnames(xm) <- c("Estimate", "Std.Error", "t.Value", "pValue")
      colnames(my) <- c("Estimate", "Std.Error", "t.Value", "pValue")
      
      ACME$feat <- colnames(M)
      ADE$feat <- colnames(M)
      PM$feat <- colnames(M)
      TE$feat <- colnames(M)
      xm$feat <- colnames(M)
      my$feat <- colnames(M)
      
      
      
      
      
      ### Compute OIE by bootstrap
      
      # bootstrap
      acme_sum <- matrix(nrow = 1, ncol = boots)
      
      for (i in 1:ncol(acme_sum)) {
        
        if (is.data.frame(X)||is.matrix(X)){
          samp <- sample(dim(X)[1], replace = T)
        } else if (is.vector(X)) {
          samp <- sample(length(X), replace = T)
        }
        data_samp <- data.frame(X = X[samp], m= m[samp, ], covars = covars[samp, ])
        
        # effect A X -> M
        mod1 <- stats::lm(m ~ X + . , data = data_samp)
        A <- t(sapply(summary(mod1), function(x) x$coeff[2, ]))
        A <- data.frame(feat = rownames(A), A)
        
        # effect B m -> Y
        if(Y_type=="binary"){
          dat.1 <- data.frame(Y = Y[samp], X = X[samp], m=m[samp,], covars = covars[samp,])
          #check that glm converges to include this iteration
          mod2 <- tryCatch(
            stats::glm(Y ~ .,family = "binomial" , data = dat.1),
            warning = function(w) {
              message(paste0("glm.fit produces warning, iteration ", i, " of bootstrap is not used"))
              return(NULL)
            }
          )
          
          if (!is.null( mod2)) { 
            B <- as.data.frame(summary(mod2)$coeff[3:(ncol(m) + 2), ])
          }
        } 
        
        if(Y_type=="continuous"){
          dat.1 <- data.frame(Y = Y[samp], X = X[samp], m=m[samp,], covars = covars[samp,])
          mod2 <- stats::lm(Y ~ ., data = dat.1)
          B <- as.data.frame(summary(mod2)$coeff[3:(ncol(m) + 2), ])
        }
        
        colnames(B) <- c("B", "B_sd", "B_tv", "B_pv")
        colnames(A)[2:5] <- c("A", "A_sd", "A_tv", "A_pv")
        
        ab <- cbind(A, B)
        rownames(ab) <- NULL
        
        # effect A*B
        ab$AB <- ab$A * ab$B
        acme_sum[i] <- sum(ab$AB)
      } # end of bootstrap
      
      ### Compute ODE and OTE for the given model
      
      if(Y_type == "continuous") {
        message("Computing ODE and OTE for continuous outcome.")
         data_total = data.frame(X =X, Y= Y, covars = covars)
        mod_total_effect = stats::lm(Y ~ . , data =  data_total)
         data_direct = data.frame(X =X, Y= Y, M =M ,covars = covars)
        mod_direct_effect = stats::lm(Y ~ ., data =  data_direct)
      }
      
      if(Y_type == "binary") {
        message("Computing ODE and OTE for binary outcome.")
         data_total = data.frame(X =X, Y= Y, covars = covars)
        mod_total_effect = stats::glm(Y ~ . , family = "binomial",  data =  data_total)
         data_direct = data.frame(X =X, Y= Y, M = M ,covars = covars)
        mod_direct_effect = stats::glm(Y ~ . , family = "binomial", , data =  data_direct)
      }
      
      ote = summary(mod_total_effect)$coefficients[2,]
      ode = summary(mod_direct_effect)$coefficients[2,]
      
      oie = as.vector(acme_sum)
      oie_med = median(as.vector(acme_sum))
      oie_sd = sd(as.vector(acme_sum))
      
      # oies[[expo_var_id]]= oie
      # oies_med[[expo_var_id]] = oie_med
      # oies_sd[[expo_var_id]] = oie_sd
      # otes[[expo_var_id]] = ote
      # odes[[expo_var_id]] = ode
      
      tmp = list(ACME=ACME, ADE=ADE, PM=PM, TE=TE, xm=xm, my=my, oie=oie, oie_med=oie_med, oie_sd=oie_sd, ote=ote, ode=ode)
      
      
      
      effects[[expo_var_id]] = tmp
      
    }
  }
  
  obj = list(effects = effects,
             input = object$input)
  
  class(obj) = "hdmax2_step2"
  
  return(obj)
}