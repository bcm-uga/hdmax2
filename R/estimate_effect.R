##' Estimate effects for a set of mediation markers
##'
##' Estimate various quantities for causal mediation analysis for a set of
##' markers, including average causal mediation effects
##' (indirect effect), average direct effects, proportions mediated,
##' and total effect.
##'
##' @param M a response variable matrix with n rows and p columns.
##' Response variables must be encoded as numeric. No NAs allowed.
##' @param X Exposure. An explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Exposure).
##' Explanatory variables must be encoded as numeric variables.
##' @param Y Outcome. An explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Outcome).
##' Explanatory variables must be encoded as numeric variables.
##' @param U set of latent factors from runAS function (need include covariable)
##' @param FDR FDR threshold to pass markers in mediation analysis
##' @param sims number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.
##' 10000 is recommended.
##' @param covar covariables
##' @param boots number of bootstrap
##' @param mod2_type second regression type "linear", "logistic", "surv_Cox"
##' @param ... argument of the mediate function from the mediation package
##'
##' @return
##' 
##' Tables of results of mediation analyzes for markers with a Q-value below the FDR threshold.
##' Composition of tables: estimated effect, confidence interval and mediation pValue.
##'  - ACME, estimation of the average causal mediation effect (the indirect effect)
##'  - ADE, estimation average direct effect
##'  - PM, estimation of the proportion mediated
##'  - TE, estimation of the total effect
##'  
##' Regressions:
##'  - xm, regression X on M
##'  - my, regression M on Y
##' 
##' @details
##'
##' We use the mediate function of the mediation package on the set of markers having Q-value lower
##' than the FDR threshold. It estimates their indirect effects and 
##' tests their significance.
##'
##' @export
##' @author 
##' @examples 


estimate_effect <- function(X, Y, m, covar, U , boots = 100, sims = 3,  mod2_type= "linear", ...) {
  
  if (is.null(colnames(m))) {
    colnames(m) <- 1:ncol(m)
  }
  
  M = m
  covar = cbind(covar, U)
  
  ### Compute ACME, ADE, PM and TE from package mediation
  
  # from package mediation
  ACME <- matrix(ncol = 4, nrow = ncol(M))
  ADE <- matrix(ncol = 4, nrow = ncol(M))
  PM <- matrix(ncol = 4, nrow = ncol(M))
  TE <- matrix(ncol = 4, nrow = ncol(M))
  
  # from linear models
  xm <- matrix(ncol = 4, nrow = ncol(M))
  my <- matrix(ncol = 4, nrow = ncol(M))
  
  for (i in 1:ncol(M)) {
    
    dat.x <- data.frame(X = X, Mi = M[, i], covar = covar)
    dat.y <- data.frame(X = X, Mi = M[, i], covar = covar, Y = Y)
    
    # ici cas de deux reg lineaires pour les deux associations
    # TODO refaire pour les autres regressions
    
    
    mod1 <- stats::lm(Mi ~ X + ., data = dat.x)
    
    if(mod2_type=="linear"){
      mod2 <- stats::lm(Y ~ X + Mi + ., data = dat.y)
    }
    
    if(mod2_type=="logistic"){
      mod2 <- stats::glm(Y ~ X + Mi + ., data = dat.y)
    }
    
    # if(mod2_type=="surv_Cox"){
    # mod2 = survival::survreg(Y ~ X + Mi , dist='exponential', data = dat.y)
    # }
    
    
    # for linear models
    xm[i, ] <- summary(mod1)$coeff[2, ] # effect of X
    my[i, ] <- summary(mod2)$coeff[3, ] # effect of M
    
    
    
    med <- mediation::mediate(mod1, mod2, sims = sims, treat = "X", mediator = "Mi", ...)
    
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
    samp <- sample(length(X), replace = T)
    
    if(mod2_type=="logistic"){
      dat.1 <- data.frame(X, m, covar = covar)
      mod1 <- glm(Y[samp] ~ ., data = dat.1[samp, ])
      B <- as.data.frame(summary(mod1)$coeff[3:(ncol(m) + 2), ])
    }
    
    if(mod2_type=="linear"){
      # effet B m -> Y
      dat.1 <- data.frame(X, m, covar = covar)
      mod1 <- lm(Y[samp] ~ ., data = dat.1[samp, ])
      B <- as.data.frame(summary(mod1)$coeff[3:(ncol(m) + 2), ])
      #B <- as.data.frame(summary(mod1)$coeff[3:10, ])
    }
    # effet A X -> M
    mod2 <- lm(m[samp, ] ~ X[samp] + covar[samp, ])
    A <- t(sapply(summary(mod2), function(x) x$coeff[2, ]))
    A <- data.frame(feat = rownames(A), A)
    # A <- separate(A, CpG, c("0", "CpG"), " ")[, -1]
    
    colnames(B) <- c("B", "B_sd", "B_tv", "B_pv")
    colnames(A)[2:5] <- c("A", "A_sd", "A_tv", "A_pv")
    
    ab <- cbind(A, B)
    rownames(ab) <- NULL
    
    # effet A*B
    ab$AB <- ab$A * ab$B
    
    acme_sum[i] <- sum(ab$AB)
  }
  
  ### Compute ODE and OTE for the given model
  
  mod_total_effect = lm(Y~X+covar)
  ote = mod_total_effect$coefficients

  mod_direct_effect = lm(Y~X+mediators_top10+covar)
  ode = mod_direct_effect$coefficients
  
  obj = list(ACME = ACME,
              ADE = ADE,
              PM = PM,
              TE = TE,
              xm = xm,
              my = my,
              oie = as.vector(acme_sum),
              oie_med = median(as.vector(acme_sum)),
              oie_sd = sd(as.vector(acme_sum)),
			  ote = ote,
			  ode = ode
              )
	
  class(obj) = "hdmax2"
   
  return(obj)
  
}