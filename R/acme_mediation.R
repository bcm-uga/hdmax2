##' Estimate effects for a set of mediation markers
##'
##' Estimate various quantities for causal mediation analysis for a set of
##' markers, including average causal mediation effects
##' (indirect effect), average direct effects, proportions mediated,
##' and total effect.
##'
##' @param qval set of qValues from max2 function
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
##' @author Basile Jumentier
##' @examples
##'
##' library(hdmax2)
##'
##'
##'
##' #
##' res <- runAS(X_matrix = example$X, Y_matrix = example$Y, M_matrix = example$M, X_type = "binary", Y_type = "continuous", K = 5)
##'
##' res <- acme_mediation(qval = res$max2$qval,
##'                             X = example$X,
##'                             Y = example$Y,
##'                             M = example$M,
##'                             U = res$mod1$U, sims = 3,
##'                             FDR = 0.5)
##'
##'

acme_mediation <- function(qval, X, Y, M, covar = NULL, U = NULL, FDR = 0.1, sims = 3, mod2_type, ...) {
  
  if (is.null(colnames(M))) {
    colnames(M) <- 1:ncol(M)
  }
  
  M <- M[, qval <= FDR]
  
  ##' Tables of results of mediation analyzes for markers with a Q-value below the FDR threshold.
  ##' Composition of tables: estimated effect, confidence interval and mediation pValue.
  ##'  - ACME, estimation of the average causal mediation effect (the indirect effect)
  ##'  - ADE, estimation average direct effect
  ##'  - PM, estimation of the proportion mediated
  ##'  - TE, estimation of the total effect
  # from package mediation
  ACME <- matrix(ncol = 4, nrow = ncol(M))
  ADE <- matrix(ncol = 4, nrow = ncol(M))
  PM <- matrix(ncol = 4, nrow = ncol(M))
  TE <- matrix(ncol = 4, nrow = ncol(M))
  
  # from linear models
  xm <- matrix(ncol = 4, nrow = ncol(M))
  my <- matrix(ncol = 4, nrow = ncol(M))
  
  for (i in 1:ncol(M)) {
    
    dat.x <- data.frame(X = X, Mi = M[, i], covar = cbind(covar, U))
    dat.y <- data.frame(X = X, Mi = M[, i], covar = cbind(covar, U), Y = Y)
    
    # ici cas de deux reg lineaires pour les deux associations
    # TODO refaire pour les autres regressions
    
    
    mod1 <- stats::lm(Mi ~ X + ., data = dat.x)
    
    if(mod2_type=="continuous"){
    mod2 <- stats::lm(Y ~ X + Mi + ., data = dat.y)
    }
    
    if(mod2_type=="surv_Cox"){
    mod2 = survival::survreg(Y ~ X + Mi , dist='exponential', data = dat.y)
    }
    
    
    # # for linear models
    # xm[i, ] <- summary(mod1)$coeff[2, ] # effect of X
    # my[i, ] <- summary(mod2)$coeff[3, ] # effect of M
    # 
    
    
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
  
  return(list(ACME = ACME,
              ADE = ADE,
              PM = PM,
              TE = TE,
              xm = xm,
              my = my))
  
}


plot_summary_ACME <- function(ACME) {
  
  # for check problem
  res_med$ACME = ACME
  p <- ggplot(res_med$ACME, aes(est, stats::reorder(feat, est), color = pval <= 0.05, shape = pval <= 0.05)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = CI_2.5, xmax = CI_97.5)) +
    geom_point(size = 0.8) +
    theme_bw() +
    xlab("ACME (Average Causal Mediation Effect)") +
    ylab("CpG") +
    theme(panel.border = element_blank(),
          panel.spacing = unit(0.01, "lines"),
          axis.ticks = element_blank()) +
    scale_color_manual(values = c("black", "red"))
  
  print(p)
}