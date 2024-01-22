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
    TE[i, ] <- c(med$tau.coef, med$tau0.ci[1], med$tau0.ci[2], med$tau.p)
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