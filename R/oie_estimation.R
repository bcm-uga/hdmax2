oie_estimation <- function(X, m, Y, C, boots = 100) {
  
  # bootstrap
  acme_sum <- matrix(nrow = 1, ncol = boots)
  
  for (i in 1:ncol(acme_sum)) {
    samp <- sample(length(X), replace = T)
    
    # effet B m -> Y
    dat.1 <- data.frame(X, m, C)
    mod1 <- lm(Y[samp] ~ ., data = dat.1[samp, ])
    B <- as.data.frame(summary(mod1)$coeff[3:(ncol(m) + 2), ])
    
    # effet A X -> M
    mod2 <- lm(m[samp, ] ~ X[samp] + C[samp, ])
    A <- t(sapply(summary(mod2), function(x) x$coeff[2, ]))
    A <- data.frame(CpG = rownames(A), A)
    # A <- separate(A, CpG, c("0", "CpG"), " ")[, -1]
    
    colnames(B) <- c("B", "B_sd", "B_tv", "B_pv")
    colnames(A)[2:5] <- c("A", "A_sd", "A_tv", "A_pv")
    
    ab <- cbind(A, B)
    rownames(ab) <- NULL
    
    # effet A*B
    ab$AB <- ab$A * ab$B
    
    acme_sum[i] <- sum(ab$AB)
  }
  
  return(oie = as.vector(acme_sum))
}