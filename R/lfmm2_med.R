##' lfmm2 adapted function for hdmax2 from LEA::lfmm2
##' @param input a response variable matrix with n rows and p columns
##' @param env An explanatory variable matrix with n rows and d columns.
##' @param K latent factor number
##' @param lambda ridge penalization parameter
##' @param effect.sizes true or false to obtain effect sizes
##' @return an object with the following attributes 
##' @export 
##' @author Florence Pittion
##' @examples
##' data(simu_data)
##' mod.lfmm1 = lfmm2_med(input = simu_data$M, 
##' env = simu_data$X_binary, 
##' K = 5,
##' effect.sizes = FALSE)


lfmm2_med = function(input,
                  env, 
                  K, 
                  lambda = 1e-5,
                  effect.sizes = FALSE){
  
  ## Check input response matrix 
  ## LEA  
  if (is.character(input)){
    Y <- read_lfmm(input)
    lst.unique <- unique(as.numeric(Y))
    if (9 %in% lst.unique){
      stop("'input' file contains missing data (9's). Use the 'impute()' function to impute them.")
    }
    if (-9 %in% lst.unique){
      stop("'input' file contains missing data (-9's). Use the 'impute()' function to impute them.")
    }
  } else {
    ## Y is an R object       
    if (is.null(input)){
      stop("NULL value for argument 'input'.")
    }
    Y <- as.matrix(input)
    Y[Y == 9] <- NA
    Y[Y == -9] <- NA
    if (anyNA(Y)) {
      stop("The input matrix contains missing values: NA, 9 or -9 not allowed.")
    }
  }
  
  ## Check independent/covariate env matrix  
  ## LEA 
  if (is.character(env)){
    X <- read_env(env)
    if (anyNA(X)){
      stop("'env' file contains missing data (NA).")
    }
  } else {
    if (is.null(env)){
      stop("NULL value for argument 'env'.")
    }
    X <- as.matrix(env)
    if (anyNA(X)) {
      stop("The environmental matrix contains NA.")
    }
  }
  
  if (length(K) > 1){
    stop("Multiple values of K not allowed.")
  }
  if (lambda <= 0){
    stop("The ridge regularization parameter must be positive.")
  }
  
  d <-  ncol(X) #number of environmental variables
  n <-  nrow(X) #number of individuals
  
  if (nrow(Y) != n){
    stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
  }
  
  if (n < d) {
    stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
  }
  
  # centering  
  Xs <- scale(X, scale = FALSE)
  Ys <- scale(Y, scale = FALSE)
  
  # run SVD of X: X = Q Sigma R
  
  svx <- svd(x = Xs, nu = n)
  Q <- svx$u
  
  d_lambda <- c(sqrt(lambda/(lambda + svx$d)), rep(1, n-d))
  d_lambda_inv <- c(sqrt((lambda + svx$d)/lambda), rep(1, n-d))
  D_inv <- diag(d_lambda_inv)
  D  <- diag(d_lambda)
  
  # run SVD of modified Y    
  svk <- svd(D %*% t(Q) %*% Ys, nu = K)
  
  if (K > 1) {
    Sigma_k <- diag(svk$d[1:K])
  } else {
    Sigma_k <- as.matrix(svk$d[1])
  }
  
  # compute the latent matrix W
  W <- Q %*% D_inv %*% tcrossprod(svk$u %*% Sigma_k, svk$v[,1:K])
  
  # compute LFMM factors U and loadings V
  # Non orthogonal factors
  U <- crossprod(t(Q %*% D_inv), svk$u %*% Sigma_k)
  #U <- Q %*% D_inv %*% svk$u %*% Sigma_k
  V <- svk$v[,1:K]
  
  # compute environmental effect sizes 
  if (effect.sizes){
    B <- (t(Ys - W) %*% Xs) %*% solve(t(Xs) %*% Xs + diag(lambda, nrow = d, ncol = d))
    B <- as.matrix(B)
  } else
  {B <-  matrix(NA)}
  
  
  
  obj= list()
  # obj$K <- as.integer(K)
  # obj$lambda <- as.numeric(lambda)
  # obj$B <- as.matrix(B)
  obj$U <- as.matrix(U)
  obj$V <- as.matrix(V)
  
  ## LEA 
  class(obj) = "lfmm2"
  return(obj)
}