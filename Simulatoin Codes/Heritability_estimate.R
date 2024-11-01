# y is the vector of phenotypes from the sample
# X is the matrix of clinical variables with n rows and q columns 
# (n is the number of sample and q is the number of clinical variables)
# K is the kernel matrix generated from the genotypes
# lambda is a positive tuning parameter that controls the performance (by default it is 0.1)
RKHS_ls <- function(y, X, K, lambda)
{
  n <- length(y) # extract the sample size using the length of the vector y
  Id <- diag(1, n)
  XtX <- t(X) %*% X  # an alternative is XtX <- crossprod(X)
  P <- 0#X %*% solve(XtX) %*% t(X)
  
  #alpha_hat <- solve(K%*%(Id - P)%*%K + lambda * K) %*% K%*%(Id - P)%*%y
  alpha_hat <- solve((Id - P)%*%K + lambda * Id) %*% (Id-P) %*% y
  beta_hat <- 0#solve(XtX) %*% t(X) %*% (y - K %*% alpha_hat)
  
  var_list <- list(
    alpha_hat = alpha_hat,
    beta_hat = beta_hat
  )
  return(var_list)
}


# Function for estimating (broad-sense) heritability
RKHS_herit <- function(y, X, K, lambda)
{
  X<-0
  temp <- RKHS_ls(y, 0, K, lambda)  # call the function RKHS_ls and store the results 
  # into the variable temp
  alpha_hat <- temp$alpha_hat
  beta_hat  <- temp$beta_hat
  
  n  <- length(y)
  q  <- 0#ncol(X)
  Id <- diag(1, n)
  J  <- matrix(1, n, n)
  
  var_g <- 1/(n-1) * t(alpha_hat) %*% K %*% (Id - (1/n) * J) %*% K %*% alpha_hat
  #var_e <- (1/(n-q)) * t(y - X %*% beta_hat - K %*% alpha_hat) %*% 
  #                  (y - X %*% beta_hat - K %*% alpha_hat)
  #var_e <- (1/(n-q)) * crossprod(y - X %*% beta_hat - K %*% alpha_hat)
  var_e <- (1/(n-q)) * crossprod(y  - K %*% alpha_hat)
  
  herit <- var_g / (var_g + var_e)
  
  return(c(herit,var_g,var_e))
  
}




