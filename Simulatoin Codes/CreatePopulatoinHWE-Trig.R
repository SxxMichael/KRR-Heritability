#a function to be used to generate population files based on the HWE
#g_sigma is used to generate the effects for the causal snps (need to choose a suitable value)
#so that heritability is between 20% and 80%
#err_sigma is used to generate the random error
#maf_lb is the "smallest" value for the maf
#c, deg are the parameters for calculating the polynomial kernel
#sigma is the parameter in the gaussian kernel

create_population_file <- function (n, p, g_sigma, err_sigma, maf_lb = 0.01, K=1, deg=2, sigma = 1) {
  #generate genotype source
  #maf_lb = maf lower bound
  G1 <- matrix(0, n, p)
  
  #c<-5
  maf <- runif (p, maf_lb, 0.5)
  for(i in 1:n) {
    for(j in 1:p) {
      G1[i,] <- sample(c(0,1,2), p, replace = T, 
                       prob = c((1-maf[j])^2, 2 * (1-maf[j]) * maf[j], maf[j]^2))
    }
  }
  
  G_scale <- apply(G1, 2, scale)
  #simulate genotype function
  G_genotypes <- sin((G_scale %*% rnorm(p,0,g_sigma)))+ (G_scale %*% rnorm(p,0,g_sigma)) * 2
  x1 <- sample(c(0,1), n, replace = T)
  x2 <- runif(n)
  beta <- rep(0,3)
  F_clinical_variables <- cbind(1, x1, x2) %*% beta
  e <- rnorm(n, 0, err_sigma)
  phe <- F_clinical_variables + G_genotypes + e
  herit <- var(G_genotypes) / var(phe)
  #calculate the linear, polynomial and gaussian kernel
  K_lin <- (1/p) * G_scale %*% t(G_scale)
  c<-1
  K_poly <- ( c + K_lin )^ deg
  Distance <- as.matrix(dist(G_scale, diag = T, upper = T))
  K_gauss <- exp(-Distance^2/(2*p*sigma^2))
  var_list <- list(
    herit = herit,
    data = cbind(phe, x1, x2, G_genotypes, G_scale),
    K_lin = K_lin,
    K_poly = K_poly,
    K_gauss = K_gauss
  )
  #return(herit)
  return(var_list)
}
