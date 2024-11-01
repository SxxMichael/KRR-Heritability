#a function to generate population file from the ped and info file
create_population_file_ped <- function(ped, info, p, g_sigma, err_sigma, maf_lb = 0.01, c=1, deg=2, sigma = 1) {
  maf <- info$maf
  #c<-5
  #sample the causal snps from the ped file
  causal_snps_id <- sample(1:ncol(ped), p)
  causal_snps <- ped[, causal_snps_id]
  G_scale <- apply(causal_snps, 2, scale)
  
  
  #try to finish the remaining part
  #simulate genotype function
  #G_genotypes <- ((G_scale %*% rnorm(p,0,g_sigma)))^3  + (G_scale %*% rnorm(p,0,g_sigma))*5
  #G_genotypes <- sin(G_scale%*% rnorm(p,0,g_sigma)) + 5
  G_genotypes <- sin((G_scale %*% rnorm(p,0,g_sigma))) + (G_scale %*% rnorm(p,0,g_sigma)) * 2
  #generate clinical variables
  n <- nrow(ped)
  x1 <- sample(c(0,1), n, replace=T)
  x2 <- runif(n)
  beta<- rep(0,3)
  F_clinical_variables <- cbind(1,x1,x2) %*% beta
  e <- rnorm(n, 0, err_sigma)
  phe <- F_clinical_variables + G_genotypes + e
  herit <- var(G_genotypes)/ var(phe)
  K_lin <- (1/p) * G_scale %*% t(G_scale)
  K_poly <- (c + K_lin )^ deg
  Distance <- as.matrix(dist(G_scale, diag = T, upper = T))
  K_gauss <- exp(-Distance^2/(2*p*sigma^2))
  #lambda<- 10
  var_list <- list(
    herit = herit,
    data = cbind(phe, x1, x2, G_genotypes, G_scale),
    K_lin = K_lin,
    K_poly = K_poly,
    K_gauss = K_gauss
    
  )
  return(var_list)
}