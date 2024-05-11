library(FragmentCP)
library(changepoints)

cov1 <- function(mu = 0, r = 5, sigma, n = 100, m = 5, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3){
  Ly <- matrix(0, n, m)
  L <- domain[2] - domain[1]
  a_vec <- runif(n, min = domain[1], max = domain[2] - delta * L)
  Lt <- t(sapply(1:n, function(i) {
    a <- a_vec[i]
    a + sort(runif(m, min = 0, max = delta * L))
  }))
  Ly0 <- t(apply(Lt, 1, function(tobs) {
    m = length(tobs)
    if(is.function(mu)) 
      mui <- mu(tobs)
    else if(length(mu) == 1) 
      mui <- rep(mu, length(tobs))
    else stop("mu must be a scalar or a function.")
    C_mat = matrix(0, nrow = r, ncol = r)
    diag(C_mat) = seq(from = 2, by = -0.6, length.out = r)
    # print(C_mat)
    temp = evaluate_basis(r, c(0,1), tobs)
    Sigma_mat = temp %*% C_mat %*% t(temp)
    y0 <- mui + MASS::mvrnorm(n = 1, mu = rep(0, m), Sigma = Sigma_mat*sigma^2)
    return(c(y0))
  }))
  Ly <- Ly0 + matrix(rnorm(n = n*m, sd = sigma_epsilon), n, m)
  Lr <- NULL
  for(i in 1:n){
    Lr = cbind(Lr, tcrossprod(Ly[i,]))
  }
  R <- list(t = Lt, y = Ly, r = Lr)
  return(R)
}

cov2 <- function(mu = 0, r = 5, sigma, n = 100, m = 5, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3){
  Ly <- matrix(0, n, m)
  L <- domain[2] - domain[1]
  a_vec <- runif(n, min = domain[1], max = domain[2] - delta * L)
  Lt <- t(sapply(1:n, function(i) {
    a <- a_vec[i]
    a + sort(runif(m, min = 0, max = delta * L))
  }))
  Ly0 <- t(apply(Lt, 1, function(tobs) {
    m = length(tobs)
    if(is.function(mu)) 
      mui <- mu(tobs)
    else if(length(mu) == 1) 
      mui <- rep(mu, length(tobs))
    else stop("mu must be a scalar or a function.")
    C_mat = matrix(0, nrow = r, ncol = r)
    diag(C_mat) = c(seq(from = 5, by = -0.8, length.out = r))
    # print(C_mat)
    temp = evaluate_basis(r, c(0,1), tobs)
    Sigma_mat = temp %*% C_mat %*% t(temp)
    y0 <- mui + MASS::mvrnorm(n = 1, mu = rep(0, m), Sigma = Sigma_mat*sigma^2)
    return(c(y0))
  }))
  Ly <- Ly0 + matrix(rnorm(n = n*m, sd = sigma_epsilon), n, m)
  Lr <- NULL
  for(i in 1:n){
    Lr = cbind(Lr, tcrossprod(Ly[i,]))
  }
  R <- list(t = Lt, y = Ly, r = Lr)
  return(R)
}

lambda = 0.00001
ext = 0.1
maxIt = 1
r1 = 3
r2 = 3
sigma1 = 1
sigma2 = 1
alpha1 = 0.1
B=1000
n1 = 650
M1 = 250

# xi_set = c(650, 700, 800, 900, 1000)
# lambda_set = c(0.0003, 0.0005, 0.001)

xi_set = c(650, 700)
lambda_set = c(0.0003, 0.0005)

# generate data
data1 = cov1(mu = 0, r = r1, sigma1, n = 100, m = 30, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.6)
data2 = cov2(mu = 0, r = r2, sigma2, n = 100, m = 30, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.6)
data = list("t"= rbind(data1$t, data2$t), "y" = rbind(data1$y, data2$y), "r" = cbind(data1$r, data2$r))

# FFDP
CV_result <- CV_search_DP_fragment_lambda_xi(data$t, data$y, data$r, r = 3, lambda_set, xi_set, ext, maxIt)
min_idx = which.min(CV_result$test_error) 

temp <- min_idx %% length(xi_set)
if(temp == 0) {temp <-  length(xi_set)}
selected_xi <- xi_set[temp]
selected_lambda <- lambda_set[ceiling(min_idx/length(xi_set))]
selected_lambda <- lambda_set[ceiling(min_idx/length(xi_set))]

cpt_init = unlist(CV_result$cpt_hat[min_idx])
print(cpt_init)


#local refinement
refined_eta <- local_refine_fragment(data$t, data$y, data$r, cpt_init,
                                     r=3, selected_lambda, ext = 0.1, maxIt = 1,
                                     w = 0.75)

# confidence interval construction
kappa2_list <-  kappa2_fragment(data$t, data$y, data$r, cpt_init, r=3, selected_lambda, ext = 0.1, maxIt = 1)
kappa2 <- kappa2_list$kappa2
C <- kappa2_list$C
sigma2_list <- sigma2_fragment(cpt_init, kappa2, C, data$t, data$y, data$r)
sigma_before <- sigma2_list[[1]]
sigma_after <- sigma2_list[[2]]
sigma2 <- mean(c(sigma_before, sigma_after))
drift <- omega_fragment(cpt_init, kappa2, C, data$t, data$y, data$r)

d1 <- rep(0,B)
for (b in 1:B) {
  set.seed(100+10*b)
  d1[b] <- seq(-n1*M1, n1*M1)[which.min(simu.2BM_Drift_1LRV(n1, M1, drift, sigma2))]/n1
}

d2 <- rep(0,B)
for (b in 1:B) {
  set.seed(100+10*b)
  d2[b] <- seq(-n1*M1, n1*M1)[which.min(simu.2BM_Drift(n1, M1, drift, sigma_before, sigma_after))]/n1
}

# CI by step 4
CI_matrix_lower =  refined_eta - ceiling(quantile(d1, probs = 1-alpha1/2)/kappa2) 
CI_matrix_upper =  refined_eta - floor(quantile(d1, probs = alpha1/2)/kappa2)

# CI by step 4'
CI_matrix_lower =  refined_eta - ceiling(quantile(d2, probs = 1-alpha1/2)/kappa2) 
CI_matrix_upper =  refined_eta - floor(quantile(d2, probs = alpha1/2)/kappa2)















