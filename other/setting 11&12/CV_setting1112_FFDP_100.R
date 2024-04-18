library(FragmentCP)
library(changepoints)

cov_function1 <- function(t){
  temp <- sqrt(5)*(6*t^2-6*t+1)
  return(temp)
}

cov_function2 <- function(t){
  temp <- sqrt(2)*log(t+0.1)
  return(temp)
}

cov_function3 <- function(t){
  temp <- (t-0.5)^2/0.2
}

evaluate_cov_function1 <- function(t){
  evaluate <- matrix(0, nrow = length(t), ncol = 2)
  evaluate[, 1] <- cov_function1(t)
  evaluate[, 2] <- cov_function2(t)
  return(evaluate)
}

evaluate_cov_function2 <-  function(t){
  evaluate <- matrix(0, nrow = length(t), ncol = 2)
  evaluate[, 1] <- cov_function1(t)
  evaluate[, 2] <- cov_function3(t)
  return(evaluate)
}

temp_fragment_data11 <- function(mu = 0, sigma, n = 100, m = 5, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3){
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
    C_mat = matrix(0, nrow = 2, ncol = 2)
    diag(C_mat) = c(1,0.8)
    temp = evaluate_cov_function1(tobs)
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

temp_fragment_data12 <- function(mu = 0, sigma, n = 100, m = 5, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3){
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
    C_mat = matrix(0, nrow = 2, ncol = 2)
    diag(C_mat) = c(1,1.5)
    temp = evaluate_cov_function2(tobs)
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
sigma1 = 1
sigma2 = 5

iteration <- 100
true_K <- 1
true_change <- 100
true_K_mat <- matrix(true_K, nrow = 1, ncol = iteration)


K_matrix <- matrix(0, nrow = 1, ncol = iteration)
Hausdroff_distance <- matrix(0, nrow = 1, ncol = iteration)


for (i in 1:iteration) {
  set.seed(1000+10*i)
  print(c("iteration",i))
  data1 = temp_fragment_data11(mu = 0, sigma1, n = 100, m = 40, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.6)
  data2 = temp_fragment_data12(mu = 0, sigma2, n = 100, m = 40, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.6)
  data = list("t"= rbind(data1$t, data2$t), "y" = rbind(data1$y, data2$y), "r" = cbind(data1$r, data2$r))
  # m = 10
  # xi_set = c(2200, 2300, 2500, 2700, 2900, 3000)
  
  # m =20
  # xi_set = c(3500, 4000, 4300, 4500, 4900)
  
  # m =30
  # xi_set = c(4300, 4500, 4900, 5500, 5900)
  
  # m = 40
  xi_set = c(4300, 4500, 4900, 5500, 5900)
  
  CV_cpt_result = CV_search_DP_fragment(data$t, data$y, data$r, r = 3, lambda, xi_set, ext, maxIt)
  min_idx = which.min(CV_cpt_result$test_error) 
  xi_set[min_idx]
  cpt_init = unlist(CV_cpt_result$cpt_hat[min_idx])
  if(length(cpt_init) == 0) {
    cpt_init <- 0
  }
  print(cpt_init)
  K_init = unlist(CV_cpt_result$K_hat[min_idx])
  Hausdroff_distance[1,i] <-  Hausdorff.dist(cpt_init, true_change)
  K_matrix[1,i] <- K_init
}

print(K_matrix)

K_proportion_large <- sum(K_matrix > true_K_mat)/iteration
K_proportion_small <- sum(K_matrix < true_K_mat)/iteration
K_proportion_equal <- sum(K_matrix == true_K_mat)/iteration
mean.diff_K <- mean(abs(K_matrix - true_K_mat))
var.diff_K <- var(abs(K_matrix - true_K_mat)[1,])

mean.Hausdroff_distance <- mean(Hausdroff_distance)
var.Hausdroff.distance <- var(Hausdroff_distance[1,])

result = list(cpt_dist_mean = mean.Hausdroff_distance,
              cpt_dist_var = var.Hausdroff.distance,
              K_diff_mean = mean.diff_K,
              K_diff_var = var.diff_K,
              K_proportion_large = K_proportion_large,
              K_proportion_small = K_proportion_small,
              K_proportion_equal = K_proportion_equal)
print(result)




