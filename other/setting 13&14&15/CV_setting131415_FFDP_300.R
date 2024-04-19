library(FragmentCP)
library(changepoints)

temp_fragment_data13 <- function(mu = 0, r = 5, sigma, n = 100, m = 5, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3){
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
    diag(C_mat) = seq(from = 1.8, by = -0.5, length.out = r)
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

temp_fragment_data14 <- function(mu = 0, r = 5, sigma, n = 100, m = 5, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3){
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
    diag(C_mat) = seq(from = 8, by = -1.5, length.out = r)
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

temp_fragment_data15 <- function(mu = 0, r = 5, sigma, n = 100, m = 5, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3){
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
    diag(C_mat) = seq(from = 1.8, by = -0.2, length.out = r)
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
sigma1 = 1
sigma2 = 1

iteration <- 100
true_K <- 2
true_change <- c(0,75,187,300)
true_K_mat <- matrix(true_K, nrow = 1, ncol = iteration)


K_matrix <- matrix(0, nrow = 1, ncol = iteration)
Hausdroff_distance <- matrix(0, nrow = 1, ncol = iteration)


for (i in 1:iteration) {
  set.seed(1000+10*i)
  print(c("iteration",i))
  data1 = temp_fragment_data13(mu = 0, r = 2, 1, n = 75, m = 15, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.5)
  data2 = temp_fragment_data14(mu = 0, r = 2, 1, n = 112, m = 15, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.5)
  data3 = temp_fragment_data15(mu = 0, r = 2, 1, n = 113, m = 15, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.5)
  
  data = list("t"= rbind(data1$t, data2$t, data3$t), "y" = rbind(data1$y, data2$y, data3$y), "r" = cbind(data1$r, data2$r, data3$r))
  
  # xi_set = c(650, 700, 750, 800)
  # xi_set = c(750, 800, 850, 900)
  xi_set = c(900, 1200, 1300, 1400, 1500)
  
  CV_cpt_result = CV_search_DP_fragment(data$t, data$y, data$r, r = 2, lambda, xi_set, ext, maxIt)
  min_idx = which.min(CV_cpt_result$test_error) 
  xi_set[min_idx]
  cpt_init = c(0, unlist(CV_cpt_result$cpt_hat[min_idx]), 300)
  
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
sd.diff_K <- sqrt(var(abs(K_matrix - true_K_mat)[1,]))

mean.Hausdroff_distance <- mean(Hausdroff_distance)/300
sd.Hausdroff.distance <- sqrt(var(Hausdroff_distance[1,]/300))

result = list(cpt_dist_mean = mean.Hausdroff_distance,
              cpt_dist_sd = sd.Hausdroff.distance,
              K_diff_mean = mean.diff_K,
              K_diff_sd = sd.diff_K,
              K_proportion_large = K_proportion_large,
              K_proportion_small = K_proportion_small,
              K_proportion_equal = K_proportion_equal)
print(result)

