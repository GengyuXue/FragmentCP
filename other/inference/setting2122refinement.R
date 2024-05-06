library(FragmentCP)

local_refine_fragment_truth = function(Lt, Ly, Lr, cpt_init_hat, C1, C2, r, lambda, ext, maxIt, w = 0.9){
  if(length(cpt_init_hat) == 0){return(cpt_init_hat)}
  else{
    n = nrow(Lt)
    m = ncol(Lt)
    Khat = length(cpt_init_hat)
    cpt_hat_long = c(0, cpt_init_hat, n)
    cpt_refined = rep(0, Khat+1)
    C_hat_list <- vector("list", Khat+1)
    C_hat_list[[1]] <- C1
    C_hat_list[[2]] <- C2
    
    for (k in 1:Khat){
      s = ceiling(w*cpt_hat_long[k] + (1-w)*cpt_hat_long[k+1])
      e = floor((1-w)*cpt_hat_long[k+1] + w*cpt_hat_long[k+2])
      lower = s + 2
      upper = e - 2
      
      b = sapply(lower:upper, function(eta) (error_test_fragment(Lt, Lr, r, s, eta,C_hat_list[[k]]) 
                                             + error_test_fragment(Lt, Lr, r, eta+1, e, C_hat_list[[k+1]])))
      cpt_refined[k+1] = s + which.min(b)
    }
    return(cpt_refined[-1])
  }
}

CV_search_DP_fragment_lambda_xi = function(Lt, Ly, Lr, r, lambda_set, xi_set, ext, maxIt){
  change_point_list <- vector("list", length(lambda_set)+ length(xi_set))
  test_error <- rep(0, length(lambda_set)*length(xi_set))
  index <- 1
  for (i in 1:length(lambda_set)) {
    for (j in 1:length(xi_set)) {
      output_temp <- CV_fragment(Lt, Ly, Lr, r, lambda_set[i], xi_set[j], ext, maxIt)
      test_error[index] <- output_temp$test_error
      change_point_list[[index]] <- output_temp$cpt_hat
      index <- index + 1
    }
  }
  return(list(cpt_hat = change_point_list, test_error = test_error))
}


temp_fragment_data21 <- function(mu = 0, r = 5, sigma, n = 100, m = 5, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3){
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
    diag(C_mat) = c(0.5, seq(from = 2, by = -0.6, length.out = r-1))
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

temp_fragment_data22 <- function(mu = 0, r = 5, sigma, n = 100, m = 5, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3){
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
    diag(C_mat) = c(3.5, seq(from = 4, by = -0.5, length.out = r-1))
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


C1 = matrix(0, nrow = 3, ncol = 3)
diag(C1) = c(0.5, seq(from = 2, by = -0.6, length.out = 3-1))

C2 = matrix(0, nrow = 3, ncol = 3)
diag(C2) = c(3.5, seq(from = 4, by = -0.5, length.out = 3-1))
############################
ext = 0.1
maxIt = 1


for (l in 1:10) {
  print(c("iteration", l))
  set.seed(1000+5*l)
  # set.seed(90)
  data1 = temp_fragment_data21(mu = 0, r = 3, sigma=1, n = 300, m = 20, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.5)
  data2 = temp_fragment_data22(mu = 0, r = 3, sigma=1, n = 300, m = 20, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.5)
  data = list("t"= rbind(data1$t, data2$t), "y" = rbind(data1$y, data2$y), "r" = cbind(data1$r, data2$r))
  # xi_set = c(650, 690)
  xi_set = c(1650, 1750)
  # lambda_set = c(0.0005, 0.001)
  lambda_set = c(0.0003, 0.0005, 0.001)
  
  # CV_cpt_result = CV_search_DP_fragment(data$t, data$y, data$r, r = 3, lambda, xi_set, ext, maxIt)
  # min_idx = which.min(CV_cpt_result$test_error) 
  # xi_set[min_idx]
  # cpt_init = unlist(CV_cpt_result$cpt_hat[min_idx])
  # print(cpt_init)
  
  CV_result <- CV_search_DP_fragment_lambda_xi(data$t, data$y, data$r, r = 3, lambda_set, xi_set, ext, maxIt)
  min_idx = which.min(CV_result$test_error) 
  temp <- min_idx %% length(xi_set)
  if(temp == 0) {temp <-  length(xi_set)}
  selected_xi <- xi_set[temp]
  slected_lambda <- lambda_set[ceiling(min_idx/length(xi_set))]
  
  cpt_init = unlist(CV_result$cpt_hat[min_idx])
  print(cpt_init)
  
  if(length(cpt_init) != 1){next}
  else {
    refined_eta <- local_refine_fragment(data$t, data$y, data$r, cpt_init,
                                         r=3, slected_lambda, ext = 0.1, maxIt = 1,
                                         w = 0.75)
    refined_eta_with_truth <- local_refine_fragment_truth(data$t, data$y, data$r, cpt_init, C1, C2,
                                                          r=3, slected_lambda, ext = 0.1, maxIt = 1,
                                                          w = 0.75)
    print(c("refined", refined_eta, "refined with truth", refined_eta_with_truth))
  }
}



# 
# 
# data1 = temp_fragment_data10(mu = 0, r = 3, sigma=1, n = 100, m = 10, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.5)
# data2 = temp_fragment_data13(mu = 0, r = 3, sigma=1, n = 100, m = 10, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.5)
# data = list("t"= rbind(data1$t, data2$t), "y" = rbind(data1$y, data2$y), "r" = cbind(data1$r, data2$r))
# xi_set = 
# lambda_set = c(0.0007, 0.001, 0.002, 0.003)
# 
# CV_result <- CV_search_DP_fragment_lambda_xi(data$t, data$y, data$r, r = 3, lambda_set, xi_set, ext, maxIt)
# min_idx = which.min(CV_result$test_error) 
# temp <- min_idx %% length(xi_set)
# if(temp == 0) {temp <-  length(xi_set)}
# selected_xi <- xi_set[temp]
# slected_lambda <- lambda_set[ceiling(min_idx/length(xi_set))]
#   
# cpt_init = unlist(CV_result$cpt_hat[min_idx])





