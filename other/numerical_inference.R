library(FragmentCP)

temp_fragment_data7 <- function(mu = 0, r = 5, sigma, n = 100, m = 5, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3){
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

temp_fragment_data9 <- function(mu = 0, r = 5, sigma, n = 100, m = 5, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3){
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


iteration <- 300
r1 = 3
r2 = 3
lambda = 0.001
ext = 0.1
maxIt = 1
B=1000
alpha1 = 0.001
alpha2 = 0.05
CI_matrix1 <- matrix(NaN, nrow = iteration, ncol = 2)
CI_matrix2 <- matrix(NaN, nrow = iteration, ncol = 2)

initial <- rep(0, iteration)
refined <- rep(0, iteration)

for (l in 1:iteration) {
  print(c("iteration", l))
  set.seed(1000+100*l)
  data1 = temp_fragment_data10(mu = 0, r = 3, sigma=1, n = 100, m = 30, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.6)
  data2 = temp_fragment_data13(mu = 0, r = 3, sigma=1, n = 100, m = 30, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.6)
  data = list("t"= rbind(data1$t, data2$t), "y" = rbind(data1$y, data2$y), "r" = cbind(data1$r, data2$r))
  xi_set = c(450, 500)
  CV_cpt_result = CV_search_DP_fragment(data$t, data$y, data$r, r = 3, lambda, xi_set, ext, maxIt)
  min_idx = which.min(CV_cpt_result$test_error) 
  xi_set[min_idx]
  cpt_init = unlist(CV_cpt_result$cpt_hat[min_idx])
  print(cpt_init)
  
  if(length(cpt_init) != 1){
    next
  }
  else{
    #local refinement
    initial[l] <- cpt_init
    refined_eta <- local_refine_fragment(data$t, data$y, data$r, cpt_init,
                                         r=3, lambda, ext = 0.1, maxIt = 1,
                                         w = 0.75)
    refined[l] <- refined_eta
    print(c("refined", refined_eta))
    # estimate jump size
    kappa2_list <-  kappa2_fragment(data$t, data$y, data$r, cpt_init, r=3, lambda, ext = 0.1, maxIt = 1)
    kappa2 <- kappa2_list$kappa2
    print(c("kappa", kappa2))
    C <- kappa2_list$C
    sigma2_list <- sigma2_fragment(cpt_init, kappa2, C, data$t, data$y, data$r)
    sigma_before <- sigma2_list[[1]]
    sigma_after <- sigma2_list[[2]]
    drift <- omega_fragment(cpt_init, kappa2, C, data$t, data$y, data$r)
    
    d <- rep(0,B)
    for (b in 1:B) {
      set.seed(100+10*b+10*l)
      d[b] <- seq(-280, 280)[which.min(simu.2BM_Drift1(280, drift, sigma_before, sigma_after))]
    }
    print(quantile(d, probs = c(alpha1/2, 1-alpha1/2))/kappa2)
    print(quantile(d, probs = c(alpha2/2, 1-alpha2/2))/kappa2)
    
    CI_matrix1[l,] = quantile(d, probs = c(alpha1/2, 1-alpha1/2))/kappa2 + refined_eta
    CI_matrix2[l,] = quantile(d, probs = c(alpha2/2, 1-alpha2/2))/kappa2 + refined_eta
  }
}
CI_matrix1 <- na.omit(CI_matrix1)
CI_matrix2 <- na.omit(CI_matrix2)

width1 <- rep(0, nrow(CI_matrix))
width2 <- rep(0, nrow(CI_matrix))

true_change <- 100

count1 <- 0
count2 <- 0


for(i in 1:nrow(CI_matrix1)){
  if (true_change >= CI_matrix1[i,1] & true_change <= CI_matrix1[i,2]){
    count1 <- count1 +1
  }
  width1[i] <- CI_matrix1[i,2] - CI_matrix1[i,1]
}

for(i in 1:nrow(CI_matrix2)){
  if (true_change >= CI_matrix2[i,1] & true_change <= CI_matrix2[i,2]){
    count2 <- count2 +1
  }
  width2[i] <- CI_matrix2[i,2] - CI_matrix2[i,1]
}


initial <- na.omit(initial)
refined <- na.omit(refined)
truth <- rep(true_change, length(initial))

print(nrow(CI_matrix))
dist_initial.mean <- mean(abs(initial - truth)/(2*true_change))
dist_initial.sd <- sd(abs(initial - truth)/(2*true_change))

dist_refined.mean <- mean(abs(refined - truth)/(2*true_change))
dist_refined.sd <- sd(abs(refined - truth)/(2*true_change))

coverage1 <- count1/nrow(CI_matrix1)
coverage2 <- count2/nrow(CI_matrix2)

width1.mean <- mean(width1)
sd.wideth1.mean <- sd(width1)

width2.mean <- mean(width2)
sd.wideth2.mean <- sd(width2)

result <- data.frame(dist_initial.mean, dist_initial.sd, dist_refined.mean, dist_refined.sd,
                     coverage1,coverage2, width1.mean,sd.wideth1.mean,width2.mean, sd.wideth2.mean)

print(result)
  
  
  
  
