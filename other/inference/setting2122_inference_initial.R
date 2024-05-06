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

B=1000
alpha1 = 0.1
alpha2 = 0.15
# alpha3 = 0.1
# alpha4 = 0.2

n1 <- 800
M1 <- 370

# n2 <- 650
# M2 <- 300

iteration <- 200
CI_matrix1 <- matrix(NaN, nrow = iteration, ncol = 2)
CI_matrix2 <- matrix(NaN, nrow = iteration, ncol = 2)
# CI_matrix3 <- matrix(NaN, nrow = iteration, ncol = 2)
# CI_matrix4 <- matrix(NaN, nrow = iteration, ncol = 2)

true_change <- 150

initial <- rep(NaN, iteration)
refined <- rep(NaN, iteration)

# C1 = matrix(0, nrow = 3, ncol = 3)
# diag(C1) = c(0.8, seq(from = 2, by = -0.6, length.out = 3-1))
# 
# C2 = matrix(0, nrow = 3, ncol = 3)
# diag(C2) = c(3.5, seq(from = 2.5, by = -0.6, length.out = 3-1))

for (l in 1:iteration) {
  print(c("iteration", l))
  set.seed(1000+10*l)
  # set.seed(90)
  data1 = temp_fragment_data21(mu = 0, r = 3, sigma=1, n = 150, m = 20, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.5)
  data2 = temp_fragment_data22(mu = 0, r = 3, sigma=1, n = 150, m = 20, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.5)
  data = list("t"= rbind(data1$t, data2$t), "y" = rbind(data1$y, data2$y), "r" = cbind(data1$r, data2$r))
  xi_set = c(760, 860)
  lambda_set = c(0.0003, 0.0005, 0.001)
  
  # lambda_set = c(0.0005, 0.001, 0.0015)
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
  selected_lambda <- lambda_set[ceiling(min_idx/length(xi_set))]
  
  cpt_init = unlist(CV_result$cpt_hat[min_idx])
  print(cpt_init)
  
  if(length(cpt_init) != 1){next}
  else {
    refined_eta <- local_refine_fragment(data$t, data$y, data$r, cpt_init,
                                         r=3, selected_lambda, ext = 0.1, maxIt = 1,
                                         w = 0.75)
    
    refined_eta_with_truth <- local_refine_fragment_truth(data$t, data$y, data$r, cpt_init, C1, C2,
                                                          r=3, slected_lambda, ext = 0.1, maxIt = 1,
                                                          w = 0.75)
    print(c("refined", refined_eta, "refined with truth", refined_eta_with_truth))
    
    initial[l] <- cpt_init
    refined[l] <- refined_eta
    
    # estimate jump size
    kappa2_list <-  kappa2_fragment(data$t, data$y, data$r, cpt_init, r=3, selected_lambda, ext = 0.1, maxIt = 1)
    kappa2 <- kappa2_list$kappa2
    print(c("kappa", kappa2))
    C <- kappa2_list$C
    sigma2_list <- sigma2_fragment(cpt_init, kappa2, C, data$t, data$y, data$r)
    sigma_before <- sigma2_list[[1]]
    sigma_after <- sigma2_list[[2]]
    drift <- omega_fragment(cpt_init, kappa2, C, data$t, data$y, data$r)
    
    d1 <- rep(0,B)
    for (b in 1:B) {
      set.seed(100+10*b+10*l)
      d1[b] <- seq(-n1*M1, n1*M1)[which.min(simu.2BM_Drift(n1, M1, drift, sigma_before, sigma_after))]/n1
    }
    
    
    # d2 <- rep(0,B)
    # for (b in 1:B) {
    #   set.seed(100+10*b+10*l)
    #   d2[b] <- seq(-n2*M2, n2*M2)[which.min(simu.2BM_Drift(n2, M2, drift, sigma_before, sigma_after))]/n2
    # }
    # 
    
    print(quantile(d1, probs = c(alpha1/2, 1-alpha1/2))/kappa2)
    print(quantile(d1, probs = c(alpha2/2, 1-alpha2/2))/kappa2)
    # print(quantile(d2, probs = c(alpha3/2, 1-alpha3/2))/kappa2)
    # print(quantile(d2, probs = c(alpha4/2, 1-alpha4/2))/kappa2)
    
    
    CI_matrix1[l,1] = floor(quantile(d1, probs = alpha1/2)/kappa2) + refined_eta
    CI_matrix1[l,2] = ceiling(quantile(d1, probs = 1-alpha1/2)/kappa2) + refined_eta
    
    CI_matrix2[l,1] = floor(quantile(d1, probs = alpha2/2)/kappa2) + refined_eta
    CI_matrix2[l,2] = ceiling(quantile(d1, probs = 1-alpha2/2)/kappa2) + refined_eta
    
    # CI_matrix3[l,1] = floor(quantile(d2, probs = alpha3/2)/kappa2) + refined_eta
    # CI_matrix3[l,2] = ceiling(quantile(d2, probs = 1-alpha3/2)/kappa2) + refined_eta
    # 
    # CI_matrix4[l,1] = floor(quantile(d2, probs = alpha4/2)/kappa2) + refined_eta
    # CI_matrix4[l,2] = ceiling(quantile(d2, probs = 1-alpha4/2)/kappa2) + refined_eta
    
    # CI_matrix1[l,] = ceiling(quantile(d, probs = c(alpha1/2, 1-alpha1/2))/kappa2 + refined_eta)
    # CI_matrix2[l,] = ceiling(quantile(d, probs = c(alpha2/2, 1-alpha2/2))/kappa2 + refined_eta)
    # CI_matrix3[l,] = ceiling(quantile(d, probs = c(alpha3/2, 1-alpha3/2))/kappa2 + refined_eta)
    # CI_matrix4[l,] = ceiling(quantile(d, probs = c(alpha4/2, 1-alpha4/2))/kappa2 + refined_eta)
  }
}

CI_matrix1 <- na.omit(CI_matrix1)
CI_matrix2 <- na.omit(CI_matrix2)
# CI_matrix3 <- na.omit(CI_matrix3)
# CI_matrix4 <- na.omit(CI_matrix4)

width1 <- rep(0, nrow(CI_matrix1))
width2 <- rep(0, nrow(CI_matrix2))
# width3 <- rep(0, nrow(CI_matrix3))
# width4 <- rep(0, nrow(CI_matrix4))



count1 <- 0
count2 <- 0
# count3 <- 0
# count4 <- 0


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

# for(i in 1:nrow(CI_matrix3)){
#   if (true_change >= CI_matrix3[i,1] & true_change <= CI_matrix3[i,2]){
#     count3 <- count3 +1
#   }
#   width3[i] <- CI_matrix3[i,2] - CI_matrix3[i,1]
# }
# 
# for(i in 1:nrow(CI_matrix4)){
#   if (true_change >= CI_matrix4[i,1] & true_change <= CI_matrix4[i,2]){
#     count4 <- count4 +1
#   }
#   width4[i] <- CI_matrix4[i,2] - CI_matrix4[i,1]
# }


initial <- na.omit(initial)
refined <- na.omit(refined)
truth <- rep(true_change, length(initial))

print(nrow(CI_matrix1))
abs_dist_initial.mean <- mean(abs(initial - truth)/(2*true_change))
abs_dist_initial.sd <- sd(abs(initial - truth)/(2*true_change))

abs_dist_refined.mean <- mean(abs(refined - truth)/(2*true_change))
abs_dist_refined.sd <- sd(abs(refined - truth)/(2*true_change))

prop_large_init <- sum(initial > true_change)/(length(initial))
prop_small_init <- sum(initial < true_change)/(length(initial))

prop_large_refined <- sum(refined > true_change)/(length(initial))
prop_small_refined <- sum(refined < true_change)/(length(initial))

dist_initial.mean <- mean((initial - truth)/(2*true_change))
dist_initial.sd <- sd((initial - truth)/(2*true_change))

diff_refined.mean <- mean((refined - truth)/(2*true_change))
diff_refined.sd <- sd((refined - truth)/(2*true_change))




coverage1 <- count1/nrow(CI_matrix1)
coverage2 <- count2/nrow(CI_matrix2)
# coverage3 <- count3/nrow(CI_matrix3)
# coverage4 <- count4/nrow(CI_matrix4)

width1.mean <- mean(width1)
sd.width1.mean <- sd(width1)
width1.med <- median(width1)

width2.mean <- mean(width2)
sd.width2.mean <- sd(width2)
width2.med <- median(width2)

# width3.mean <- mean(width3)
# sd.width3.mean <- sd(width3)
# width3.med <- median(width3)
# 
# width4.mean <- mean(width4)
# sd.width4.mean <- sd(width4)
# width4.med <- median(width4)

local_result <- data.frame(abs_dist_initial.mean, abs_dist_initial.sd, abs_dist_refined.mean, abs_dist_refined.sd, prop_large_init,
                           prop_small_init, prop_large_refined, prop_small_refined, dist_initial.mean, dist_initial.sd, diff_refined.mean,diff_refined.sd)


# inf_result <- data.frame(coverage1,coverage2,coverage3,coverage4, width1.mean,sd.width1.mean,width1.med, width2.mean, sd.width2.mean,
#                      width2.med, width3.mean, sd.width3.mean, width3.med, width4.mean, sd.width4.mean, width4.med)

inf_result <- data.frame(coverage1,coverage2, width1.mean,sd.width1.mean,width1.med, width2.mean, sd.width2.mean,
                         width2.med)

print(local_result)
print(inf_result)

save.image(file = "n150_m20_3.RData")






