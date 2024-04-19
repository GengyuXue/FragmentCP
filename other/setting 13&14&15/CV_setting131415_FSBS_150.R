library(fda)
library(ks)
library(mvtnorm)
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

# parameter used to generate seeded interval
C=3
#dimension of data
d=1
iteration = 100


true_K <- 2
true_change <- c(0, 37, 93, 150)
true_K_mat <- matrix(true_K, nrow = 1, ncol = iteration)
K_matrix <- matrix(0, nrow = 1, ncol = iteration)
Hausdroff_distance <- matrix(0, nrow = 1, ncol = iteration)

for (l in 1:100) {
  print(c("iteration", l))
  set.seed(1000+10*l)
  data1 = temp_fragment_data13(mu = 0, r = 2, 1, n = 37, m = 15, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.5)
  data2 = temp_fragment_data14(mu = 0, r = 2, 1, n = 56, m = 15, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.5)
  data3 = temp_fragment_data15(mu = 0, r = 2, 1, n = 57, m = 15, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.5)
  
  data = list("t"= rbind(data1$t, data2$t, data3$t), "y" = rbind(data1$y, data2$y, data3$y), "r" = cbind(data1$r, data2$r, data3$r))
  
  m = ncol(data$y)
  Tb = nrow(data$y)
  nts<-matrix(NaN,nrow=Tb, ncol=1)
  for(i in 1:Tb){
    nts[i] = m
  }
  N_t=sum(nts)
  # N_t
  
  dat <- matrix(NaN, nrow = N_t, ncol = 3)
  
  labels_Ts=matrix(NaN,nrow = N_t,ncol=1)
  aux1=0
  for (i in 1:Tb) {
    aux1=nts[i]+aux1
    for (j in 1:nts[i]) {
      labels_Ts[aux1-nts[i]+j]=i
    }
  }
  # labels_Ts
  
  t_vector <- matrix(NaN, nrow = N_t, ncol =1)
  y_vector <- matrix(NaN, nrow = N_t, ncol =1)
  for (i in 1:Tb) {
    for (j in 1:m) {
      t_vector[(i-1)*m+j,1] <- c(data$t[i,j])
      y_vector[(i-1)*m+j,1] <- c(data$y[i,j])
    }
  }
  
  # create the data matrix containing x_t_is, y_t_is and labels
  dat[, 1] <- labels_Ts
  dat[, 2] <- t_vector
  dat[, 3] <- y_vector
  # dat
  
  #The choice of bandwidth h_bar
  H_bar_stimated=hpi(dat[,2])
  # H_bar_stimated
  h_bar=H_bar_stimated
  
  ######################################################
  ######################### CV #########################
  ######################################################
  #Creating testing and training sets.
  dat_train=dat[(m+1):(2*m),]
  dat_test=dat[1:m,]
  for(i in 1:(Tb/2-1)){
    dat_test=rbind(dat_test,dat[(2*i*m+1):((2*i+1)*m),])
    dat_train=rbind(dat_train,dat[((2*i+1)*m+1):((2*(i+1))*m),])
  }
  dat_train[,1]=ceiling(dat_train[,1]/2)
  
  #We create the seeded intervals to work with training data
  s.inter_train<-seeded.intervals(Tb/2, C)
  
  
  
  #We create the estimator functions 
  g_hat_i=function(eta1,eta2,h,x,dat_train,h_bar,phat)
  {
    res=0
    for (i in (eta1+1):eta2) {
      res=res+statistic(x,dat_train,i,h,h_bar,d,phat)
    }
    res=res/(eta2-eta1)
    return(res)
  }
  
  #We create the possible values for h and tau
  h_int= 0.25
  tau_int=c(9.5, 10)
  l_tau_int=length(tau_int)
  
  #We compute errors of estimation
  
  errors=matrix(NaN,1,l_tau_int)
  K_cv <- matrix(NaN,1,l_tau_int)
  CP_cv <- list()
  #start_time <- Sys.time()
  for (ind in 1:l_tau_int) {
    S=seedBS(dat_train, s.inter_train, h_int, h_bar, tau_int[ind], d, m)
    K_cv[ind] <- length(S)
    
    S=as.vector(S)
    S=append(S,0)
    S=append(S,Tb/2)
    S=sort(S)
    CP_cv[[ind]] <- S*2
    error=0
    for(j in 1:(length(S)-1)){
      for(t in (S[j]+1):S[j+1]){
        for(i in 1:m){
          indices=which(dat_test[,1]==(2*t-1))
          X_t=dat_test[indices,2:(d+1)]
          Y_t=dat_test[indices,(d+2):(d+2)]
          phat=p_hat(X_t[i],dat_train,h_bar,d)
          error=error+(g_hat_i(S[j],S[j+1],h_int,X_t[i],dat_train,h_bar,phat)-Y_t[i])^2
        }
      }
    }
    errors[1,ind]=error
    # }  
  }
  min_error=min(errors)
  tau_min=min(which(errors==min_error))
  # tau_min=tau_int[tau_min]
  K_matrix[1,l] <- K_cv[tau_min]
  print(CP_cv[tau_min])
  Hausdroff_distance[1,l] <-  Hausdorff.dist(CP_cv[[tau_min]], true_change)
}

print(K_matrix)
print(Hausdroff_distance)
result <- list(K_diff = K_matrix, dist = Hausdroff_distance)
# write.csv(result, "CV_setting1112_100_m40_FSBS_61to100.csv")
result

K_proportion_large <- sum(K_matrix > true_K_mat)/iteration
K_proportion_small <- sum(K_matrix < true_K_mat)/iteration
K_proportion_equal <- sum(K_matrix == true_K_mat)/iteration
mean.diff_K <- mean(abs(K_matrix - true_K_mat))
var.diff_K <- sqrt(var(abs(K_matrix - true_K_mat)[1,]))

mean.Hausdroff_distance <- mean(Hausdroff_distance/150)
var.Hausdroff.distance <- sqrt(var(Hausdroff_distance[1,]/150))

result = list(cpt_dist_mean = mean.Hausdroff_distance,
              cpt_dist_var = var.Hausdroff.distance,
              K_diff_mean = mean.diff_K,
              K_diff_var = var.diff_K,
              K_proportion_large = K_proportion_large,
              K_proportion_small = K_proportion_small,
              K_proportion_equal = K_proportion_equal)
print(result)