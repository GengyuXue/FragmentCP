library(FragmentCP)
library(changepoints)
library(data.table)
library(MASS)
library(fda)
library(tictoc)
library(rainbow)
library(sde)
library(xtable)
library(mvtnorm)
library(tseries)
library(expm)
library(tensorA)


### Theorem 2.1/2.2
# xdm - finite realization of DEMEANed functional time series data, where curves are stored in columns.
# u - a fraction index over the interval [0, 1]
ZNstat <- function(xdm, u){
  grid_point = nrow(xdm)
  N = ncol(xdm)
  k = floor(N*u)
  prek = matrix(rowSums(apply(as.matrix(xdm[,1:k]),2,function(x){x%o%x})),grid_point,grid_point)
  fullk = matrix(rowSums(apply(as.matrix(xdm),2,function(x){x%o%x})),grid_point,grid_point)
  ZNu = N^(-1/2) * (prek - (k/N)*fullk)
  return(ZNu)
}

ZNstat_cp <- function(xdm, u){
  grid_point = nrow(xdm)
  N = ncol(xdm)
  k = floor(N*u)
  prek = matrix(rowSums(apply(as.matrix(xdm[,1:k]),2,function(x){x%o%x})),grid_point,grid_point)
  fullk = matrix(rowSums(apply(as.matrix(xdm),2,function(x){x%o%x})),grid_point,grid_point)
  ZNu = (prek - (k/N)*fullk)
  return(ZNu)
}


# T_N statistic introduced after Theorem 2.1
# xf - finite realization of functional time series data, where curves are stored in columns.
TNstat <- function(xf){
  int_approx=function(x){
    temp_n=NROW(x)
    return((1/temp_n)*sum(x))}	
  grid_point = nrow(xf)
  N = ncol(xf)
  xdm = apply(xf,2,function(x,xmean){x-xmean}, xmean = rowMeans(xf))
  uind = seq(0,1,length = N+1)[2:(N+1)]; zn2=list()
  for (i in 1:N){
    zn2[[i]] = (ZNstat(xdm, uind[i]))^2
  }
  inm = Reduce(`+`, zn2)/N
  return((1/grid_point)^2*sum(inm))
}


# T_N(\kappa) statistic introduced after Theorem 2.3
weight_TNstat <- function(xf,kappa){
  int_approx_tensor<-function(x){# x is a 4-dimensional tensor
    dt=length(dim(x))
    temp_n=nrow(x)
    return(sum(x)/(temp_n^dt))}
  
  grid_point = nrow(xf)
  N = ncol(xf)
  xdm = apply(xf,2,function(x,xmean){x-xmean}, xmean = rowMeans(xf))
  uind = seq(0,1,length = N+1)[2:(N+1)]; 
  zn2 = list(); zn_cp = c(rep(0,N))
  for (i in 1:(N-1)){
    zn2[[i]] = (ZNstat(xdm, uind[i]))^2 / ((uind[i]*(1-uind[i]))^(2*kappa))### kappa = 1/4
    
    zn_cp[i] = (N/(i*(N-i)))^(kappa) *  int_approx_tensor( (ZNstat_cp(xdm, uind[i]))^2 )
  }
  inm = Reduce(`+`, zn2)/N
  stat = (1/grid_point)^2*sum(inm)
  
  mcp = max(zn_cp[ (0.1*N):(0.9*N)])
  changepoint = which(zn_cp == mcp)
  
  return(list (stat, changepoint) )
}


## Critical values

## useful functions for computing critical values
## long-run covariance operator
# dat - an array with dimension (grid_point,grid_point,N)
long_run_covariance_4tensor <- function (dat)
{
  grid_point = dim(dat)[1]
  T = dim(dat)[3]
  datmean = apply(dat, c(1,2), mean)
  center_dat = sweep(dat, 1:2, datmean)
  
  cov_l <- function(band, nval) {
    cov_sum = gamma_l(0, nval)
    
    for (ik in 1:(nval - 1)) {
      cov_sum = cov_sum + kweights(ik/band, kernel = "Bartlett") * (2*gamma_l(ik, nval))# + gamma_l(ik,nval))    ##aperm(gamma_l(ik,nval),c(2,1,3,4)))
    }
    return(cov_sum)
  }
  
  gamma_l <- function(lag, T) {
    gamma_lag_sum = 0
    if (lag >= 0) {
      for (ij in 1:(T - lag)) {
        gamma_lag_sum = gamma_lag_sum + center_dat[,,ij] %o% center_dat[,,(ij + lag)]
      }
    }
    else {
      for (ij in 1:(T + lag)) {
        gamma_lag_sum = gamma_lag_sum + center_dat[,,(ij - lag)] %o% center_dat[, ij]
      }
    }
    return(gamma_lag_sum/(T-lag))
  }
  hat_h_opt = T^(1/4)
  lr_covop = cov_l(band = hat_h_opt, nval = T)
  
  return(lr_covop)
}

kweights <- function (x, kernel = c("Truncated", "Bartlett", "Parzen", "Tukey-Hanning", 
                                    "Quadratic Spectral"), normalize = FALSE) 
{
  kernel <- match.arg(kernel)
  if (normalize) {
    ca <- switch(kernel, Truncated = 2, Bartlett = 2/3, Parzen = 0.539285, 
                 `Tukey-Hanning` = 3/4, `Quadratic Spectral` = 1)
  }
  else ca <- 1
  switch(kernel, Truncated = {
    ifelse(ca * x > 1, 0, 1)
  }, Bartlett = {
    ifelse(ca * x > 1, 0, 1 - abs(ca * x))
  }, Parzen = {
    ifelse(ca * x > 1, 0, ifelse(ca * x < 0.5, 1 - 6 * (ca * 
                                                          x)^2 + 6 * abs(ca * x)^3, 2 * (1 - abs(ca * x))^3))
  }, `Tukey-Hanning` = {
    ifelse(ca * x > 1, 0, (1 + cos(pi * ca * x))/2)
  }, `Quadratic Spectral` = {
    y <- 6 * pi * x/5
    ifelse(x < 1e-04, 1, 3 * (1/y)^2 * (sin(y)/y - cos(y)))
  })
}  

# compute critical values for TNstat (T_N)
criticalvalueMC <-function(xf,len){
  grid_point = nrow(xf)
  N = ncol(xf)
  
  rref = runif(len, 0, 1)
  rref = c(sort(rref), 1)
  rrefind = round(rref * dim(xf)[1])
  rrefind[which(rrefind==0)] = 1
  xfMC = xf[rrefind,]
  
  xdm = apply(xfMC,2,function(x,xmean){x-xmean}, xmean = rowMeans(xfMC))
  zi = zm = array(0, c((len+1), (len+1), N))
  for (i in 1:N){
    zi[,,i] = xdm[,i]%o%xdm[,i]
  }
  zimean = apply(zi,c(1,2),mean)
  for (i in 1:N){
    zm[,,i] = zi[,,i]-zimean
  }
  lrcov = long_run_covariance_4tensor(zm) ##23.883 sec elapsed
  lrcov = as.tensor(lrcov/(len+1)^2)
  eigvals=svd.tensor(lrcov,c(3,4),by="e")
  eigmat=as.vector(eigvals$d)
  
  lim_sum = 0
  for (ell in 1:length(eigmat)){
    klim = 0
    for (k in 1:1000){
      Nm=rnorm(2000,mean=0,sd=1)
      klim = klim + eigmat[ell]/((pi*k)^2)*Nm^2
    }
    lim_sum=lim_sum+klim
  }
  
  #lim_sum= rowSums(apply(matrix(seq(1,length(eigmat),1),1),2, function(x){ frac = eigmat[x]/((pi*seq(1,k,1))^2);
  #  rowSums(apply(matrix(seq(1,k,1),1),2,function(xx){frac[xx]*rnorm(5000,mean=0,sd=1)^2}))} ) )
  #klim = rowSums(t(frac*t(munor)))
  cv=quantile(lim_sum,probs=c(0.90,0.95,0.99))  
  return(cv)
}

# compute critical values for weight_TNstat ( T_N(\kappa) )

weight_criticalvalueMC <-function(xf,len,kappa){
  grid_point = nrow(xf)
  N = ncol(xf)
  
  ## cov weight function
  times = 1:grid_point/grid_point
  wmat = matrix(NA,grid_point-2,grid_point-2)
  for (i in 2:(grid_point-1)){
    for (j in 2:(grid_point-1)){
      wmat[i-1,j-1]= (min(times[i],times[j])-times[i]*times[j])/( (times[i]*(1-times[i]))^kappa * (times[j]*(1-times[j]))^kappa )
    }
  }
  weig = as.vector(svd(wmat/grid_point)$d)
  
  ## cov operators
  rref = runif(len, 0, 1)
  rref = c(sort(rref), 1)
  rrefind = round(rref * dim(xf)[1])
  rrefind[which(rrefind==0)] = 1
  xfMC = xf[rrefind,]
  
  xdm = apply(xfMC,2,function(x,xmean){x-xmean}, xmean = rowMeans(xfMC))
  zi = zm = array(0, c((len+1), (len+1), N))
  for (i in 1:N){
    zi[,,i] = xdm[,i]%o%xdm[,i]
  }
  zimean = apply(zi,c(1,2),mean)
  for (i in 1:N){
    zm[,,i] = zi[,,i]-zimean
  }
  
  lrcov = long_run_covariance_4tensor(zm)
  lrcov = as.tensor(lrcov/(len+1)^2)
  eigvals=svd.tensor(lrcov,c(3,4),by="e")
  eigmat=as.vector(eigvals$d)
  
  lim_sum = 0
  for (ell in 1:length(eigmat)){
    klim = 0
    for (k in 1:length(weig)){
      Nm=rnorm(2000,mean=0,sd=1)
      klim = klim + eigmat[ell]*weig[k]*Nm^2
    }
    lim_sum=lim_sum+klim
  }
  
  cv=quantile(lim_sum,probs=c(0.90,0.95,0.99))  
  return(cv)
}

##### The results in Section 4. Theorem 4.1
################## some useful functions for the limit distribution of estimator
mkfun <- function(kappa,theta,tim){
  if (tim < 0){
    mk = (1-kappa)*(1-theta)+kappa*theta
  }
  else if (tim > 0){
    mk = (1-kappa)*theta+kappa*(1-theta)
  }
  else {
    mk = 0
  }
  return(mk)
}


twowiener <- function(lim_N){
  times1 = seq(-1,0,length=lim_N/2)
  times2 = seq(0,1,length=lim_N/2)
  times = c(times1[1:(lim_N/2-1)], times2)
  
  w1 = BM(x=0, t0=0, T=1, N=lim_N/2-1)
  w2 = BM(x=0, t0=0, T=1, N=lim_N/2-1)
  w = c(rev(w1[2:(lim_N/2)]), w2)
  return(list(times,w))
}


lim_cp <- function(reps,kappa,theta){
  N = 10000
  cpseq = c(rep(NA,reps))
  for (j in 1:reps){
    wtime = twowiener(N)
    times = wtime[[1]]
    w = wtime[[2]]
    wseq = c(rep(NA,N-1))
    for (i in 1:(N-1)){
      wseq[i] = w[i]-abs(times[i])*mkfun(kappa,theta,times[i])
    }
    mcp = max(wseq)
    limcp = which(wseq == mcp)
    cpseq[j] = limcp
  }
  return(cpseq)
}




###########  estimator


##### estimate size of change
sizechange <- function(xd, kstar){
  N=ncol(xd)
  
  sample_cov<-function(data){
    N=ncol(data)
    varmtx = 0
    for (i in 1:N){
      varmtx = varmtx + data[,i]%o%data[,i]
    }
    return(varmtx/N)
  }
  
  error = apply(xd,2,function(x,xmean){x-xmean}, xmean = rowMeans(xd))
  error_before = error[,1:kstar]
  error_after = error[,(kstar+1):N]
  
  var_before = sample_cov(error_before)
  var_after = sample_cov(error_after)
  var_change = var_before - var_after
  return(var_change)
}

########

l2norm <-function(vec){
  return(sqrt(sum(vec^2)))
}

tau_est <-function(xd, kstar, len){
  grid_point = nrow(xd)
  N = ncol(xd)
  
  rref = runif(len, 0, 1)
  rref = c(sort(rref), 1)
  rrefind = round(rref * grid_point)
  rrefind[which(rrefind==0)] = 1
  xdmc = xd[rrefind,]
  
  sample_cov<-function(data){
    N=ncol(data)
    varmtx = 0
    for (i in 1:N){
      varmtx = varmtx + data[,i]%o%data[,i]
    }
    return(varmtx/N)
  }
  
  error = apply(xdmc,2,function(x,xmean){x-xmean}, xmean = rowMeans(xdmc))
  error_before = error[,1:kstar]
  error_after = error[,(kstar+1):N]
  
  var_before = sample_cov(error_before)
  var_after = sample_cov(error_after)
  var_change = var_before - var_after
  
  ## change star
  var_1 = var_2 = 0
  
  for (i in 1:kstar){
    var_1 = var_1 + (xdmc[,i]-rowMeans(xdmc))%o%(xdmc[,i]-rowMeans(xdmc))
  }
  var_1 = 1/kstar * var_1
  
  for (i in (kstar+1):N){
    var_2 = var_2 + (xdmc[,i]-rowMeans(xdmc))%o%(xdmc[,i]-rowMeans(xdmc))
  }
  var_2 = 1/(N-kstar) * var_2
  
  var_star = (var_1-var_2)/l2norm(var_1-var_2)
  
  ## longrun cov
  
  zi = zm = array(0, c((len+1), (len+1), N))
  for (i in 1:N){
    zi[,,i] = error[,i]%o%error[,i]
  }
  
  v_dat=array(0, c(len+1, len+1, N))
  for (i in 1:N){
    if (i <= kstar){
      v_dat[,,i] = zi[,,i]-var_1
    }else{
      v_dat[,,i] = zi[,,i]-var_2
    }
  }
  
  int_approx_tensor<-function(x){# x is a 4-dimensional tensor
    dt=length(dim(x))
    temp_n=nrow(x)
    return((1/temp_n)^dt * sum(x))}
  
  longd = long_run_covariance_4tensor(v_dat)
  
  frontvs = rearvs = 0
  for (i in 1:21){
    for (j in 1:21){
      frontvs = frontvs + var_star%o%longd[i,,j,]
    }
  }
  for (i in 1:21){
    for (j in 1:21){
      rearvs = rearvs + frontvs[,i,,j]%o%var_star
    }
  }
  tau = int_approx_tensor(rearvs)
  
  return(list(var_change,tau))
}

##########################################################################
##########################################################################

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
##########################################################################
##########################################################################
lambda = 0.00001
ext = 0.1
maxIt = 1
r1 = 3
r2 = 3
sigma1 = 1
sigma2 = 1

true_K <- 1
true_change <- 100
true_K_mat <- matrix(true_K, nrow = 1, ncol = 100)
K_matrix <- matrix(0, nrow = 1, ncol = 100)
Hausdroff_distance <- matrix(0, nrow = 1, ncol = 100)

change_est <- matrix(0, nrow = 1, ncol = 100)

smooth_function <- function(data, num_basis = 5, rangeval = c(0,1), no_grids = 30){
  # change format of data
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
  
  # create the data matrix containing t, y and labels
  dat[, 1] <- labels_Ts
  dat[, 2] <- t_vector
  dat[, 3] <- y_vector
  
  # create data used for smoothing
  dat_y=matrix(0,m,Tb)
  aux_count=0
  for(t in 1:Tb){
    aux_count=m+aux_count
    for(i in 1:m){
      dat_y[i,t]=y_vector[aux_count-m+i]
    }
  }
  # dat_y
  
  dat_t=matrix(0,m,Tb)
  aux_count=0
  for(t in 1:Tb){
    aux_count=m+aux_count
    for(i in 1:m){
      dat_t[i,t]=t_vector[aux_count-m+i]
    }
  }
  # dat_t
  
  # smoothing using spline basis
  basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=num_basis, norder=4)
  smooth_func <- smooth.basis(dat_t, dat_y, basis)$fd
  eval_points <- seq(0, 1, length.out = no_grids)
  eval <- as.matrix(eval.fd(eval_points , smooth_func))
  return(eval)
}

kappa = 0.25
lens = 9
reps = 10
sigma1 = 1
sigma2 = 5

for (i in 1:100){
  print(c("iteration", i))
  set.seed(1000+10*i)
  
  data1 = temp_fragment_data11(mu = 0, sigma1, n = 100, m = 40, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.6)
  data2 = temp_fragment_data12(mu = 0, sigma2, n = 100, m = 40, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.6)
  data = list("t"= rbind(data1$t, data2$t), "y" = rbind(data1$y, data2$y), "r" = cbind(data1$r, data2$r))
  
  data_smooth <- smooth_function(data, num_basis = 4, rangeval = c(0,1), no_grids = 50)
  
  stat_d0 = weight_TNstat(data_smooth, kappa)
  cv_d0 = weight_criticalvalueMC(data_smooth,len = lens,kappa)
  # print(stat_d0[[2]])
  if (stat_d0[[1]]> cv_d0[2]){
    change_est[1,i] = stat_d0[[2]]
    K_matrix[1,i] = 1
    
  }else{
    change_est[1,i] = 0
  }
  print(change_est[1,i])
}

print(K_matrix)

for (i in 1:100) {
  Hausdroff_distance[1,i] <-  Hausdorff.dist(change_est[1,i], true_change)
}

K_proportion_large <- sum(K_matrix > true_K_mat)/100
K_proportion_small <- sum(K_matrix < true_K_mat)/100
K_proportion_equal <- sum(K_matrix == true_K_mat)/100
mean.diff_K <- mean(abs(K_matrix - true_K_mat))
var.diff_K <- sqrt(var(abs(K_matrix - true_K_mat)[1,]/300))

mean.Hausdroff_distance <- mean(Hausdroff_distance)/300
var.Hausdroff.distance <- sqrt(var(Hausdroff_distance[1,]/300))

result = list(cpt_dist_mean = mean.Hausdroff_distance,
              cpt_dist_var = var.Hausdroff.distance,
              K_diff_mean = mean.diff_K,
              K_diff_var = var.diff_K,
              K_proportion_large = K_proportion_large,
              K_proportion_small = K_proportion_small,
              K_proportion_equal = K_proportion_equal)
print(result)

write.csv(result, "setting1112_WCUSUM_100_m40.csv")


