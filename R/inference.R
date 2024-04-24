# simulate two-sided Brownian motion with drift
#' @export
simu.2BM_Drift = function(n, drift, LRV1, LRV2){
  # 1: before 2:after
  z_vec = rnorm(2*n)
  w_vec1 = rev(cumsum(z_vec[n:1])/sqrt(1:n))
  w_vec2 = cumsum(z_vec[(n+1):(2*n)])/sqrt(1:n)
  
  v_vec1 = drift*abs(seq(-n, -1)) + sqrt(LRV1)*w_vec1
  v_vec2 = drift*abs(seq(1, n)) + sqrt(LRV2)*w_vec2
  
  v_vec = c(v_vec1, 0, v_vec2)
  # w_vec = c(rev(cumsum(z_vec[n:1])/sqrt(1:n)), 0, cumsum(z_vec[(n+1):(2*n)])/sqrt(1:n))
  # v_vec = drift*abs(seq(-n, n)) + sqrt(LRV)*w_vec
  return(v_vec)
}

# refined change point estimator
#' @export
local_refine_fragment = function(Lt, Ly, Lr, cpt_init_hat, r, lambda, ext, maxIt, w = 0.9){
  n = nrow(Lt)
  m = ncol(Lt)
  Khat = length(cpt_init_hat)
  cpt_hat_long = c(0, cpt_init_hat, n)
  cpt_refined = rep(0, Khat+1)
  C_hat_list <- vector("list", Khat+1)
  for(i in 1:(Khat+1)){
    Lt_se = Lt[(cpt_hat_long[i]+1):cpt_hat_long[i+1],]
    Ly_se = Ly[(cpt_hat_long[i]+1):cpt_hat_long[i+1],]
    Lr_se = Lr[,as.vector(sapply((cpt_hat_long[i]+1):cpt_hat_long[i+1], function(i){(m*(i-1)+1):(m*i)}))]
    C_hat_list[[i]] = cov_basis(Lt_se, Ly_se, Lr_se, r, lambda, ext, maxIt)$C
  }
  for (k in 1:Khat){
    s = ceiling(w*cpt_hat_long[k] + (1-w)*cpt_hat_long[k+1])
    e = floor((1-w)*cpt_hat_long[k+1] + w*cpt_hat_long[k+2])
    lower = s + 2
    upper = e - 2
    # Lt_se = Lt[s:eta,]
    # Ly_se = Ly[s:eta,]
    # Lr_se = Lr[,as.vector(sapply(s:eta, function(i){(m*(i-1)+1):(m*i)}))]
    # b = sapply(lower:upper, function(eta)(cov_basis(Lt[s:eta,], Ly[s:eta,], Lr[,as.vector(sapply(s:eta, function(i){(m*(i-1)+1):(m*i)}))], r, lambda, ext, maxIt)$error 
    # + cov_basis(Lt[(eta+1):e,], Ly[(eta+1):e,], Lr[,as.vector(sapply((eta+1):e, function(i){(m*(i-1)+1):(m*i)}))], r, lambda, ext, maxIt)$error))
    
    b = sapply(lower:upper, function(eta) (error_test_fragment(Lt, Lr, r, s, eta,C_hat_list[[k]]) 
                                           + error_test_fragment(Lt, Lr, r, eta+1, e, C_hat_list[[k+1]])))
    cpt_refined[k+1] = s + which.min(b)
  }
  return(cpt_refined[-1])
}

# estimate the jump size
#' @export
kappa2_fragment = function(Lt, Ly, Lr, cpt_hat, r, lambda, ext, maxIt){
  n = nrow(Lt)
  m = ncol(Lt)
  Khat = length(cpt_hat)
  cpt_hat_long = c(0, cpt_hat, n)
  C_hat_list = vector("list", Khat+1)
  kappa2_hat_vec = rep(NA, Khat)
  for(i in 1:(Khat+1)){
    Lt_se = Lt[(cpt_hat_long[i]+1):cpt_hat_long[i+1],]
    Ly_se = Ly[(cpt_hat_long[i]+1):cpt_hat_long[i+1],]
    Lr_se = Lr[,as.vector(sapply((cpt_hat_long[i]+1):cpt_hat_long[i+1], function(i){(m*(i-1)+1):(m*i)}))]
    C_hat_list[[i]] = cov_basis(Lt_se, Ly_se, Lr_se, r, lambda, ext, maxIt)$C
  }
  #regular_grid = seq(0.05, 0.95, length.out = 400)
  #basis_mat = evaluate_basis(r, c(0,1), Lt[i,])
  for(i in 1:Khat){
    kappa2_hat_vec[i] = sum((C_hat_list[[i+1]] - C_hat_list[[i]])^2) #(1/400^2)*predict_cov(regular_grid, C_hat_list[[i+1]] - C_hat_list[[i]], ext = 0)
  }
return(list(kappa2=kappa2_hat_vec, C=C_hat_list))
}

# estimate sigma2
#' @export
sigma2_fragment = function(cpt_hat, kappa2_hat, C_list, Lt, Ly, Lr){
  n = nrow(Lt)
  m = ncol(Lt)
  r = dim(C_list[[1]])[1]
  Khat = length(cpt_hat)
  cpt_long = c(0, cpt_hat, n)
  sig2_before_hat = rep(NA, Khat)
  sig2_after_hat = rep(NA, Khat)
  for (k in 1:Khat){
    s = cpt_long[k]
    t = cpt_long[k+1]
    e = cpt_long[k+2]
    kappa2 = kappa2_hat[k]
    C1_hat = C_list[[k]]
    C2_hat = C_list[[k+1]]
    Lp1 = matrix(NA, m, m)
    Lp2 = matrix(NA, m, m)
    Z1 = matrix(NA, m, m)
    Z2 = matrix(NA, m, m)
    z_before_vec = rep(NA, t-s)
    z_after_vec = rep(NA, e-t)
    for(j in (s+1):t){
      basis_mat = evaluate_basis(r, c(0,1), Lt[j,])
      Lp1 = basis_mat %*% C1_hat %*% t(basis_mat)
      Lp2 = basis_mat %*% C2_hat %*% t(basis_mat)
      Z1 = Lr[,(m*(j-1)+1):(m*j)] - Lp1
      z_before_vec[j-s] = sum(Z1*(Lp2-Lp1)/sqrt(kappa2))
    }
    sig2_before_hat[k] = 4*var(z_before_vec)
    for(j in (t+1):e){
      basis_mat = evaluate_basis(r, c(0,1), Lt[j,])
      Lp1 = basis_mat %*% C1_hat %*% t(basis_mat)
      Lp2 = basis_mat %*% C2_hat %*% t(basis_mat)
      Z2 = Lr[,(m*(j-1)+1):(m*j)] - Lp2
      z_after_vec[j-t] = sum(Z2*(Lp2-Lp1)/sqrt(kappa2))
    }
    sig2_after_hat[k] = 4*var(z_after_vec)
  }
  return(list(sig2_before = sig2_before_hat, sig2_after = sig2_after_hat))
}

# estimate omega
#' @export
omega_fragment = function(cpt_hat, kappa2_hat, C_list, Lt, Ly, Lr){
  n = nrow(Lt)
  m = ncol(Lt)
  r = dim(C_list[[1]])[1]
  Khat = length(cpt_hat)
  cpt_long = c(0, cpt_hat, n)
  omega_hat = rep(NA, Khat)
  for (k in 1:Khat){
    s = cpt_long[k]
    t = cpt_long[k+1]
    e = cpt_long[k+2]
    kappa2 = kappa2_hat[k]
    C1_hat = C_list[[k]]
    C2_hat = C_list[[k+1]]
    Lp1 = matrix(NA, m, m)
    Lp2 = matrix(NA, m, m)
    z_vec = rep(NA, e-s)
    for(j in (s+1):e){
      basis_mat = evaluate_basis(r, c(0,1), Lt[j,])
      Lp1 = basis_mat %*% C1_hat %*% t(basis_mat)
      Lp2 = basis_mat %*% C2_hat %*% t(basis_mat)
      z_vec[j-s] = sum((Lp2-Lp1)^2/kappa2)
    }
    omega_hat[k] = mean(z_vec)
  }
  return(omega_hat)
}




