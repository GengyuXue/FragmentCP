#' @export
CUSUM_fragment = function(Lt, Ly, Lr, r, s, e, t, lambda, ext, maxIt){
  m = ncol(Lt)
  Lt_se = Lt[(s+1):e,]
  Ly_se = Ly[(s+1):e,]
  Lr_se = Lr[,as.vector(sapply((s+1):e, function(i){(m*(i-1)+1):(m*i)}))]
  Lt_st = Lt[(s+1):t,]
  Ly_st = Ly[(s+1):t,]
  Lr_st = Lr[,as.vector(sapply((s+1):t, function(i){(m*(i-1)+1):(m*i)}))]
  Lt_te = Lt[(t+1):e,]
  Ly_te = Ly[(t+1):e,]
  Lr_te = Lr[,as.vector(sapply((t+1):e, function(i){(m*(i-1)+1):(m*i)}))]
  error_se = cov_basis(Lt_se, Ly_se, Lr_se, r, lambda, ext, maxIt)$error
  error_st = cov_basis(Lt_st, Ly_st, Lr_st, r, lambda, ext, maxIt)$error
  error_te = cov_basis(Lt_te, Ly_te, Lr_te, r, lambda, ext, maxIt)$error
  result = error_se - error_st - error_te
  return(result)
}

## seeded intervals
#' @export
seeded_intervals = function(n, Delta){
  M = ceiling(log2(n/Delta))+1
  n_vec = c(2^c(1:M) - 1)
  s_vec = n/2^c(1:M)
  l_vec = 2*s_vec
  output = NULL
  for(k in 1:M){
    for(i in 1:n_vec[k]){
      output = rbind(output, c(ceiling((i-1)*s_vec[k]), floor((i-1)*s_vec[k]+l_vec[k])))
    }
  }
  return(output)
}

## SBS
#' @export
SBS_fragment = function(Lt, Ly, Lr, r, lambda, ext, maxIt, s, e, Alpha, Beta, Delta, level = 0){
  # print(paste0("SBS at level: ", level))
  Alpha_new = pmax(Alpha, s)
  Beta_new = pmin(Beta, e)
  idx = which(Beta_new - Alpha_new > 2*Delta)
  Alpha_new = Alpha_new[idx]
  Beta_new = Beta_new[idx]
  M = length(Alpha_new)
  S = NULL
  Dval = NULL
  Level = NULL
  Parent = NULL
  if(M == 0){
    return(list(S = S, Dval = Dval, Level = Level, Parent = Parent))
  }else{
    level = level + 1
    parent = matrix(c(s, e), nrow = 2)
    a = rep(0, M)
    b = rep(0, M)
    for(m in 1:M){
      s_star = Alpha_new[m] + Delta
      e_star = Beta_new[m] - Delta
      temp = rep(0, e_star - s_star)
      for(t in s_star:e_star){
        temp[t-s_star] = CUSUM_fragment(Lt, Ly, Lr, r, Alpha_new[m], Beta_new[m], t, lambda, ext, maxIt)
      }
      best_value = max(temp)
      best_t = which.max(temp) + s_star
      a[m] = best_value
      b[m] = best_t
    }
    m_star = which.max(a)
  }
  temp1 = SBS_fragment(Lt, Ly, Lr, r, lambda, ext, maxIt, s, b[m_star]-1, Alpha, Beta, Delta, level)
  temp2 = SBS_fragment(Lt, Ly, Lr, r, lambda, ext, maxIt, b[m_star], e, Alpha, Beta, Delta, level)
  S = c(temp1$S, b[m_star], temp2$S)
  Dval = c(temp1$Dval, a[m_star], temp2$Dval)
  Level = c(temp1$Level, level, temp2$Level)
  Parent = cbind(temp1$Parent, parent, temp2$Parent)
  result = list(S = S, Dval = Dval, Level = Level, Parent = Parent)
  class(result) = "BS"
  return(result)
}


# cross validation
#' @export
CV_SBS_fragment = function(Lt, Ly, Lr, r, lambda, ext, maxIt, zeta, Delta = 20){
  n = nrow(Lt)
  m = ncol(Lt)
  even_indexes = seq(2, n, 2)
  odd_indexes = seq(1, n, 2)
  train_Lt = Lt[odd_indexes,]
  train_Ly = Ly[odd_indexes,]
  train_Lr = Lr[, as.vector(sapply(odd_indexes, function(i){(m*(i-1)+1):(m*i)}))]
  validation_Lt = Lt[even_indexes,]
  validation_Ly = Ly[even_indexes,]
  validation_Lr = Lr[, as.vector(sapply(even_indexes, function(i){(m*(i-1)+1):(m*i)}))]
  s_intervals = seeded_intervals(length(odd_indexes), 30)
  init_train = SBS_fragment(train_Lt, train_Ly, train_Lr, r, lambda, ext, maxIt, 0, length(odd_indexes), s_intervals[,1], s_intervals[,2], Delta)
  init_train_cpt = thresholdBS(init_train, zeta)$cpt_hat[,1]
  if(length(init_train_cpt) >= 1){
    init_train_cpt_long = c(0, init_train_cpt, length(odd_indexes))
    train_error = 0
    test_error = 0
    init_train_beta = NULL
    for(k in 1:(length(init_train_cpt)+1)){
      train_Lt_temp = train_Lt[(init_train_cpt_long[k]+1):(init_train_cpt_long[k+1]),]
      train_Ly_temp = train_Ly[(init_train_cpt_long[k]+1):(init_train_cpt_long[k+1]),]
      train_Lr_temp = train_Lr[,as.vector(sapply((init_train_cpt_long[k]+1):(init_train_cpt_long[k+1]), function(i){(m*(i-1)+1):(m*i)}))]
      test_Lt_temp = validation_Lt[(init_train_cpt_long[k]+1):(init_train_cpt_long[k+1]),]
      test_Ly_temp = validation_Ly[(init_train_cpt_long[k]+1):(init_train_cpt_long[k+1]),]
      test_Lr_temp = validation_Lr[,as.vector(sapply((init_train_cpt_long[k]+1):(init_train_cpt_long[k+1]), function(i){(m*(i-1)+1):(m*i)}))]
      train_temp = cov_basis(train_Lt_temp, train_Ly_temp, train_Lr_temp, r, lambda, ext, maxIt)
      train_error = train_error + train_temp$error
      test_error = test_error + error_test_fragment(test_Lt_temp, test_Lr_temp, r, 1, nrow(test_Lt_temp), train_temp$C)
    }
    init_cpt = odd_indexes[init_train_cpt]
    K_hat = length(init_train_cpt)
  }else{
    init_cpt = init_train_cpt
    K_hat = 0
    train_temp = cov_basis(train_Lt, train_Ly, train_Lr, r, lambda, ext, maxIt)
    train_error = train_temp$error
    test_error = error_test_fragment(validation_Lt, validation_Lr, r, 1, ncol(validation_Lt), train_temp$C)
  }
  result = list(cpt_hat = init_cpt, K_hat, test_error = test_error, train_error = train_error)
  return(result)
}

#' @export
CV_search_SBS_fragment = function(Lt, Ly, Lr, r, lambda, ext, maxIt, zeta_set, Delta = 20){
  output = sapply(1:length(zeta_set), function(j) CV_SBS_fragment(Lt, Ly, Lr, r, lambda, ext, maxIt, zeta_set[j], Delta))
  cpt_hat = output[seq(1,4,4),]## estimated change points
  K_hat = output[seq(2,4,4),]## number of estimated change points
  test_error = output[seq(3,4,4),]## validation loss
  train_error = output[seq(4,4,4),]## training loss                                                      
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
  return(result)
}
