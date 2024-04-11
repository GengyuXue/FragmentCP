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
  print(paste0("SBS at level: ", level))
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