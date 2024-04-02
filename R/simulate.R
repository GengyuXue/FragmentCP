# function to generate functional data observed on fragments
#' @export
fragment_data <- function(mu = 0, cov = periodic_cov_fct, r = 5, sigma, n = 100, m = 5, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3){
  Ly <- matrix(0, n, m)
  L <- domain[2] - domain[1]
  a_vec <- runif(n, min = domain[1], max = domain[2] - delta * L)
  Lt <- t(sapply(1:n, function(i) {
    a <- a_vec[i]
    a + sort(runif(m, min = 0, max = delta * L))
  }))
  Ly0 <- t(apply(Lt, 1, function(tobs) {
    if(is.function(mu)) 
      mui <- mu(tobs)
    else if(length(mu) == 1) 
      mui <- rep(mu, length(tobs))
    else stop("mu must be a scalar or a function.")
    y0 <- mui + gaussian_process(cov, 1, tobs, r, sigma)
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


#dist_mat = dist(1:r, method = "manhattan", diag = T, upper = T)
#C_mat = as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r)))

#' @export
periodic_cov_fct1 <- function(tObs, r = 5){
  dist_mat = dist(1:r, method = "manhattan", diag = T, upper = T)
  C_mat = as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r)))
  temp = evaluate_basis(r, c(0,1), tObs)
  C = temp %*% C_mat %*% t(temp)
  return(C)
}

#' @export
periodic_cov_fct2 <- function(tObs, r = 5){
  dist_mat = dist(1:r, method = "manhattan", diag = T, upper = T)
  C_mat = as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r)))
  temp = evaluate_basis(r, c(0,1), tObs)
  C = temp %*% C_mat %*% t(temp)
  return(C)
}


shrink <- function(Lt,ext){
  ext + (1-ext*2)*Lt
}

#' @export
predict_cov <- function(grids, C_mat, ext){
  r <- dim(C_mat)[1]
  if(ext > 0) grids <- shrink(grids, ext)
  B <- evaluate_basis(r, c(0,1), grid = grids)
  return(B %*% C_mat %*% t(B))
}

# white_fct <- function(tObs){
#   C = diag(rep(1, length(tObs)))
#   return(C)
# }


#' @export
gaussian_process <- function(cov_fct, n, tObs, r, sigma){
  m = length(tObs)
  Sigma <- cov_fct(tObs, r)
  result = MASS::mvrnorm(n = n, mu = rep(0, m), Sigma = Sigma*sigma^2)
  return(result)
}

#' @export
part2local = function(parti_vec){
  N = length(parti_vec)
  localization = c()
  r = N
  l = parti_vec[r]
  localization = c(l, localization)
  while(r > 0){
    r = l
    l = parti_vec[r]
    localization = c(l, localization)
  }
  return(localization[-1])
}
