lambda = 0.0001
ext = 0
maxIt = 1
r1 = 5
r2 = 3
Delta = 10

data1 = fragment_data(mu = 0, cov = periodic_cov_fct, r = r1, n = 100, m = 3, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3)
data2 = fragment_data(mu = 0, cov = periodic_cov_fct, r = r2, n = 100, m = 3, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3)

data$t = rbind(data1$t, data2$t)
data$y = rbind(data1$y, data2$y)
data$r = cbind(data1$r, data2$r)

#cov.obj = cov_basis(data$t, data$y, data$r, r, lambda, ext, maxIt)
#C_est <- cov.obj$C

DP_fragment(data$t, data$y, data$r, r = 5, lambda, xi = 1, ext, maxIt, Delta)
