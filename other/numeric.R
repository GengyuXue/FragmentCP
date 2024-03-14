library(synfd)
X_fct = centered.process(name='wiener',dispersion=0.001)
data = fragment_data(mu = 0, X_fct = X_fct, n = 200, m = 3, sigma_epsilon = 1, domain = c(0, 1), delta = 0.3)
r = 5
lambda = 0.0001
ext = 0.05
cov.obj = cov_basis(data$t, data$y, data$r, r, lambda, ext, 2)
C_hat <- cov.obj$C
auxmat = auxiliary_mat(r, data$t)
DP_fragment(data$t, data$y, data$r, 5, 0.0001, 1, 0.05, 2, 10)
