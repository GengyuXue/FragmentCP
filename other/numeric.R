lambda = 0.00001
ext = 0.05
maxIt = 1
r1 = 5
r2 = 7
sigma1 = 1
sigma2 = 1
Delta = 20

data1 = fragment_data(mu = 0, cov = periodic_cov_fct, r = r1, sigma1, n = 75, m = 5, sigma_epsilon = 0.005, domain = c(0, 1), delta = 0.3)
data2 = fragment_data(mu = 0, cov = periodic_cov_fct, r = r2, sigma2, n = 75, m = 5, sigma_epsilon = 0.005, domain = c(0, 1), delta = 0.3)

data$t = rbind(data1$t, data2$t)
data$y = rbind(data1$y, data2$y)
data$r = cbind(data1$r, data2$r)

cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
C_est1 <- cov.obj1$C


dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))

cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
C_est2 <- cov.obj2$C
dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
 

cpt_result = DP_fragment(data$t, data$y, data$r, r = 5, lambda, xi = 20, ext, maxIt, Delta)
part2local(cpt_result$partition)
