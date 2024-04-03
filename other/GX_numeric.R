lambda = 0.001
ext = 0.05
maxIt = 1
r1 = 2
r2 = 2
sigma1 = 1
sigma2 = 4
Delta = 50

# simulate data with one change point
data1 = temp_fragment_data3(mu = 0, r = r1, sigma1, n = 100, m = 10, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
data2 = temp_fragment_data4(mu = 0, r = r2, sigma2, n = 100, m = 10, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
data = list("t"= rbind(data1$t, data2$t), "y" = rbind(data1$y, data2$y), "r" = cbind(data1$r, data2$r))


regular_grid = seq(0.05, 0.95, length.out = 50)
# obtain the estimated covariance function
cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt) # estimate covariance function
C_est1 <- cov.obj1$C
C_est1
cov_est1 = predict_cov(regular_grid, C_est1)
# obtain the true covariance function
dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
C_1 <-  matrix(0, nrow = r1, ncol = r1)
diag(C_1) <- (1:r1)^(-1)
#C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
C_1
cov1 = predict_cov(regular_grid, C_1)
# compare
mean((cov_est1 - cov1)^2)

cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
C_est2 <- cov.obj2$C
C_est2
cov_est2 = predict_cov(regular_grid, C_est2)
dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
C_2 <-  matrix(0, nrow = r2, ncol = r2)
diag(C_2) <- (r2:1)^(-1)
C_2
cov2 = predict_cov(regular_grid, C_2)
mean((cov_est2 - cov2)^2)
#sum((cov_est2 - cov2)^2)/sum(cov2^2)

mean((cov1 - cov2)^2)
mean((cov_est1 - cov_est2)^2)


cov.obj = cov_basis(data$t, data$y, data$r, r1, lambda, ext, maxIt)
C_est <- cov.obj$C
cov_est = predict_cov(regular_grid, C_est, ext)
mean((cov_est - cov1)^2)
mean((cov_est - cov2)^2)
mean((cov_est - (cov1+cov2)/2)^2)

cpt_result = DP_fragment(data$t, data$y, data$r, r = 2, lambda, xi = 0.85, ext, maxIt, Delta=55)
part2local(cpt_result$partition)
