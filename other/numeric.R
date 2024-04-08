lambda = 0.000001
ext = 0.05
maxIt = 1
r1 = 3
r2 = 3
sigma1 = 1
sigma2 = 1
Delta = 50

# simulate data with one change point
data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 200, m = 200, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.6)
data2 = fragment_data2(mu = 0, r = r2, sigma2, n = 200, m = 200, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.6)
data = list("t"= rbind(data1$t, data2$t), "y" = rbind(data1$y, data2$y), "r" = cbind(data1$r, data2$r))


regular_grid = seq(0.05, 0.95, length.out = 100)
# obtain the estimated covariance function
cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt) # estimate covariance function
C_est1 <- cov.obj1$C
cov_est1 = predict_cov(regular_grid, C_est1)
# obtain the true covariance function
dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
cov1 = predict_cov(regular_grid, C_1)
# compare
mean((C_est1 - C_1)^2)
mean((cov_est1 - cov1)^2)
#cov.obj <- covfunc(as.list(data.frame(t(data1$t))),as.list(data.frame(t(data1$y))),newt=NULL,mu=NULL,method='FOURIER', tuning='lle',weig=NULL)
#cov1.hat <- predict(cov.obj,regular_grid)
#mean((cov.obj$C[1:3, 1:3] - C_1)^2)
#mean((cov1.hat - cov1)^2)

cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
C_est2 <- cov.obj2$C
cov_est2 = predict_cov(regular_grid, C_est2)
dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(r2:1)))
cov2 = predict_cov(regular_grid, C_2)
mean((C_est2 - C_2)^2)
mean((cov_est2 - cov2)^2)
#sum((cov_est2 - cov2)^2)/sum(cov2^2)
#cov.obj <- covfunc(as.list(data.frame(t(data2$t))),as.list(data.frame(t(data2$y))),newt=NULL,mu=NULL,method='FOURIER', tuning='lle',weig=NULL)
#cov2.hat <- predict(cov.obj,regular_grid)
#mean((cov.obj$C[1:3, 1:3] - C_2)^2)
#mean((cov2.hat - cov2)^2)

mean((C_1 - C_2)^2)
mean((cov_est1 - cov_est2)^2)
#mean((cov1.hat - cov2.hat)^2)
mean((cov1 - cov2)^2)
#mean((cov_est1 - cov_est2)^2)
#mean((cov1.hat - cov2.hat)^2)


cov.obj = cov_basis(data$t, data$y, data$r, r1, lambda, ext, maxIt)
C_est <- cov.obj$C
cov_est = predict_cov(regular_grid, C_est, ext)
mean((cov_est - cov1)^2)
mean((cov_est - cov2)^2)
mean((cov_est - (cov1+cov2)/2)^2)

cpt_result = DP_fragment(data$t, data$y, data$r, r = 3, lambda, xi = 0.5, ext, maxIt, Delta)
part2local(cpt_result$partition)


B = 50
lambda = 0.00001
regular_grid = seq(0.05, 0.95, length.out = 50)
results_n100_m5 = matrix(NA, B, 2)
results2_n100_m5 = matrix(NA, B, 2)
for(b in 1:B){
  print(b)
  data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 100, m = 5, sigma_epsilon = 0.0001, domain = c(0, 1), delta = 0.3)
  cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
  C_est1 <- cov.obj1$C
  cov_est1 = predict_cov(regular_grid, C_est1)
  dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
  C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
  cov1 = predict_cov(regular_grid, C_1)
  results_n100_m5[b,1] = cov.obj1$error
  results_n100_m5[b,2] = mean((cov_est1 - cov1)^2)
  cov.obj <- covfunc(as.list(data.frame(t(data1$t))),as.list(data.frame(t(data1$y))),newt=NULL,mu=NULL,method='FOURIER', tuning='lle',weig=NULL)
  cov.hat <- predict(cov.obj,regular_grid)
  results2_n100_m5[b,2] = mean((cov.hat - cov1)^2)
}


results_n100_m10 = matrix(NA, B, 2)
for(b in 1:B){
  data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 100, m = 10, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
  cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
  C_est1 <- cov.obj1$C
  cov_est1 = predict_cov(regular_grid, C_est1)
  dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
  C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
  cov1 = predict_cov(regular_grid, C_1)
  results_n100_m10[b,1] = cov.obj1$error
  results_n100_m10[b,2] = mean((cov_est1 - cov1)^2)
}
results_n100_m15 = matrix(NA, B, 2)
for(b in 1:B){
  data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 100, m = 15, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
  cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
  C_est1 <- cov.obj1$C
  cov_est1 = predict_cov(regular_grid, C_est1)
  dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
  C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
  cov1 = predict_cov(regular_grid, C_1)
  results_n100_m15[b,1] = cov.obj1$error
  results_n100_m15[b,2] = mean((cov_est1 - cov1)^2)
}
results_n100_m20 = matrix(NA, B, 2)
for(b in 1:B){
  data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 100, m = 20, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
  cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
  C_est1 <- cov.obj1$C
  cov_est1 = predict_cov(regular_grid, C_est1)
  dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
  C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
  cov1 = predict_cov(regular_grid, C_1)
  results_n100_m20[b,1] = cov.obj1$error
  results_n100_m20[b,2] = mean((cov_est1 - cov1)^2)
}
results_n150_m5 = matrix(NA, B, 2)
for(b in 1:B){
  data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 150, m = 5, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
  cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
  C_est1 <- cov.obj1$C
  cov_est1 = predict_cov(regular_grid, C_est1)
  dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
  C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
  cov1 = predict_cov(regular_grid, C_1)
  results_n150_m5[b,1] = cov.obj1$error
  results_n150_m5[b,2] = mean((cov_est1 - cov1)^2)
}
results_n150_m10 = matrix(NA, B, 2)
for(b in 1:B){
  data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 150, m = 10, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
  cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
  C_est1 <- cov.obj1$C
  cov_est1 = predict_cov(regular_grid, C_est1)
  dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
  C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
  cov1 = predict_cov(regular_grid, C_1)
  results_n150_m10[b,1] = cov.obj1$error
  results_n150_m10[b,2] = mean((cov_est1 - cov1)^2)
}
results_n150_m15 = matrix(NA, B, 2)
for(b in 1:B){
  data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 150, m = 15, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
  cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
  C_est1 <- cov.obj1$C
  cov_est1 = predict_cov(regular_grid, C_est1)
  dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
  C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
  cov1 = predict_cov(regular_grid, C_1)
  results_n150_m15[b,1] = cov.obj1$error
  results_n150_m15[b,2] = mean((cov_est1 - cov1)^2)
}
results_n150_m20 = matrix(NA, B, 2)
for(b in 1:B){
  data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 150, m = 20, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
  cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
  C_est1 <- cov.obj1$C
  cov_est1 = predict_cov(regular_grid, C_est1)
  dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
  C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
  cov1 = predict_cov(regular_grid, C_1)
  results_n150_m20[b,1] = cov.obj1$error
  results_n150_m20[b,2] = mean((cov_est1 - cov1)^2)
}
results_n200_m5 = matrix(NA, B, 2)
for(b in 1:B){
  data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 200, m = 5, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
  cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
  C_est1 <- cov.obj1$C
  cov_est1 = predict_cov(regular_grid, C_est1)
  dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
  C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
  cov1 = predict_cov(regular_grid, C_1)
  results_n200_m5[b,1] = cov.obj1$error
  results_n200_m5[b,2] = mean((cov_est1 - cov1)^2)
}
results_n200_m10 = matrix(NA, B, 2)
for(b in 1:B){
  data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 200, m = 10, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
  cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
  C_est1 <- cov.obj1$C
  cov_est1 = predict_cov(regular_grid, C_est1)
  dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
  C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
  cov1 = predict_cov(regular_grid, C_1)
  results_n200_m10[b,1] = cov.obj1$error
  results_n200_m10[b,2] = mean((cov_est1 - cov1)^2)
}
results_n200_m15 = matrix(NA, B, 2)
for(b in 1:B){
  data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 200, m = 15, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
  cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
  C_est1 <- cov.obj1$C
  cov_est1 = predict_cov(regular_grid, C_est1)
  dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
  C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
  cov1 = predict_cov(regular_grid, C_1)
  results_n200_m15[b,1] = cov.obj1$error
  results_n200_m15[b,2] = mean((cov_est1 - cov1)^2)
}
results_n200_m20 = matrix(NA, B, 2)
for(b in 1:B){
  data1 = fragment_data1(mu = 0, r = r1, sigma1, n = 200, m = 20, sigma_epsilon = 0.1, domain = c(0, 1), delta = 0.3)
  cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda, ext, maxIt)
  C_est1 <- cov.obj1$C
  cov_est1 = predict_cov(regular_grid, C_est1)
  dist_mat = dist(1:r1, method = "manhattan", diag = T, upper = T)
  C_1 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r1)))
  cov1 = predict_cov(regular_grid, C_1)
  results_n200_m20[b,1] = cov.obj1$error
  results_n200_m20[b,2] = mean((cov_est1 - cov1)^2)
}
boxplot(results_n100_m5[,1]/(100), results_n100_m10[,1]/(100), results_n100_m15[,1]/(100), results_n100_m20[,1]/(100), results_n150_m5[,1]/(150), results_n150_m10[,1]/(150), results_n150_m15[,1]/(150), results_n150_m20[,1]/(150), results_n200_m5[,1]/(200), results_n200_m10[,1]/(200), results_n200_m15[,1]/(200), results_n200_m20[,1]/(200))
boxplot(results_n100_m5[,2]/50, results_n100_m10[,2]/50, results_n100_m15[,2]/50, results_n100_m20[,2]/50, results_n150_m5[,2]/50, results_n150_m10[,2]/50, results_n150_m15[,2]/50, results_n150_m20[,2]/50, results_n200_m5[,2]/50, results_n200_m10[,2]/50, results_n200_m15[,2]/50, results_n200_m20[,2]/50)








B = 200
lambda = 0.00001
results2_n100_m5 = matrix(NA, B, 2)
for(b in 1:B){
  data2 = fragment_data(mu = 0, cov = periodic_cov_fct2, r = r2, sigma2, n = 100, m = 5, sigma_epsilon = 0.5, domain = c(0, 1), delta = 0.3)
  cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
  C_est2 <- cov.obj2$C
  dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
  C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
  cov_est2 = predict_cov(seq(0.05, 0.95, length.out = 50), C_est2, ext)
  results2_n100_m5[b,1] = cov.obj2$error
  results2_n100_m5[b,2] = sum((cov_est2 - cov2)^2)
}
results2_n100_m10 = matrix(NA, B, 2)
for(b in 1:B){
  data2 = fragment_data(mu = 0, cov = periodic_cov_fct2, r = r2, sigma2, n = 100, m = 10, sigma_epsilon = 0.5, domain = c(0, 1), delta = 0.3)
  cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
  C_est2 <- cov.obj2$C
  dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
  C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
  cov_est2 = predict_cov(seq(0.05, 0.95, length.out = 50), C_est2, ext)
  results2_n100_m10[b,1] = cov.obj2$error
  results2_n100_m10[b,2] = sum((cov_est2 - cov2)^2)
}
results2_n100_m15 = matrix(NA, B, 2)
for(b in 1:B){
  data2 = fragment_data(mu = 0, cov = periodic_cov_fct2, r = r2, sigma2, n = 100, m = 15, sigma_epsilon = 0.5, domain = c(0, 1), delta = 0.3)
  cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
  C_est2 <- cov.obj2$C
  dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
  C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
  cov_est2 = predict_cov(seq(0.05, 0.95, length.out = 50), C_est2, ext)
  results2_n100_m15[b,1] = cov.obj2$error
  results2_n100_m15[b,2] = sum((cov_est2 - cov2)^2)
}
results2_n100_m20 = matrix(NA, B, 2)
for(b in 1:B){
  data2 = fragment_data(mu = 0, cov = periodic_cov_fct2, r = r2, sigma2, n = 100, m = 20, sigma_epsilon = 0.5, domain = c(0, 1), delta = 0.3)
  cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
  C_est2 <- cov.obj2$C
  dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
  C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
  cov_est2 = predict_cov(seq(0.05, 0.95, length.out = 50), C_est2, ext)
  results2_n100_m20[b,1] = cov.obj2$error
  results2_n100_m20[b,2] = sum((cov_est2 - cov2)^2)
}
results2_n150_m5 = matrix(NA, B, 2)
for(b in 1:B){
  data2 = fragment_data(mu = 0, cov = periodic_cov_fct2, r = r2, sigma2, n = 150, m = 5, sigma_epsilon = 0.5, domain = c(0, 1), delta = 0.3)
  cov.obj2 = cov_basis(data1$t, data1$y, data1$r, r2, lambda, ext, maxIt)
  C_est2 <- cov.obj2$C
  dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
  C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
  cov_est2 = predict_cov(seq(0.05, 0.95, length.out = 50), C_est2, ext)
  results2_n150_m5[b,1] = cov.obj2$error
  results2_n150_m5[b,2] = sum((cov_est2 - cov2)^2)
}
results2_n150_m10 = matrix(NA, B, 2)
for(b in 1:B){
  data2 = fragment_data(mu = 0, cov = periodic_cov_fct2, r = r2, sigma2, n = 150, m = 10, sigma_epsilon = 0.5, domain = c(0, 1), delta = 0.3)
  cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
  C_est2 <- cov.obj2$C
  dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
  C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
  cov_est2 = predict_cov(seq(0.05, 0.95, length.out = 50), C_est2, ext)
  results2_n150_m10[b,1] = cov.obj2$error
  results2_n150_m10[b,2] = sum((cov_est2 - cov2)^2)
}
results2_n150_m15 = matrix(NA, B, 2)
for(b in 1:B){
  data2 = fragment_data(mu = 0, cov = periodic_cov_fct2, r = r2, sigma2, n = 150, m = 15, sigma_epsilon = 0.5, domain = c(0, 1), delta = 0.3)
  cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
  C_est2 <- cov.obj2$C
  dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
  C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
  cov_est2 = predict_cov(seq(0.05, 0.95, length.out = 50), C_est2, ext)
  results2_n150_m15[b,1] = cov.obj2$error
  results2_n150_m15[b,2] = sum((cov_est2 - cov2)^2)
}
results2_n150_m20 = matrix(NA, B, 2)
for(b in 1:B){
  data2 = fragment_data(mu = 0, cov = periodic_cov_fct2, r = r2, sigma2, n = 150, m = 20, sigma_epsilon = 0.5, domain = c(0, 1), delta = 0.3)
  cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
  C_est2 <- cov.obj2$C
  dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
  C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
  cov_est2 = predict_cov(seq(0.05, 0.95, length.out = 50), C_est2, ext)
  results2_n150_m20[b,1] = cov.obj2$error
  results2_n150_m20[b,2] = sum((cov_est2 - cov2)^2)
}
results2_n200_m5 = matrix(NA, B, 2)
for(b in 1:B){
  data2 = fragment_data(mu = 0, cov = periodic_cov_fct2, r = r2, sigma2, n = 200, m = 5, sigma_epsilon = 0.5, domain = c(0, 1), delta = 0.3)
  cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
  C_est2 <- cov.obj2$C
  dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
  C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
  cov_est2 = predict_cov(seq(0.05, 0.95, length.out = 50), C_est2, ext)
  results2_n200_m5[b,1] = cov.obj2$error
  results2_n200_m5[b,2] = sum((cov_est2 - cov2)^2)
}
results2_n200_m10 = matrix(NA, B, 2)
for(b in 1:B){
  data2 = fragment_data(mu = 0, cov = periodic_cov_fct2, r = r2, sigma2, n = 200, m = 10, sigma_epsilon = 0.5, domain = c(0, 1), delta = 0.3)
  cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
  C_est2 <- cov.obj2$C
  dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
  C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
  cov_est2 = predict_cov(seq(0.05, 0.95, length.out = 50), C_est2, ext)
  results2_n200_m10[b,1] = cov.obj2$error
  results2_n200_m10[b,2] = sum((cov_est2 - cov2)^2)
}
results2_n200_m15 = matrix(NA, B, 2)
for(b in 1:B){
  data2 = fragment_data(mu = 0, cov = periodic_cov_fct2, r = r2, sigma2, n = 200, m = 15, sigma_epsilon = 0.5, domain = c(0, 1), delta = 0.3)
  cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
  C_est2 <- cov.obj2$C
  dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
  C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
  cov_est2 = predict_cov(seq(0.05, 0.95, length.out = 50), C_est2, ext)
  results2_n200_m15[b,1] = cov.obj2$error
  results2_n200_m15[b,2] = sum((cov_est2 - cov2)^2)
}
results2_n200_m20 = matrix(NA, B, 2)
for(b in 1:B){
  data2 = fragment_data(mu = 0, cov = periodic_cov_fct2, r = r2, sigma2, n = 200, m = 20, sigma_epsilon = 0.5, domain = c(0, 1), delta = 0.3)
  cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda, ext, maxIt)
  C_est2 <- cov.obj2$C
  dist_mat = dist(1:r2, method = "manhattan", diag = T, upper = T)
  C_2 <- as.matrix(2^(-dist_mat - 5/2)) + diag(1.5^(1-(1:r2)))
  cov_est2 = predict_cov(seq(0.05, 0.95, length.out = 50), C_est2, ext)
  results2_n200_m20[b,1] = cov.obj2$error
  results2_n200_m20[b,2] = sum((cov_est2 - cov2)^2)
}
boxplot(results2_n100_m5[,1]/(100), results2_n100_m10[,1]/(100), results2_n100_m15[,1]/(100), results2_n100_m20[,1]/(100), results2_n150_m5[,1]/(150), results2_n150_m10[,1]/(150), results2_n150_m15[,1]/(150), results2_n150_m20[,1]/(150), results2_n200_m5[,1]/(200), results2_n200_m10[,1]/(200), results2_n200_m15[,1]/(200), results2_n200_m20[,1]/(200))
boxplot(results2_n100_m5[,2]/50, results2_n100_m10[,2]/50, results2_n100_m15[,2]/50, results2_n100_m20[,2]/50, results2_n150_m5[,2]/50, results2_n150_m10[,2]/50, results2_n150_m15[,2]/50, results2_n150_m20[,2]/50, results2_n200_m5[,2]/50, results2_n200_m10[,2]/50, results2_n200_m15[,2]/50, results2_n200_m20[,2]/50)
