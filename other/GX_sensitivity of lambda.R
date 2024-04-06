# without change point
ext = 0.05
maxIt = 1
r = 3
sigma = 1
lambda <- c(0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001)
iteration <- 100
error_matrix <- matrix(0, nrow = iteration, ncol = length(lambda))
for (j in 1:iteration) {
  data = temp_fragment_data9(mu = 0, r = r, sigma, n = 100, 
                             m = 30, sigma_epsilon = 0.00001,
                             domain = c(0, 1), delta = 0.5)
  for (i in 1:length(lambda)) {
    cov.obj = cov_basis(data$t, data$y, data$r, r1, lambda[i], ext, maxIt)
    error_matrix[j,i] <- cov.obj$error
  }
}

error_matrix

boxplot(error_matrix[,1], error_matrix[,2], error_matrix[,3], error_matrix[,4],
        error_matrix[,5], error_matrix[,6], error_matrix[,7], error_matrix[,8])

# boxplot(error_matrix[,1], error_matrix[,2], error_matrix[,3], error_matrix[,4],
#         error_matrix[,5], error_matrix[,6], error_matrix[,7],
#         ylim = c(100000, 500000))




# with one change point
ext = 0.05
maxIt = 1
r1 = 3
r2 = 3
r = 3
sigma1 = 1
sigma2 = 1
lambda <- c(0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001)
iteration <- 100

error_matrix1 <- matrix(0, nrow = iteration, ncol = length(lambda))
error_matrix2 <- matrix(0, nrow = iteration, ncol = length(lambda))
error_matrix_whole <- matrix(0, nrow = iteration, ncol = length(lambda))
error_matrix_diff <- matrix(0, nrow = iteration, ncol = length(lambda))


for (j in 1:iteration) {
  data1 = temp_fragment_data7(mu = 0, r = r1, sigma1, n = 100, m = 30, sigma_epsilon = 0.00001, domain = c(0, 1), delta = 0.5)
  data2 = temp_fragment_data9(mu = 0, r = r2, sigma2, n = 100, m = 30, sigma_epsilon = 0.00001, domain = c(0, 1), delta = 0.5)
  data = list("t"= rbind(data1$t, data2$t), "y" = rbind(data1$y, data2$y), "r" = cbind(data1$r, data2$r))
  for (i in 1:length(lambda)) {
    cov.obj1 = cov_basis(data1$t, data1$y, data1$r, r1, lambda[i], ext, maxIt)
    cov.obj2 = cov_basis(data2$t, data2$y, data2$r, r2, lambda[i], ext, maxIt)
    cov.obj_whole = cov_basis(data$t, data$y, data$r, r, lambda[i], ext, maxIt)
    
    error_matrix1[j,i] <- cov.obj1$error
    error_matrix2[j,i] <- cov.obj2$error
    error_matrix_whole[j,i] <- cov.obj_whole$error
    error_matrix_diff[j,i] <- cov.obj_whole$error - cov.obj1$error - cov.obj2$error
  }
}


boxplot(error_matrix1[,1], error_matrix1[,2], error_matrix1[,3], error_matrix1[,4],
        error_matrix1[,5], error_matrix1[,6], error_matrix1[,7], error_matrix1[,8])

boxplot(error_matrix2[,1], error_matrix2[,2], error_matrix2[,3], error_matrix2[,4],
        error_matrix2[,5], error_matrix2[,6], error_matrix2[,7], error_matrix2[,8])

boxplot(error_matrix_whole[,1], error_matrix_whole[,2], error_matrix_whole[,3], error_matrix_whole[,4],
        error_matrix_whole[,5], error_matrix_whole[,6], error_matrix_whole[,7], error_matrix_whole[,8])

boxplot(error_matrix_diff[,1], error_matrix_diff[,2], error_matrix_diff[,3], error_matrix_diff[,4],
        error_matrix_diff[,5], error_matrix_diff[,6], error_matrix_diff[,7], error_matrix_diff[,8])























