ext = 0.05
maxIt = 1
r = 3
sigma = 1
Delta = 50
regular_grid = seq(0.05, 0.95, length.out = 50)



lambda <- c(0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001)
iteration <- 100
error_matrix <- matrix(0, nrow = iteration, ncol = length(lambda))
for (j in 1:iteration) {
  data = temp_fragment_data1(mu = 0, r = r, sigma, n = 100, 
                             m = 30, sigma_epsilon = 0.00001,
                             domain = c(0, 1), delta = 0.6)
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








