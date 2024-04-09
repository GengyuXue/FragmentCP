# function to generate functional data observed on fragments
#' @export
CV_fragment = function(Lt, Ly, Lr, r, lambda, xi, ext, maxIt){
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
  train_result = DP_fragment(train_Lt, train_Ly, train_Lr, r, lambda, xi, ext, maxIt, Delta=0)
  init_cpt_train = part2local(train_result$partition)
  init_cpt_train.long = c(0, init_cpt_train, nrow(train_Lt))
  diff_point = diff(init_cpt_train.long)
  init_cpt_update = 0
  flag = 1
  for(i in 1:length(diff_point)){
    if(flag == 1){
      gap = diff_point[i]
    }else{
      gap = gap + diff_point[i]
    }
    if(gap >= 10){
      cpt_temp = init_cpt_update[length(init_cpt_update)] + gap
      init_cpt_update = c(init_cpt_update, cpt_temp)
      flag = 1
    }else{
      flag = 0
    }
  }
  init_cpt = odd_indexes[init_cpt_update[c(-1, -length(init_cpt_update))]]
  len = length(init_cpt)
  init_cpt_long = c(init_cpt_update[c(-1, -length(init_cpt_update))], floor(n/2))
  interval = matrix(0, nrow = len+1, ncol = 2)
  interval[1,] = c(1, init_cpt_long[1])
  if(len > 0){
    for(j in 2:(1+len)){
      interval[j,] = c(init_cpt_long[j-1]+1, init_cpt_long[j])
    }
  }
  trainmat = sapply(1:(len+1), function(index) error_seg_fragment(train_Lt, train_Ly, train_Lr, r, interval[index,1], interval[index,2], lambda, ext, maxIt))
  C_mat = vector("list", len+1)
  training_loss = matrix(0, nrow = 1, ncol = len+1)
  for(col in 1:(len+1)){
    C_mat[[col]] = trainmat[1,col]$C
    training_loss[,col] = as.numeric(trainmat[2,col]$error)
  }
  validationmat = sapply(1:(len+1), function(index) error_test_fragment(validation_Lt, validation_Lr, interval[index,1], interval[index,2], C_mat[[index]]))
  result = list(cpt_hat = init_cpt, K_hat = len, test_error = sum(validationmat), train_error = sum(training_loss))
  return(result)
}



#' @noRd
error_test_fragment = function(Lt, Lr, s, e, C_mat){
  m = ncol(Lt)
  Lr_new = Lr[,as.vector(sapply(s:e, function(i){(m*(i-1)+1):(m*i)}))]
  Lp_new = matrix(NA, m, (e-s+1)*m)
  for(i in s:e){
    basis_mat = evaluate_basis(r, c(0,1), Lt[i,])
    Lp_new[,((i-s)*m+1):((i-s+1)*m)] = basis_mat %*% C_mat %*% t(basis_mat)
  }
  norm = sum((Lr_new - Lp_new)^2)
  error = norm/(m*m)
  return(error)
} 


#' @export
CV_search_DP_fragment = function(Lt, Ly, Lr, r, lambda, xi_set, ext, maxIt){
  output = sapply(1:length(xi_set), function(j) CV_fragment(Lt, Ly, Lr, r, lambda, xi_set[j], ext, maxIt))
  cpt_hat = output[seq(1,4,4),]## estimated change points
  K_hat = output[seq(2,4,4),]## number of estimated change points
  test_error = output[seq(3,4,4),]## validation loss
  train_error = output[seq(4,4,4),]## training loss                                                      
  result = list(cpt_hat = cpt_hat, K_hat = K_hat, test_error = test_error, train_error = train_error)
  return(result)
}