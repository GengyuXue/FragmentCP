lambda = 0.00001
ext = 0.1
maxIt = 1
r1 = 3
r2 = 3
sigma1 = 1
sigma2 = 1

# parameter used to generate seeded interval
C=3
#dimension of data
d=1

data1 = temp_fragment_data7(mu = 0, r = r1, sigma1, n = 100, m = 30, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.6)
data2 = temp_fragment_data9(mu = 0, r = r2, sigma2, n = 100, m = 30, sigma_epsilon = 0.01, domain = c(0, 1), delta = 0.6)
data = list("t"= rbind(data1$t, data2$t), "y" = rbind(data1$y, data2$y), "r" = cbind(data1$r, data2$r))

m = ncol(data$y)
Tb = nrow(data$y)
nts<-matrix(NaN,nrow=Tb, ncol=1)
for(i in 1:Tb){
  nts[i] = m
}
N_t=sum(nts)
N_t

dat <- matrix(NaN, nrow = N_t, ncol = 3)

labels_Ts=matrix(NaN,nrow = N_t,ncol=1)
aux1=0
for (i in 1:Tb) {
  aux1=nts[i]+aux1
  for (j in 1:nts[i]) {
    labels_Ts[aux1-nts[i]+j]=i
  }
}
# labels_Ts

t_vector <- matrix(NaN, nrow = N_t, ncol =1)
y_vector <- matrix(NaN, nrow = N_t, ncol =1)
for (i in 1:Tb) {
  for (j in 1:m) {
    t_vector[(i-1)*m+j,1] <- c(data$t[i,j])
    y_vector[(i-1)*m+j,1] <- c(data$y[i,j])
  }
}

# create the data matrix containing x_t_is, y_t_is and labels
dat[, 1] <- labels_Ts
dat[, 2] <- t_vector
dat[, 3] <- y_vector
dat

#The choice of bandwidth h_bar
H_bar_stimated=hpi(dat[,2])
H_bar_stimated
h_bar=H_bar_stimated

######################################################
######################### CV #########################
######################################################
#Creating testing and training sets.
dat_train=dat[(m+1):(2*m),]
dat_test=dat[1:m,]
for(i in 1:(Tb/2-1)){
  dat_test=rbind(dat_test,dat[(2*i*m+1):((2*i+1)*m),])
  dat_train=rbind(dat_train,dat[((2*i+1)*m+1):((2*(i+1))*m),])
}
dat_train[,1]=ceiling(dat_train[,1]/2)

#We create the seeded intervals to work with training data
s.inter_train<-seeded.intervals(Tb/2, C)



#We create the estimator functions 
g_hat_i=function(eta1,eta2,h,x,dat_train,h_bar,phat)
{
  res=0
  for (i in (eta1+1):eta2) {
    res=res+statistic(x,dat_train,i,h,h_bar,d,phat)
  }
  res=res/(eta2-eta1)
  return(res)
}

#We create the possible values for h and tau
h_int= 0.2
# h_int=seq(0.79,0.81, 0.01)
tau_int=c(10,11,12)
# tau_int=c(1.08)
l_tau_int=length(tau_int)

#We compute errors of estimation

errors=matrix(NaN,1,l_tau_int)
K_cv <- matrix(NaN,1,l_tau_int)
CP_cv <- list()
#start_time <- Sys.time()
  for (ind in 1:l_tau_int) {
    S=seedBS(dat_train, s.inter_train, h_int, h_bar, tau_int[ind], d, m)
    K_cv[ind] <- length(S)
    # if(is.null(S)){
    #   errors[ind1,ind2]=10000000000
    #   next
    # }
    # else{
      S=as.vector(S)
      S=append(S,0)
      S=append(S,Tb/2)
      S=sort(S)
      CP_cv[[ind]] <- S
      error=0
      for(j in 1:(length(S)-1)){
        for(t in (S[j]+1):S[j+1]){
          for(i in 1:m){
            indices=which(dat_test[,1]==(2*t-1))
            X_t=dat_test[indices,2:(d+1)]
            Y_t=dat_test[indices,(d+2):(d+2)]
            phat=p_hat(X_t[i],dat_train,h_bar,d)
            error=error+(g_hat_i(S[j],S[j+1],h_int[ind1],X_t[i],dat_train,h_bar,phat)-Y_t[i])^2
          }
        }
      }
      errors[ind1,ind2]=error
    # }  
  }
min_error=min(errors)
tau_min=which(errors==min_error)
tau_min=tau_int[tau_min]




S=seedBS(dat_train, s.inter_train, 0.01, 5, tau_int[ind2], d, m)









