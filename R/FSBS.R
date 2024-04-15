##Gaussian Kernel
#' @export
K_h<-function(x,h,d){
  res<-(1/(h)^d)*(1/sqrt(2*pi)^d)*exp((-1/2)*(norm(x/h, type = "2"))^2)
  return(res)
}


## Seeded intervals
#' @export
seeded.intervals <- function(Tb,C, decay = 2, unique.int = T){
  aux <- log(Tb)
  depth <- C*log(aux, base = decay)
  depth <- ceiling(depth)
  M <- sum(2^(1:depth)-1)
  
  boundary_mtx <- matrix(NA, ncol = 3)
  colnames(boundary_mtx) <- c("layer","st", "end")
  boundary_mtx[1, ] <- c(1,0, Tb) #first interval when k=1 in the definition
  
  if(depth>=2)
  {
    for(i in 2:depth){
      int_length <- Tb * (1/decay)^(i-1) #my lk in def seeded intervals
      n_int <- ceiling(Tb/int_length)*2-1 #my nk bar in def seeded intervals sometimes very slight numerical inaccuracies
      aux=floor(seq(0, Tb-int_length, length.out = (n_int)))
      lab=matrix(i,nrow=length(aux),ncol=1)
      aux=cbind(lab,aux)
      boundary_mtx <- rbind(boundary_mtx,cbind(aux, ceiling(seq(int_length, Tb, length.out = (n_int)))))
    }
  }
  
  if(unique.int){return(unique(boundary_mtx))}
  boundary_mtx
}

### dat is the data matrix. sum(n_t) observations, 1+d+n columns, first column= time label, 2:d+1 col= x_t,i vector, d+2:1+d+n col= y_t,i vec
#x is a d-dim vector

#' @export
p_hat<-function(x, dat,h_bar,d){
  N_t<-nrow(dat)
  aux<-matrix(NaN,nrow=N_t, ncol=1)
  for(i in 1:N_t){
    aux[i]=K_h(x-dat[i,2:(d+1)],h_bar,d)
  }
  aux<-as.vector(aux)
  res<-mean(aux)
  return(res)
}

#Statistic F_t,h
#' @export
statistic<-function(x,dat,t,h,h_bar,d,phat){
  indices=which(dat[,1]==t)
  n_t=length(indices)
  X_t=dat[indices, 2:(d+1)]
  Y_t=dat[indices,(d+2):(1+d+1) ]
  aux<-matrix(NaN,nrow=n_t, ncol=1)
  for (i in 1:n_t) {
    aux[i]=Y_t[i]*K_h(x-X_t[i],h,d)
  }
  res<-mean(aux)/phat
  return(res)
}


#CUSUM statistic
#' @export
CUSUM<-function(s,e,x,t,h,h_bar,dat,d,phat){
  aux_1=0
  aux=s+1
  if(aux<=t){
    for(i in aux:t){
      aux_1=aux_1+statistic(x,dat,i,h,h_bar,d,phat)
    }
  }
  aux_1=aux_1 * sqrt((e-t)/((e-s)*(t-s)))
  aux_2=0
  aux1=t+1
  if(aux1<=e)
  {
    for(i in aux1:e){
      aux_2=aux_2+statistic(x,dat,i,h,h_bar,d,phat)
    }
  }
  aux_2=aux_2 *sqrt((t-s)/((e-s)*(e-t)))
  return(aux_1-aux_2)
}

# Algorithm
#' @export
FSBS<-function(s,e,h,h_bar,tau,rho,dat,s.inter,samp,d,CUSUM_mat){
  n_inter=length(s.inter[,1])
  samp_size=length(samp)
  AD=matrix(ncol=2)
  AD[1,1]=-1
  AD[1,2]=0
  for(i in 1:n_inter){
    alpha=s.inter[i,2]
    beta=s.inter[i,3]
    alp_bet=beta-alpha
    if((s<= alpha) && (beta<=e) && (alp_bet > 2*rho) && (ceiling(alpha+rho) <= floor(beta-rho) )  ){
      indices1=which(CUSUM_mat[,1]==alpha)
      mat_aux=CUSUM_mat[indices1,]
      indices2=which(mat_aux[,2]==beta)
      mat_aux=mat_aux[indices2,]
      A_aux=max(mat_aux[,5])
      D_aux=mat_aux[which(mat_aux[,5]==A_aux)[1],4]
    }
    else{
      A_aux=-1
      D_aux=0
    }
    AD=rbind(AD,c(A_aux,D_aux))
  }
  
  
  A_star=max(AD[,1])
  ind=which(AD[,1]==A_star)[1] #First element that reaches the maximum
  if (is.na(A_star)) {
    print('Missing')
    return(NULL)
  }
  else{
    if(A_star > tau){
      S=AD[ind,2]
      aux=AD[ind,2]
      S=append(S,FSBS(s,aux,h,h_bar,tau, rho, dat,s.inter,samp,d,CUSUM_mat))
      S=append(S,FSBS(aux,e,h,h_bar,tau, rho, dat,s.inter,samp,d,CUSUM_mat))
      return(S)
    }
  }
}


#s.inter= seeded intervals. Has three columns. The first one is a label of the layer the interval belongs to
#d= dim of x_t,i's
#S empty vector to store changepoints

#' @export
seedBS <- function (dat, s.inter, h, h_bar, tau, d,n ){
  #Obtain params from dat
  S=0
  S=as.vector(S)
  Tb=max(dat[,1]) #number of times
  N_t=length(dat[,1]) #total number of observations
  
  #Obtain params from s.inter
  n_inter=length(s.inter[,1]) #total number of seeded intervals
  
  #Sample log(Tb) points
  samp_size=floor(log(Tb))
  set.seed(1)
  samp=sample(N_t, samp_size, replace=F)
  
  #initialization of interval
  s=0
  e=Tb
  
  #set rho
  rho=log(Tb)/(n*h^{d})
  
  #Creating statistic matrix
  statistic_matrix=matrix(NaN,nrow=Tb,ncol=samp_size)
  for(i in 1:Tb){
    for(j in 1:samp_size)
    {
      x=dat[samp[j],2:(d+1)]
      phat=p_hat(x,dat,h_bar,d)
      statistic_matrix[i,j]=statistic(x,dat,i,h,h_bar,d,phat)
    }
  }
  
  
  
  #Compute CUSUM statistics 
  n_inter=length(s.inter[,1])
  CUSUM_mat=matrix(NaN, ncol=6)
  aux_count=0
  for(i in 1:n_inter){
    alpha=s.inter[i,2]
    beta=s.inter[i,3]
    if((beta-alpha) > (2*rho)){
      t.lo=ceiling(alpha+rho)
      t.up=floor(beta-rho)
      if(t.lo <= t.up){
        for(j in 1:samp_size)
        {
          for(t in t.lo:t.up)
          {
            aux_count=1+aux_count
            
            val_from_statis_matrix_1=statistic_matrix[(alpha+1):t,j]
            sum_val_from_statis_matrix_1=sum(val_from_statis_matrix_1)
            sum_val_from_statis_matrix_1=sum_val_from_statis_matrix_1*sqrt((beta-t)/((beta-alpha)*(t-alpha)))
            val_from_statis_matrix_2=statistic_matrix[(t+1):beta,j]
            sum_val_from_statis_matrix_2=sum(val_from_statis_matrix_2)
            sum_val_from_statis_matrix_2=sum_val_from_statis_matrix_2*sqrt((t-alpha)/((beta-alpha)*(beta-t)))
            CUSUM_x_alpha_beta_t=sum_val_from_statis_matrix_1-sum_val_from_statis_matrix_2
            CUSUM_mat=rbind(CUSUM_mat, c( alpha, beta, j, t , abs(CUSUM_x_alpha_beta_t),aux_count ) )
          }
          
        }
        
      }
      
    }
  }
  CUSUM_mat=CUSUM_mat[-1,]
  
  S=FSBS(s,e,h,h_bar,tau,rho,dat,s.inter,samp,d,CUSUM_mat)
  return(S)
}




