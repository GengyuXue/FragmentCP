devtools::install_github("linulysses/mcfda")
library(dplyr)
library(mcfda)
require(lubridate)
require(ggplot2)
library(gridExtra)
library(FragmentCP)
library(changepoints)
library(data.table)
library(MASS)
library(fda)
library(tictoc)
library(rainbow)
library(sde)
library(xtable)
library(mvtnorm)
library(tseries)
library(expm)
library(tensorA)

generate_pairs <- function(data, colnum){
  pairs <- matrix(NaN, nrow = 1, ncol = 2)
  for (i in 1:12) {
    for (j in 1:31) {
      temp_demand <- data %>%
        filter(month_index == i & day_index == j) %>%
        select(colnum)
      
      if(nrow(temp_demand) == 0){
        next
      }
      
      else{
        temp_pairs <- as.matrix(expand.grid(as.vector(temp_demand[,1]), as.vector(temp_demand[,1])))
        pairs <- rbind(pairs, temp_pairs)
      }
    }
  }
  pairs <- pairs[-1,]
  return(pairs)
}

# Load electricity price
whole_price <- read.csv("Germany_price.csv")

whole_price$Datetime..UTC. <- ymd_hms(whole_price$Datetime..UTC.)

whole_price$year_index <- year(whole_price$Datetime..UTC.)
whole_price$month_index <- month(whole_price$Datetime..UTC.)
whole_price$day_index <- day(whole_price$Datetime..UTC.)
whole_price$hour_index <- hour(whole_price$Datetime..UTC.)


# slect the price for year 2022
price_2022 <- whole_price %>% 
  filter(year_index == 2022) %>%
  select(-Datetime..UTC., -Datetime..Local., -ISO3.Code, -Country, -year_index) %>%
  rename(price = Price..EUR.MWhe.)

# load hourly demand
whole_demand_2022 <- read.csv("monthly_hourly_load_values_2022 .csv", sep = ";")

GE_demand_2022 <- whole_demand_2022 %>%
  filter(CountryCode == "GE") %>%
  select(-CreateDate, -Cov_ratio, -CountryCode)


GE_demand_2022[2033,2] <- "27/03/2022 02:00"

GE_demand_2022$DateUTC <- dmy_hm(GE_demand_2022$DateUTC)

GE_demand_2022$year_index <- year(GE_demand_2022$DateUTC)
GE_demand_2022$month_index <- month(GE_demand_2022$DateUTC)
GE_demand_2022$day_index <- day(GE_demand_2022$DateUTC)
GE_demand_2022$hour_index <- hour(GE_demand_2022$DateUTC)

GE_demand_2022 <- GE_demand_2022 %>%
  select(Value_ScaleTo100, month_index, day_index, hour_index) %>%
  rename(demand = Value_ScaleTo100)


#combine two dataset and remove rows with missing values
data_2022 <- left_join(price_2022, GE_demand_2022, by = c("month_index", "day_index", "hour_index"))
rows_with_na <- !complete.cases(data_2022)
rows_with_na_df <- data_2022[rows_with_na, ]

index_with_missing <- c(which(data_2022$month_index == 1 & data_2022$day_index == 17),
                        which(data_2022$month_index == 1 & data_2022$day_index == 18),
                        which(data_2022$month_index == 1 & data_2022$day_index == 19),
                        which(data_2022$month_index == 1 & data_2022$day_index == 31),
                        which(data_2022$month_index == 3 & data_2022$day_index == 5),
                        which(data_2022$month_index == 4 & data_2022$day_index == 12),
                        which(data_2022$month_index == 4 & data_2022$day_index == 25),
                        which(data_2022$month_index == 6 & data_2022$day_index == 2),
                        which(data_2022$month_index == 6 & data_2022$day_index == 21),
                        which(data_2022$month_index == 7 & data_2022$day_index == 19))
data_2022 <- data_2022[-index_with_missing,]

# peak time data
peak_data <- data_2022 %>%
  filter(hour_index<=19 & hour_index >= 8) 

# normalize demand to [0,1]
min_demand <- min(peak_data$demand)
max_demand <- max(peak_data$demand)
n_demand <- sapply(peak_data$demand, function(x) {(x-min_demand)/(max_demand- min_demand)})
peak_data <- cbind(peak_data, n_demand)

min_demand <- min(data_2022$demand)
max_demand <- max(data_2022$demand)
n_demand <- sapply(data_2022$demand, function(x) {(x-min_demand)/(max_demand- min_demand)})
data_2022 <- cbind(data_2022, n_demand)

########################## plot of domain of covariance function ##########################
pairs_whole <- generate_pairs(data_2022, 6)
df_pairs_whole <- data.frame("var1"= pairs_whole[,1], "var2" = pairs_whole[,2])

peak_pairs <- generate_pairs(peak_data, colnum = 6)
df_peak_pairs <- data.frame("var1"= peak_pairs[,1], "var2" = peak_pairs[,2])
df_diff_temp <- anti_join(df_pairs_whole, df_peak_pairs)

df_diff <- cbind(df_diff_temp, "var3" = rep("whole", nrow(df_diff_temp)))
df_pairs_peak <-  cbind(df_peak_pairs, "var3" = rep("peak", nrow(df_peak_pairs)))

df_plot_domain <- rbind(df_diff, df_pairs_peak)


gg1 <- ggplot(df_plot_domain, aes(x= var1, y=var2, col = var3))+
  geom_point()+ xlab("Electr. Demand") + ylab("Electr. Demand") + theme_classic()+
  scale_color_manual(values=c("grey40", "grey"))+
  theme(legend.position = c(0.9, 0.1))+
  labs(color = "Data type")
plot(gg1)

########################## mean estimation for data_2022 ##########################
#change the form of data
Ly <- list()
Lt <- list()
index <- 1
for (i in 1:12) {
  for (j in 1:31) {
    temp_demand <- data_2022 %>%
      filter(month_index == i & day_index == j) %>%
      arrange(n_demand)
    
    if(nrow(temp_demand) == 0){
      next
    }
    
    else{
      Ly[[index]] <- temp_demand$price
      Lt[[index]] <- temp_demand$n_demand
    }
    index <- index+1
  }
}


D <- list(t = Lt, y=Ly)

#mean estimation
mu.obj <- meanfunc(D$t,D$y,newt=NULL,method='FOURIER', tuning='cv',weig=NULL, domain=c(0,1))

Ly <- matrix(NaN, 355, 24)
Lt <- matrix(NaN, 355, 24)
index <- 1
for (i in 1:12) {
  for (j in 1:31) {
    temp_demand <- data_2022 %>%
      filter(month_index == i & day_index == j) %>%
      arrange(n_demand)
    
    if(nrow(temp_demand) == 0){
      next
    }
    
    else{
      for (k in 1:nrow(temp_demand)) {
        Ly[index, k] <- temp_demand$price[k]
        Lt[index, k] <- temp_demand$n_demand[k]
      }
      index <- index+1
    }
  }
}

rm_Ly <- matrix(NaN, 355, 24)
for (i in 1:355) {
  rm_Ly[i,] <- Ly[i,] - predict(mu.obj, Lt[i,])
}

Lr <- NULL
for(i in 1:355){
  Lr = cbind(Lr, tcrossprod(rm_Ly[i,]))
}

list_data_2022 <- list(t = Lt, y = rm_Ly, r = Lr)


########################## mean estimation for peak_data ##########################
#change the form of data
Ly <- list()
Lt <- list()
index <- 1
for (i in 1:12) {
  for (j in 1:31) {
    temp_demand <- peak_data %>%
      filter(month_index == i & day_index == j) %>%
      arrange(n_demand)
    
    if(nrow(temp_demand) == 0){
      next
    }
    
    else{
      Ly[[index]] <- temp_demand$price
      Lt[[index]] <- temp_demand$n_demand
    }
    index <- index+1
  }
}


D <- list(t = Lt, y=Ly)

#mean estimation
mu.obj <- meanfunc(D$t,D$y,newt=NULL,method='FOURIER', tuning='cv',weig=NULL, domain=c(0,1))

Ly <- matrix(NaN, 355, 12)
Lt <- matrix(NaN, 355, 12)
index <- 1
for (i in 1:12) {
  for (j in 1:31) {
    temp_demand <- peak_data %>%
      filter(month_index == i & day_index == j) %>%
      arrange(n_demand)
    
    if(nrow(temp_demand) == 0){
      next
    }
    
    else{
      for (k in 1:nrow(temp_demand)) {
        Ly[index, k] <- temp_demand$price[k]
        Lt[index, k] <- temp_demand$n_demand[k]
      }
      index <- index+1
    }
  }
}

rm_Ly <- matrix(NaN, 355, 12)
for (i in 1:355) {
  rm_Ly[i,] <- Ly[i,] - predict(mu.obj, Lt[i,])
}

Lr <- NULL
for(i in 1:355){
  Lr = cbind(Lr, tcrossprod(rm_Ly[i,]))
}

list_peak_data <- list(t = Lt, y = rm_Ly, r = Lr)


################################### Change point with list_data_2022 ##############################################
# FFDP
ext=0.1
maxIt=1
xi_set = c(60000000)
CV_cpt_result = CV_search_DP_fragment(list_data_2022$t, list_data_2022$y, list_data_2022$r, r = 2, lambda=1e-09,
                                      xi_set, ext, maxIt)
min_idx = which.min(CV_cpt_result$test_error) 
xi_set[min_idx]
cpt_init = unlist(CV_cpt_result$cpt_hat[min_idx])
cpt_init
# 51 197 225

# SBS
ext=0.1
maxIt=1
zeta_set = c(100000,150000,200000)
CV_SBS_result = CV_search_SBS_fragment(list_data_2022$t, list_data_2022$y, list_data_2022$r, r=3, lambda, 
                                       ext, maxIt, zeta_set, Delta = 10)
min_idx = which.min(CV_SBS_result$test_error) 
zeta_set[min_idx]
cpt_init = unlist(CV_SBS_result$cpt_hat[min_idx])

# WCUSUM
kappa = 0.22
lens = 14

smooth_function <- function(data, num_basis = 5, rangeval = c(0,1), no_grids = 30){
  # change format of data
  m = ncol(data$y)
  Tb = nrow(data$y)
  nts<-matrix(NaN,nrow=Tb, ncol=1)
  for(i in 1:Tb){
    nts[i] = m
  }
  N_t=sum(nts)
  # N_t
  
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
  
  # create the data matrix containing t, y and labels
  dat[, 1] <- labels_Ts
  dat[, 2] <- t_vector
  dat[, 3] <- y_vector
  
  # create data used for smoothing
  dat_y=matrix(0,m,Tb)
  aux_count=0
  for(t in 1:Tb){
    aux_count=m+aux_count
    for(i in 1:m){
      dat_y[i,t]=y_vector[aux_count-m+i]
    }
  }
  # dat_y
  
  dat_t=matrix(0,m,Tb)
  aux_count=0
  for(t in 1:Tb){
    aux_count=m+aux_count
    for(i in 1:m){
      dat_t[i,t]=t_vector[aux_count-m+i]
    }
  }
  # dat_t
  
  # smoothing using spline basis
  basis <- create.bspline.basis(rangeval=c(0, 1), nbasis=num_basis, norder=4)
  smooth_func <- smooth.basis(dat_t, dat_y, basis)$fd
  eval_points <- seq(0, 1, length.out = no_grids)
  eval <- as.matrix(eval.fd(eval_points , smooth_func))
  return(eval)
}
data_2022_smooth <- smooth_function(list_data_2022, num_basis = 4, rangeval = c(0,1), no_grids = 50)
stat_d0 = weight_TNstat(data_2022_smooth, kappa)
cv_d0 = weight_criticalvalueMC(data_2022_smooth,len = lens,kappa)
# print(stat_d0[[2]])
if (stat_d0[[1]]> cv_d0[1]){
  change_est = stat_d0[[2]]
}else{
  change_est= 0
}
change_est
#106


################################# Change point with list_peak_data ###################################
# FFDP
ext=0.1
maxIt=1
xi_set = c(200000, 300000)
CV_cpt_result = CV_search_DP_fragment(list_peak_data$t, list_peak_data$y, list_peak_data$r, r = 3, lambda=1e-09,
                                      xi_set, ext, maxIt)
min_idx = which.min(CV_cpt_result$test_error) 
xi_set[min_idx]
cpt_init = unlist(CV_cpt_result$cpt_hat[min_idx])
cpt_init
# 101, 197

# SBS
ext=0.1
maxIt=1
zeta_set = c(400000, 500000)
CV_SBS_result = CV_search_SBS_fragment(list_peak_data$t, list_peak_data$y, list_peak_data$r, r=3, lambda, 
                                       ext, maxIt, zeta_set, Delta = 10)
min_idx = which.min(CV_SBS_result$test_error) 
zeta_set[min_idx]
cpt_init = unlist(CV_SBS_result$cpt_hat[min_idx])
cpt_init
# 247 215 111 175

# WCUSUM
kappa = 0.2
lens = 15

peak_data_smooth <- smooth_function(list_peak_data, num_basis = 4, rangeval = c(0,1), no_grids = 50)
stat_d0 = weight_TNstat(peak_data_smooth, kappa)
cv_d0 = weight_criticalvalueMC(peak_data_smooth,len = lens,kappa)
# print(stat_d0[[2]])
if (stat_d0[[1]]> cv_d0[1]){
  change_est = stat_d0[[2]]
}else{
  change_est= 0
}
change_est
# 218


















