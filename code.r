
# Start from scratch
rm(list=ls(all=TRUE)) 

setwd('C:/Users/Ioanna/Documents/MSc in BUSINESS ECONOMICS WITH ANALYTICS/B SEMESTER/Econometrics - Ioannis_Dendramis/3rd')
getwd()
graphics.off()

# we are going to use the packages : expm, urca

library("expm")
library("urca")

# we are going to use the functions: lag_select, irf_fun, ols_est

source("ols_est.R")
source("irf_fun.R")
source("lag_select.R")

data<- read.csv2("datafred_raw.csv",header = T) # Load the data set

head(data,10) # view the first 10 rows of the data
tail(data,10) # view the last 10 rows of the data

str(data) # structure of data
summary(data) # summary of data. Gives descriptive statistics

dim(data) # gives the dimensions of object data

y<-data$GDPC96
p<-data$GDPCTPI
r<-data$FEDFUNDS

ygrowth <- rep(1,252)
pgrowth<- rep(1,252)

for (i in 2:252){
  ygrowth[i]<-(y[i]-y[i-1])/y[i-1]*100
  pgrowth[i]<-(p[i]-p[i-1])/p[i-1]*100
}

ygrowth<-ygrowth[2:252]
pgrowth<-pgrowth[2:252]
r<-r[2:252]

# Check stationarity of the time series variables. We want to establish that the
# variables involved in the analysis are stationary. To this end, plot the
# series, inspect their corellogram and conduct Dickey Fuller tests. Dickey Fuller
# tests are conducted with the functions ur.df from the urca package 

plot(ygrowth, type="l") # plot the series
acf(ygrowth) # Generate the corellogram

df_ygrowth_trend<-ur.df(ygrowth,type = "trend",selectlags = "AIC")
summary(df_ygrowth_trend)

# The test indicates that the series for ygrowth contains no unit roots.

plot(pgrowth, type="l")
acf(pgrowth)

df_pgrowth_trend<-ur.df(pgrowth,type="trend",selectlags = "AIC")
summary(df_pgrowth_trend)

df_pgrowth_drift<-ur.df(pgrowth,type = "drift", selectlags = "AIC")
summary(df_pgrowth_drift)

df_pgrowth_none<-ur.df(pgrowth, type = "none", selectlags = "AIC")
summary(df_pgrowth_none)

# The test indicates that the series for pgrowth contains a unit root.

plot(r, type="l")
acf(r)

df_r_trend<-ur.df(r,type = "trend", selectlags = "AIC")
summary(df_r_trend)

df_r_drift<-ur.df(r,type = "drift", selectlags = "AIC")
summary(df_r_drift)

df_r_none<-ur.df(r,type = "none", selectlags = "AIC")
summary(df_r_none)

# The test indicates that the series for r contains a unit root.

# Therefore, we need to difference the series for pgrowth and r and check
# whether or not the differenced series are stationary

ygrowth<-ygrowth[2:251]
pgrowth<-tail(pgrowth,250) - head(pgrowth,250)
r<-tail(r,250) - head(r,250)

# We are rechecking for stationarity of the differenced series

plot(pgrowth, type="l")
acf(pgrowth)

df_pgrowth_trend<-ur.df(pgrowth, type = "trend", selectlags = "AIC")
summary(df_pgrowth_trend)

# The test indicates that the differenced series for pgrowth contains no unit roots

plot(r, type="l")
acf(r)

df_r_trend<-ur.df(r,type = "trend", selectlags = "AIC")
summary(df_r_trend)

# The test indicates that the differensed series for r contains no unit roots

# All series are stationary. Therefore, we can proceed with the VAR analysis

data<-cbind(ygrowth,pgrowth,r) # This is the data set that we will work with.
                               # Observe that the first variable is ygrowht.
                               # the second is pgrowth and the third is r.

colnames(data)<-NULL #erase column names

n<- dim(data)[1];n # number of observasions in my data set
k<- dim(data)[2];k # number of variables in my data set

#Next, we determine the optimal number of lags. we are using the function lag_select
#and we are checking for up to 12 lag

optim_lag<-lag_select(12)
optim_lag # Both the SBIC and the HQIC suggest the we should set the number of
          # lags equal to 2
p<-2

# In this part of the code we calculate and plot the Impulse Response Functions

hor<-10 # The time horizon for the irf's

results<-ols_est(data,p,n)
bhat<-results$bhat;bhat
sigmahat<-results$sigmahat;sigmahat
resida<-results$resids

H<-t(chol(sigmahat));H # indentification of matrix H using Cholesky decomposition.
                       # H is the impact multiplier matrix.

bhatt<-bhat[2:(k*p+1),]
companion<-cbind(diag(k*(p-1)),matrix(0,nrow =(k*(p-1)),ncol = k))
companion<-rbind(t(bhatt),companion);companion 

irf<-irf_fun(companion,k,p,H,hor)

irf11<-irf[1,]
irf21<-irf[2,]
irf31<-irf[3,]
irf12<-irf[4,]
irf22<-irf[5,]
irf32<-irf[6,]
irf13<-irf[7,]
irf23<-irf[8,]
irf33<-irf[9,]

par(mfrow= c(3,3))

plot(0:hor,irf13,type="l",col="black",main="impulse 3, response 1", xlab="time")
abline(h=0, col="red", lty="dashed")
abline(v=0,col="black" , lty="dashed")

plot(0:hor,irf23,type="l",col="black",main="impulse 3, response 2", xlab="time")
abline(h=0, col="red", lty="dashed")
abline(v=0,col="black" , lty="dashed")

plot(0:hor,irf33,type="l",col="black",main="impulse 3, response 3", xlab="time")
abline(h=0, col="red", lty="dashed")
abline(v=0,col="black" , lty="dashed")

plot(0:hor,irf12,type="l",col="black",main="impulse 2, response 1", xlab="time")
abline(h=0, col="red", lty="dashed")
abline(v=0,col="black" , lty="dashed")

plot(0:hor,irf22,type="l",col="black",main="impulse 2, response 2", xlab="time")
abline(h=0, col="red", lty="dashed")
abline(v=0,col="black" , lty="dashed")

plot(0:hor,irf32,type="l",col="black",main="impulse 2, response 3", xlab="time")
abline(h=0, col="red", lty="dashed")
abline(v=0,col="black" , lty="dashed")

plot(0:hor,irf11,type="l",col="black",main="impulse 1, response 1", xlab="time")
abline(h=0, col="red", lty="dashed")
abline(v=0,col="black" , lty="dashed")

plot(0:hor,irf21,type="l",col="black",main="impulse 1, response 2", xlab="time")
abline(h=0, col="red", lty="dashed")
abline(v=0,col="black" , lty="dashed")

plot(0:hor,irf31,type="l",col="black",main="impulse 1, response 3", xlab="time")
abline(h=0, col="red", lty="dashed")
abline(v=0,col="black" , lty="dashed")













