# Packages ----------------------------------------------------------------

library("ggplot2")
library("gridExtra")
library("zoo")


# Data --------------------------------------------------------------------

ipi <- window(ts(read.table("Data/IPI.dat"), start = 1990, freq = 12), start = 1995)

tsggplot(ipi,title="IPI")


# Identification ----------------------------------------------------------

## Question a

m <- apply(matrix(ipi,nr=12),2,mean)
v <- apply(matrix(ipi,nr=12),2,var)
qplot(m,v,xlab="Yearly mean",ylab="Yearly variance ",main="IPI") + 
  stat_smooth(method="lm", se=FALSE)

boxplot(ipi~floor(time(ipi)))
# Variance does not seem to be constant.

logipi <- log(ipi)

m <- apply(matrix(logipi,nr=12),2,mean)
v <- apply(matrix(logipi,nr=12),2,var)
qplot(m,v,xlab="Yearly mean",ylab="Yearly variance ",main="logIPI") + 
  stat_smooth(method="lm", se=FALSE)

boxplot(logipi~floor(time(logipi)))
# It looks better for the boxplots.

plot(decompose(logipi))
monthplot(logipi)
# There is a clear seasonal pattern

d12logipi <- diff(logipi,12)

tsggplot(d12logipi,"d12logipi") + geom_hline(y=0)

d1d12logipi <- diff(d12logipi,1)

tsggplot(d1d12logipi) + geom_hline(y=0)

var(ipi)
var(logipi)
var(d12logipi)
var(d1d12logipi)

# We select d1d12logipi

ipi <- d1d12logipi


## Question b

acfggplot(ipi)
# ARMA(3,2) or ARMA(3,5) for the seasonal part. 
# AR(6) or AR(2) for the regular part.

# The two possible models are: ARIMA(6,0,0)(3,1,2)12 or ARIMA(6,0,0)(3,1,5)12
# or ARIMA(2,0,0)(3,1,2)12 or ARIMA(2,0,0)(3,1,5)12


# Estimation --------------------------------------------------------------

mod1 <- arima(ipi,order=c(2,0,0),seasonal=list(order=c(3,0,2),period=12))
mod2 <- arima(ipi,order=c(2,0,0),seasonal=list(order=c(3,0,5),period=12))
mod3 <- arima(ipi,order=c(6,0,0),seasonal=list(order=c(3,0,2),period=12)) 
mod4 <- arima(ipi,order=c(6,0,0),seasonal=list(order=c(3,0,5),period=12)) 

# Model 1
acfggplot(ipi)
par(mfrow=c(1,2))
plot(ARMAacf(mod1$model$phi,mod1$model$theta,lag.max=36),ylim=c(-1,1), 
     type="h",xlab="Lag",  ylab="", main="ACF Teoric",col=c(2,rep(1,11)))
abline(h=0)
plot(ARMAacf(mod1$model$phi,mod1$model$theta,lag.max=36, pacf=T),ylim=c(-1,1),
     type="h", xlab="Lag", ylab="", main="PACF Teoric",col=c(2,rep(1,11)))
abline(h=0)
par(mfrow=c(1,1))
# Really good for ACF, quite good for PACF.

# Model 2
acfggplot(ipi)
par(mfrow=c(1,2))
plot(ARMAacf(mod2$model$phi,mod2$model$theta,lag.max=36),ylim=c(-1,1), 
     type="h",xlab="Lag",  ylab="", main="ACF Teoric",col=c(2,rep(1,11)))
abline(h=0)
plot(ARMAacf(mod2$model$phi,mod2$model$theta,lag.max=36, pacf=T),ylim=c(-1,1),
     type="h", xlab="Lag", ylab="", main="PACF Teoric",col=c(2,rep(1,11)))
abline(h=0)
par(mfrow=c(1,1))
# Quite good for ACF and PACF.

# Model 3
acfggplot(ipi)
par(mfrow=c(1,2))
plot(ARMAacf(mod3$model$phi,mod3$model$theta,lag.max=36),ylim=c(-1,1), 
     type="h",xlab="Lag",  ylab="", main="ACF Teoric",col=c(2,rep(1,11)))
abline(h=0)
plot(ARMAacf(mod3$model$phi,mod3$model$theta,lag.max=36, pacf=T),ylim=c(-1,1),
     type="h", xlab="Lag", ylab="", main="PACF Teoric",col=c(2,rep(1,11)))
abline(h=0)
par(mfrow=c(1,1))
# Really good for ACF and PACF.

# Model 4
acfggplot(ipi)
par(mfrow=c(1,2))
plot(ARMAacf(mod4$model$phi,mod4$model$theta,lag.max=36),ylim=c(-1,1), 
     type="h",xlab="Lag",  ylab="", main="ACF Teoric",col=c(2,rep(1,11)))
abline(h=0)
plot(ARMAacf(mod4$model$phi,mod4$model$theta,lag.max=36, pacf=T),ylim=c(-1,1),
     type="h", xlab="Lag", ylab="", main="PACF Teoric",col=c(2,rep(1,11)))
abline(h=0)
par(mfrow=c(1,1))
# Really good for ACF and PACF, better than the previous one.

# WE will perform the residuals analysis on the models 3 and 4.


# Validation --------------------------------------------------------------


