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

# We will perform the residuals analysis on the models 3 and 4.


# Validation --------------------------------------------------------------

# Question a

# model 3

# Constant variance (homoscedasticity)
resid <- residplot(mod3)
scatter <- scatterggplot(mod3)
grid.arrange(resid,scatter,ncol=2)
# Pretty good

# Normal residuals
qqggplot(mod3)
# Almost ok

hist(mod3$residuals,breaks=20,freq=F)
curve(dnorm(x,mean=0,sd=sd(mod3$residuals)),col=2,add=T)
# Almost ok

# Independance of the residuals
acfts(mod3$residuals,"Residuals")
# Independant (significant lag are far away)

# Volatility
acfts(mod3$residuals,"Residuals²")
# No volatility

# White noise test (above 0.05 => wn)
ljungggplot(mod3)
# Someissues at the end


# model 4

# Constant variance (homoscedasticity)
resid <- residplot(mod4)
scatter <- scatterggplot(mod4)
grid.arrange(resid,scatter,ncol=2)
# Some outliers but seem constant

# Normal residuals
qqggplot(mod4)
# Almost ok

hist(mod4$residuals,breaks=20,freq=F)
curve(dnorm(x,mean=0,sd=sd(mod4$residuals)),col=2,add=T)
# Seem normal

# Independance of the residuals
acfts(mod4$residuals,"Residuals")
# Independant

# Volatility
acfts(mod4$residuals,"Residuals²")
# No volatility

# White noise test (above 0.05 => wn)
ljungggplot(mod4)
# OK


# Question b

cat("\nModul of AR Characteristic polynomial Roots: ", 
    Mod(polyroot(c(1,-mod3$model$phi))),"\n")
cat("\nModul of MA Characteristic polynomial Roots: ",
    Mod(polyroot(c(1,mod3$model$theta))),"\n")
# Stationary and invertible

cat("\nModul of AR Characteristic polynomial Roots: ", 
    Mod(polyroot(c(1,-mod4$model$phi))),"\n")
cat("\nModul of MA Characteristic polynomial Roots: ",
    Mod(polyroot(c(1,mod4$model$theta))),"\n")
# Stationary and invertible


# Question c

ultim=c(2013,12)
ipi <- window(ts(read.table("Data/IPI.dat"), start = 1990, freq = 12), start = 1995)

ipi2 <- window(ipi,end=ultim)
logipi2 <- log(ipi2)
d12logipi2 <- diff(logipi2,12)
d1d12logipi2 <- diff(d12logipi2,1)
ipi2 <- d1d12logipi2
# We look at the same serie without the last 12 observations.

mod3bis <- arima(ipi2,order=c(6,0,0),seasonal=list(order=c(3,0,2),period=12))
mod4bis <- arima(ipi2,order=c(6,0,0),seasonal=list(order=c(3,0,5),period=12))

mod3$coef
mod3bis$coef
# The modelis stable.

mod4$coef
mod4bis$coef
# The model is stable

# TODO: take out the non-significant parameters to see if it improves the model.
