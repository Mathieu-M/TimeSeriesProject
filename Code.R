# Packages ----------------------------------------------------------------

library("ggplot2")
library("gridExtra")
library("zoo")
source('PlotTimeSeriesFunctions.R')


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

ipi.t <- d1d12logipi


## Question b

acfts(ipi.t)
# ARMA(3,2) or ARMA(3,5) for the seasonal part. 
# AR(6) or AR(2) for the regular part.

# The two possible models are: ARIMA(6,0,0)(3,1,2)12 or ARIMA(6,0,0)(3,1,5)12
# or ARIMA(2,0,0)(3,1,2)12 or ARIMA(2,0,0)(3,1,5)12


# Estimation --------------------------------------------------------------

mod1 <- arima(ipi.t,order=c(2,0,0),seasonal=list(order=c(3,0,2),period=12))
mod2 <- arima(ipi.t,order=c(2,0,0),seasonal=list(order=c(3,0,5),period=12))
mod3 <- arima(ipi.t,order=c(6,0,0),seasonal=list(order=c(3,0,2),period=12)) 
mod4 <- arima(ipi.t,order=c(6,0,0),seasonal=list(order=c(3,0,5),period=12)) 

mod1bis <- arima(logipi,order=c(2,1,0),seasonal=list(order=c(3,1,2),period=12))
mod2bis <- arima(logipi,order=c(2,1,0),seasonal=list(order=c(3,1,5),period=12))
mod3bis <- arima(logipi,order=c(6,1,0),seasonal=list(order=c(3,1,2),period=12)) 
mod4bis <- arima(logipi,order=c(6,1,0),seasonal=list(order=c(3,1,5),period=12)) 
# We first check if the constant is needed. Since it is not we work with the models 
# without the intercept. (estimate/se<2)

## Model 1

acfmodel(mod1)
acfts(ipi.t)
acfmodel(mod1bis)
# Really good for ACF, quite good for PACF.

## Model 2

acfmodel(mod2)
acfts(ipi.t)
acfmodel(mod2bis)
# Quite good for ACF and PACF.

## Model 3

acfmodel(mod3)
acfts(ipi.t)
acfmodel(mod3bis)
# Really good for ACF and PACF.

## Model 4

acfmodel(mod4)
acfts(ipi.t)
acfmodel(mod4bis)
# Really good for ACF and PACF, better than the previous one.

# We will perform the residuals analysis on the models 3 and 4.


# Validation --------------------------------------------------------------

## Question a

# model 3

# Constant variance (homoscedasticity)
resid <- residplot(mod3bis)
scatter <- scatterggplot(mod3bis)
grid.arrange(resid,scatter,ncol=2)
# Pretty good

# Normal residuals
qqggplot(mod3bis)
# Almost ok

hist(mod3bis$residuals,breaks=20,freq=F)
curve(dnorm(x,mean=0,sd=sd(mod3bis$residuals)),col=2,add=T)
# Almost ok

# Independance of the residuals
acfts(mod3bis$residuals,"Residuals")
# Independant (significant lag are far away)

# Volatility
acfts(mod3bis$residuals,"Residuals²")
# No volatility

# White noise test (above 0.05 => wn)
ljungggplot(mod3bis)
# Someissues at the end

# model 4

# Constant variance (homoscedasticity)
resid <- residplot(mod4bis)
scatter <- scatterggplot(mod4bis)
grid.arrange(resid,scatter,ncol=2)
# Some outliers but seem constant

# Normal residuals
qqggplot(mod4bis)
# Almost ok

hist(mod4bis$residuals,breaks=20,freq=F)
curve(dnorm(x,mean=0,sd=sd(mod4bis$residuals)),col=2,add=T)
# Seem normal

# Independance of the residuals
acfts(mod4bis$residuals,"Residuals")
# Independant

# Volatility
acfts(mod4bis$residuals,"Residuals²")
# No volatility

# White noise test (above 0.05 => wn)
ljungggplot(mod4bis)
# OK


## Question b

cat("\nModul of AR Characteristic polynomial Roots: ", 
    Mod(polyroot(c(1,-mod3$model$phi))),"\n")
cat("\nModul of MA Characteristic polynomial Roots: ",
    Mod(polyroot(c(1,mod3$model$theta))),"\n")

length(Mod(polyroot(c(1,-mod3$model$phi)))[Mod(polyroot(c(1,-mod3$model$phi)))<1])
length(Mod(polyroot(c(1,mod3$model$theta)))[Mod(polyroot(c(1,mod3$model$theta)))<1])
# Stationary and invertible

cat("\nModul of AR Characteristic polynomial Roots: ", 
    Mod(polyroot(c(1,-mod4$model$phi))),"\n")
cat("\nModul of MA Characteristic polynomial Roots: ",
    Mod(polyroot(c(1,mod4$model$theta))),"\n")

length(Mod(polyroot(c(1,-mod4$model$phi)))[Mod(polyroot(c(1,-mod4$model$phi)))<1])
length(Mod(polyroot(c(1,mod4$model$theta)))[Mod(polyroot(c(1,mod4$model$theta)))<1])
# Stationary and invertible


## Question c

# Stability

ultim=c(2013,12)

ipi2 <- window(ipi,end=ultim)
logipi2 <- log(ipi2)
# We look at the same serie without the last 12 observations.

mod3bis2 <- arima(logipi2,order=c(6,1,0),seasonal=list(order=c(3,1,2),period=12))
mod4bis2 <- arima(logipi2,order=c(6,1,0),seasonal=list(order=c(3,1,5),period=12))

mod3bis$coef
mod3bis2$coef
# The model is stable.

mod4bis$coef
mod4bis2$coef
# The model is stable

# We can to take out the non-significant parameters to see if it improves the model. (estimates/se<2)

mod3bis.adjusted <- arima(logipi,order=c(6,1,0),seasonal=list(order=c(1,1,2),period=12),
                     fixed=c(NA,NA,0,0,NA,NA,NA,0,NA)) 
# After several test this model is better (aic smaller). Is it still a good fit?

acfts(ipi.t)
acfmodel(mod3bis.adjusted)
# We keep this model as the new mod3bis.

mod3bis <- mod3bis.adjusted

ar3,4,sma3
mod4bis.adjusted <- arima(logipi,order=c(6,1,0),seasonal=list(order=c(3,1,5),period=12),
                     fixed=c(NA,NA,NA,0,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)) 
# The aic is better with all the parameters.

# Prediction for model 3.

pdq <- c(6,1,0)
PDQ <- c(3,1,2)

mod3bis <- arima(logipi2,order=pdq,seasonal=list(order=PDQ,period=12))

pred <- predict(mod3bis,n.ahead=12)

pr <- window(diffinv(pred$pred,12,xi=window(ipi,start=ultim+c(-1,1),end=ultim+c(0,0))),start=ultim)
model <- mod3bis$model
varZ <- mod3bis$sigma
ma <- ARMAtoMA(ar=mod3bis$phi,ma=mod3bis$theta,lag.max=11)
ar <- ARMAtoMA(ar=mod3bis$theta,ma=mod3bis$phi,lag.max=11)
se <- c(0,sqrt((cumsum(c(1,ma,ar))^2)*varZ))

tl<-ts(pr-1.96*se,start=ultim,freq=12)
tu<-ts(pr+1.96*se,start=ultim,freq=12)
pr<-ts(pr,start=ultim,freq=12)

ts.plot(ipi,tl,tu,pr,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2011,2015),type="o",main=paste("Model ARIMA(",paste(pdq,collapse=","),")(",paste(PDQ,collapse=","),")12",sep=""))
abline(v=2010+0:5,lty=3,col=4)

pred <- predict(mod3bis,n.ahead=12)
pr <- ts(c(tail(d1d12logipi2,1),pred$pred),start=ultim,freq=12)
se <- ts(c(0,pred$se),start=ultim,freq=12)

#Intervals
tl <- ts(exp(pr-1.96*se),start=ultim,freq=12)
tu <- ts(exp(pr+1.96*se),start=ultim,freq=12)
pr <- ts(exp(pr),start=ultim,freq=12)


ts.plot(d1d12logipi,tl,tu,pr,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2010,2015),type="o",main=paste("Model ARIMA(",paste(pdq,collapse=","),")(",paste(PDQ,collapse=","),")12",sep=""))
abline(v=2007+0:8,lty=3,col=4)
