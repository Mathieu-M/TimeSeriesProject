# Packages ----------------------------------------------------------------

library("ggplot2")


# Data --------------------------------------------------------------------

ipi <- window(ts(read.table("Data/IPI.dat"), start = 1990, freq = 12), start = 1995)

plot(ipi,main="IPI",type="o")
abline(v=1990:2015,col=4,lty=3)


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

plot(d12logipi)
abline(h=0)

d1d12logipi <- diff(d12logipi,1)

plot(d1d12logipi)
abline(h=0)

var(ipi)
var(logipi)
var(d12logipi)
var(d1d12logipi)

# We select d1d12logipi

ipi <- d1d12logipi


## Question b

par(mfrow=c(1,2))
acf(ipi, ylim=c(-1,1),col=c(2,rep(1,11)),lag.max=84,main="ACF IPI")
pacf(ipi, ylim=c(-1,1),col=c(2,rep(1,11)),lag.max=84,main="PACF IPI")
par(mfrow=c(1,1))
# AR(3) for the seasonal part. 

par(mfrow=c(1,2))
acf(ipi, ylim=c(-1,1),col=c(2,rep(1,11)),lag.max=36,main="ACF IPI")
pacf(ipi, ylim=c(-1,1),col=c(2,rep(1,11)),lag.max=36,main="PACF IPI")
par(mfrow=c(1,1))
# AR(6) or AR(2) for the regular part.

# The two possible models are: ARIMA(6,0,0)(3,1,0)12 or ARIMA(2,0,0)(3,1,0)12


# Estimation --------------------------------------------------------------



