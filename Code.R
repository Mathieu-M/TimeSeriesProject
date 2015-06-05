# Packages ----------------------------------------------------------------

library("ggplot2")


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
# AR(3) for the seasonal part. 
# AR(6) or AR(2) for the regular part.

# The two possible models are: ARIMA(6,0,0)(3,1,0)12 or ARIMA(2,0,0)(3,1,0)12


# Estimation --------------------------------------------------------------



