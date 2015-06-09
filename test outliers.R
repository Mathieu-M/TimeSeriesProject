source("outlierTreatment.r")


mod.atip4=outdetec(mod4bis,dif=c(1,12),crit=2.6,LS=T)

atipics4=mod.atip4$atip[order(mod.atip4$atip[,1]),]
meses=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
data.frame(atipics4,Fecha=paste(meses[(atipics4[,1]-1)%%12+1],start(logipi)[1]+((atipics4[,1]-1)%/%12)))
mod.atip4$sigma2

data.frame(atipics4,Fecha=paste(meses[(atipics4[,1]-1)%%12+1],start(logipi)[1]+((atipics4[,1]-1)%/%12)),perc.Obs=exp(atipics4[,3])*100)



tsggplot(ipi.t)



mod.atip3=outdetec(mod3bis,dif=c(1,12),crit=2.6,LS=T)

atipics3=mod.atip3$atip[order(mod.atip3$atip[,1]),]
meses=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
data.frame(atipics3,Fecha=paste(meses[(atipics3[,1]-1)%%12+1],start(logipi)[1]+((atipics3[,1]-1)%/%12)))
mod.atip3$sigma2

data.frame(atipics3,Fecha=paste(meses[(atipics3[,1]-1)%%12+1],start(logipi)[1]+((atipics3[,1]-1)%/%12)),perc.Obs=exp(atipics3[,3])*100)


##Comparacion serie observada con la serie linealizada (sin atipicos)
logipi.lin=lineal(logipi,mod.atip3$atip) # original serie and table of outliers
ipi.lin=exp(logipi.lin)

plot(ipi)
lines(ipi.lin,col=2)

##Efecto de los atipicos en la serie de logaritmos
plot(logipi-logipi.lin)
# We plot the effects.


d1d12logipi.lin=diff(diff(logipi.lin,12))
par(mfrow=c(1,2))
acf(d1d12logipi.lin,ylim=c(-1,1),lag.max=36,col=c(2,rep(1,11)))
pacf(d1d12logipi.lin,ylim=c(-1,1),lag.max=36,col=c(rep(1,11),2))
par(mfrow=c(1,1))

acfts(d1d12logipi.lin)

##Estimaci?n del modelo para la serie linealizada
mod3bis.lin <- arima(logipi,order=c(6,1,0),seasonal=list(order=c(1,1,2),period=12),
                     fixed=c(NA,NA,0,0,NA,NA,NA,0,NA))

validation(mod3bis.lin,d1d12logipi)
# There is no problem, the model is validated.

## truncated lineal serie

ultim <- c(2013,12)
pdq <- c(6,1,0)
PDQ <- c(1,1,2)

ipi2.lin <- window(ipi.lin,end=ultim)
logipi2.lin <- log(ipi2.lin)

# Model lineal for the truncated serie.
mod3bis2.lin <- arima(logipi2.lin,order=pdq,seasonal=list(order=PDQ,period=12),
                  fixed=c(NA,NA,0,0,NA,NA,NA,0,NA))



pred <- predict(mod3bis2.lin,n.ahead=12)
pr <- ts(c(tail(logipi2.lin,1),pred$pred),start=ultim,freq=12)
se <- ts(c(0,pred$se),start=ultim,freq=12)

tl1<-ts(exp(pr-1.96*se),start=ultim,freq=12)
tu1<-ts(exp(pr+1.96*se),start=ultim,freq=12)
pr1<-ts(exp(pr),start=ultim,freq=12)

tspredggplot(ipi.lin,pred=pr1,upperb=tu1,lowerb=tl1,title="Predictions for model lineal.")
tspredggplot(ipi1,pred=pr1,upperb=tu1,lowerb=tl1,title="Predictions for model 3.")


obs <- window(ipi,start=ultim)

mod3bis2.EQM <- sqrt(sum(((obs-pr3)/obs)^2)/12)
mod3bis2.EAM <- sum(abs(obs-pr3)/obs)/12

mod3bis2.lin.EQM <- sqrt(sum(((obs-pr1)/obs)^2)/12)
mod3bis2.lin.EAM <- sum(abs(obs-pr1)/obs)/12

mod3bis2.EQM-mod3bis2.lin.EQM<0
mod3bis2.EAM-mod3bis2.lin.EAM<0

## long term

ipi1.lin <- window(ipi.lin,end=ultim+c(1,0))
logipi1.lin <- log(ipi1.lin)

pred.lin <- predict(mod3bis.lin,n.ahead=12)
pr.lin <- ts(c(tail(logipi1.lin,1),pred.lin$pred),start=ultim+c(1,0),freq=12)
se.lin <- ts(c(0,pred.lin$se),start=ultim+c(1,0),freq=12)

tl1.lin<-ts(exp(pr.lin-1.96*se.lin),start=ultim+c(1,0),freq=12)
tu1.lin<-ts(exp(pr.lin+1.96*se.lin),start=ultim+c(1,0),freq=12)
pr1.lin<-ts(exp(pr.lin),start=ultim+c(1,0),freq=12)

tspredggplot(ipi1.lin,pred=pr1.lin,upperb=tu1.lin,lowerb=tl1.lin,title="Predictions for model lineal.")