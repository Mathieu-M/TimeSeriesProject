# List of plot functions for object of class time series.

# Plot of the time serie --------------------------------------------------

tsggplot <- function(ts,title=NULL){ # plot the time serie.
  # ts must be a monthly times serie. A title can be add, must be character.
  ts.data.frame <- data.frame(date=as.Date(as.yearmon(time(ts))),as.matrix(ts))
  colnames(ts.data.frame) <- c("time","value")
  ggplot(data=ts.data.frame, mapping=aes(x=time, y=value))+geom_line() + 
    ggtitle(title) + theme(panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank()) + ylab("Values of the serie")
}


# Plot of the ACF and PACF of a time serie --------------------------------

acfts<-function(ts,title=NULL){ # acf and pacf for a ts object.
  ts.data.frame <- data.frame(c(1:length(ts)),ts)
  colnames(ts.data.frame) <- c("time","value")
  ts.acf<-acf(ts, plot=FALSE,lag.max=71)
  ts.pacf<-pacf(ts, plot=FALSE,lag.max=72)
  ci <- 0.95
  clim0 <- qnorm((1 + ci)/2)/sqrt(ts.acf$n.used)
  clim <- c(-clim0,clim0)
  hline.data <- data.frame(z=c(0,clim),type=c("base","ci","ci"))
  acfPlot <- ggplot(data.frame(lag=c(0:71),acf=ts.acf$acf)) +
    geom_hline(aes(yintercept=z,colour=type,linetype=type),hline.data) +
    geom_linerange(aes(x=lag,ymin=0,ymax=acf),
                   colour=c(rep(c("red",rep("black",11)),6))) +
    scale_colour_manual(values = c("black","blue")) +
    scale_linetype_manual(values =c("solid","dashed")) + ylab("") + 
    ggtitle("ACF") + scale_y_continuous(limits=c(-1, 1)) + 
    geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1),colour="red")
  pacfPlot <- ggplot(data.frame(lag=c(1:72),pacf=ts.pacf$acf)) +
    geom_hline(aes(yintercept=z,colour=type,linetype=type),hline.data) +
    geom_linerange(aes(x=lag,ymin=0,ymax=pacf),
                   colour=c(rep(c(rep("black",11),"red"),6))) +
    scale_colour_manual(values = c("black","blue")) +
    scale_linetype_manual(values =c("solid","dashed")) +
    ggtitle("Partial ACF") + scale_y_continuous(limits=c(-1, 1))
  grid.arrange(acfPlot,pacfPlot,ncol=2,main=title)
}


# Plot of the ACF and PACF of a ARMAacf -----------------------------------

acfmodel <- function(model){ # acf and pacf for an ARiMA model.
  modelacf <- as.data.frame(ARMAacf(model$model$phi,model$model$theta,
                                    lag.max=36))
  modelacf <- cbind(seq(0,length(modelacf[[1]])-1,by=1),modelacf)
  colnames(modelacf) <- c("lag","acf")
  acfPlot <- ggplot(modelacf,aes(lag,acf)) + 
    geom_segment(aes(x=lag,y=0,xend=lag,yend=acf),
                 colour=c(rep(c("red",rep("black",11)),3),"red")) +
    geom_hline(y=0) + scale_y_continuous(limits=c(-1, 1)) + 
    ggtitle("ACF theoric") + ylab("")
  modelpacf <- as.data.frame(ARMAacf(model$model$phi,model$model$theta,
                                     lag.max=37,pacf=TRUE))
  modelpacf <- cbind(seq(1,length(modelpacf[[1]]-1),by=1),modelpacf)
  colnames(modelpacf) <- c("lag","pacf")
  pacfPlot <- ggplot(modelpacf,aes(lag,pacf)) + 
    geom_segment(aes(x=lag,y=0,xend=lag,yend=pacf),
                 colour=c(rep(c(rep("black",11),"red"),3),"black")) +
    geom_hline(y=0) + scale_y_continuous(limits=c(-1, 1)) + 
    ggtitle("PACF theoric") + ylab("")
  grid.arrange(acfPlot,pacfPlot,ncol=2)
}


# Plot of the residuals ---------------------------------------------------

# Regular plot
residplot <- function(model){
  ci <- 3*sd(model$residuals)
  clim <- c(-ci,ci)
  hline.data <- data.frame(z=clim,type=c("ci","ci"))
  tsggplot(model$residuals) + geom_hline(aes(yintercept=z,linetype=type,
    colour=type),hline.data) + geom_hline(y=0) + ggtitle("Residuals")
}

# QQ-plot
qqggplot <- function(model){
  ts.data.frame <- data.frame(c(1:length(model$residuals)),model$residuals)
  colnames(ts.data.frame) <- c("time","residuals")
  qtype <- 7
  y <- quantile(model$residuals[!is.na(model$residuals)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  intercept <- y[1L] - slope * x[1L]
  ggplot(ts.data.frame, aes(sample = residuals)) + geom_point(stat = "qq") + 
    geom_abline(intercept=intercept,slope=slope,color="red") + 
    ggtitle("QQ-plot")
}

# Square root of the absolute residuals
scatterggplot <- function(model){
  x <- sqrt(abs(model$residuals))
  y <- NULL
  lpars=list(col=2)
  span = 2/3
  family = c("symmetric","gaussian")
  evaluation = 50
  degree = 1
  xlabel <- if (!missing(x)){
    deparse(substitute(x))
  } 
  ylabel <- if (!missing(y)){
    deparse(substitute(y))
  } 
  xy <- xy.coords(x, y, xlabel, ylabel)
  x <- xy$x
  y <- xy$y
  xlab <- xy$xlab
  ylab <- if (is.null(ylab)){
    xy$ylab
  } 
  pred <- loess.smooth(x,y, span=span, degree=degree, family=family, evalution=evaluation)
  ylim = range(y, pred$y,na.rm = TRUE)
  linh <- as.data.frame(c(list(pred),lpars))
  colnames(linh) <- c("pred","lpars","col")
  x1 <- as.data.frame(x=x,y=y)
  ggplot() + geom_point(data=x1, aes(x=x, y=y)) +  
    geom_line(data=linh,aes(x=pred,y=lpars),colour="red") + ylab(expression(sqrt(abs(residuals)))) + xlab("") + 
    ggtitle("Satter-plot with smooth")
}

# Ljung-Box test plot
ljungggplot <- function(model){
  gof.lag <- 7*frequency(get(model$series))
  rs <- model$residuals
  nlag <- gof.lag
  pval <- numeric(nlag)
  for (i in 1L:nlag){
    pval[i] <- Box.test(rs, i, type = "Ljung-Box")$p.value
  } 
  df <- data.frame(c(1:nlag),pval)
  test <- factor(df$pval<0.05)
  df <- cbind(df,test)
  colnames(df) <- c("lag","pval","test")
  ggplot(data=df,aes(x=lag, y=pval)) + geom_point(aes(colour=factor(test))) + 
    geom_hline(y=0.05,linetype="dashed",colour="blue") + geom_hline(y=0) + 
    ggtitle("p values for Ljung-Box statistic") + ylab("") + xlab("") + 
    coord_cartesian(ylim=c(-0.05, 1.05)) + 
    scale_colour_manual(values=c("TRUE"="red","FALSE"="blue"),name="Value of the p-value",
                        breaks=c("TRUE","FALSE"),labels=c("<0.05",">0.05"))
}


# Plot of the arma roots in a unite circle --------------------------------

plotarmaroots <- function(object){
  if(class(object) != "Arima")
    stop("object must be of class Arima or ar")
  if(class(object) == "Arima"){
    parvecphi <- object$model$phi
  }
  if(length(parvecphi) > 0){
    last.nonzero <- max(which(abs(parvecphi) > 1e-08))
    if (last.nonzero > 0){
      arroots <- structure(list(roots=polyroot(c(1,-parvecphi[1:last.nonzero])),
                                type="AR"), class='armaroots')
    }
  } else{
    arroots <- structure(list(roots=numeric(0),type="AR"),class='armaroots')
  }
  parvectheta <- object$model$theta
  if(length(parvectheta) > 0){
    last.nonzero <- max(which(abs(parvectheta) > 1e-08))
    if (last.nonzero > 0){
      maroots <- structure(list(roots=polyroot(c(1,parvectheta[1:last.nonzero])),
                                type="MA"), class='armaroots')
    }
  } else{
    maroots <- structure(list(roots=numeric(0),type="MA"),class='armaroots')
  }
  rootsma <- as.data.frame(1/maroots$roots)
  testma <- factor(Mod(rootsma[[1]])>1)
  rootsma <- cbind(rootsma,testma)
  colnames(rootsma) <- c("invmaroots","test")
  rootsar <- as.data.frame(1/arroots$roots)
  testar <- factor(Mod(rootsar[[1]])>1)
  rootsar <- cbind(rootsar,testar)
  colnames(rootsar) <- c("invarroots","test")
  xc <- 0
  yc <- 0
  r <- 1
  arplot <- ggplot(data=rootsar,aes(x=Re(invarroots), y=Im(invarroots))) + 
    geom_point(aes(colour=factor(testar))) +  
    ggtitle(paste("The",length(arroots$roots),"inverse",arroots$type," roots")) + 
    xlab("Real") + ylab("Imaginary") + annotate("path",x=xc+r*cos(seq(0,2*pi,
                                                                      length.out=100)),y=yc+r*sin(seq(0,2*pi,length.out=100))) + coord_fixed() + 
    scale_colour_manual(values=c("TRUE"="red","FALSE"="blue"),name="1/Roots",
                        breaks=c("TRUE","FALSE"),labels=c(">1","<1")) + 
    geom_hline(y=0,linetype="dashed") + geom_vline(x=0,linetype="dashed") 
  
  maplot <- ggplot(data=rootsma,aes(x=Re(invmaroots), y=Im(invmaroots))) + 
    geom_point(aes(colour=factor(testma))) +  
    ggtitle(paste("The",length(maroots$roots),"inverse",maroots$type," roots")) + 
    xlab("Real") + ylab("Imaginary") + annotate("path",x=xc+r*cos(seq(0,2*pi,
                                                                      length.out=100)),y=yc+r*sin(seq(0,2*pi,length.out=100))) + coord_fixed() + 
    scale_colour_manual(values=c("TRUE"="red","FALSE"="blue"),name="1/Roots",
                        breaks=c("TRUE","FALSE"),labels=c(">1","<1")) + 
    geom_hline(y=0,linetype="dashed") + geom_vline(x=0,linetype="dashed") 
  
  grid.arrange(arplot,maplot,ncol=2,main="Inverse roots of the charasteristic polynomial")
}
