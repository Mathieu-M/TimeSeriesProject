# Plot of the time serie --------------------------------------------------

tsggplot <- function(ts,title=NULL){ # plot the time serie.
  # ts must be a monthly times serie. A title can be add, must be character.
  ts.data.frame <- data.frame(date=as.Date(as.yearmon(time(ts))),as.matrix(ts))
  colnames(ts.data.frame) <- c("time","value")
  ggplot(data=ts.data.frame, mapping=aes(x=time, y=value))+geom_line() + 
    ggtitle(title) + theme(panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank()) + ylab("Values of the serie")
}


# Plot of the ACF and PACF ------------------------------------------------

acfggplot<-function(ts){ # acf and pacf for a ts object.
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
  pacfPlot <- ggplot(data.frame(lag=c(0:71),pacf=ts.pacf$acf)) +
    geom_hline(aes(yintercept=z,colour=type,linetype=type),hline.data) +
    geom_linerange(aes(x=lag,ymin=0,ymax=pacf),
                   colour=c(rep(c("red",rep("black",11)),6))) +
    scale_colour_manual(values = c("black","blue")) +
    scale_linetype_manual(values =c("solid","dashed")) +
    ggtitle("Partial ACF") + scale_y_continuous(limits=c(-1, 1))
  grid.arrange(acfPlot,pacfPlot,ncol=2)
}