wavesUI <- 
  function(id) {
  ns <- NS(id)
tagList(
  selectInput(ns("Input"), "P waves", choices = c("Control", "BTZ", "CFZ", "ATRA", "CFZATRA")),
  plotOutput(ns("plot3"))
)
  }



customPlot3 <- function(input, output, session) {
  

output$plot3 <- renderPlot({
  
  
  if(input$Input == "Control") {
    
    x2<-x
    y2<-y
    datsm<-as.data.frame(supsmu(x2,y2))
    x2=datsm$x
    y2=datsm$y
    
    smooth <- function(x, y, w=60, ...) {  
      require(zoo)
      n <- length(y)
      y.smooth<-data.detrend
      y.med <- rollapply(zoo(y.smooth), 1*w+100, FUN=max, align="center") #2*w= 2 only highest peak, 1 rising peaks, 0.5 rising and decreasing peaks
      delta <- y.med - y.smooth[-c(1:w, n+1000:w)] 
      i.max <-which(delta == 0) + w
      list(y.hat=y.smooth, x=x[i.max], i=i.max)
    }
    
    pwave <- function(w,span){
      pwave <-smooth(x2,y2,w=w,span=span)
      position_pwaves <- c() #create and empty vector to store values
      return(as.list(x2[pwave$i]))
    } 
    
    
    position_pwaves<-pwave(750, 200)
    
  wave <- function(w,span){

    wave <-smooth(x,y,w=w,span=span)
    z2 <- NULL
    get_prepeaks <- function(y2){
      z2 <- NULL
      for(i in 2:length(y2)){
        sub_ddt <- ddt[seq(as.numeric(y2[i-1]) , as.numeric(y2[i]), 1) , ]
        sub_zero <- which(abs(sub_ddt$data.detrend) <0.06)
        if(length(sub_zero) >0) z2[i-1] <- sub_ddt$time[tail(sub_zero,1)]}
      return(z2) }
    
    z3 <- NULL
    get_prepeaks2 <- function(y2){
      z3 <- NULL
      for(i in 2:length(y2)){
        sub_ddt <- ddt[seq(as.numeric(y2[i-1]) , as.numeric(y2[i]), 1) , ]
        sub_zero2 <- which(abs(sub_ddt$data.detrend) <0.06)
        if(length(sub_zero2) >0) z3[i-1] <- sub_ddt$time[head(sub_zero2,1)]}
      return(z3)}
    
    plot(x, wave$y.hat,  lwd=2, type="l", col="black", xlab = "Time (Sec)",xlim=c(295000,305000),
         main="P waves and Intervals",adj=0, bty="n",ylab="Voltage", ylim=c(-4,5), xaxt='n')
    points(x[wave$i], wave$y.hat[wave$i], col="Blue", pch=18, cex=2)
    temp_PR <- cbind(get_prepeaks(position_pwaves), 0)
    temp_QRS <- cbind(get_prepeaks2(position_pwaves), 0)
    points(temp_PR, col="Red", pch=15, cex=1)
    points(temp_QRS, col="Red", pch=15, cex=1)
  } 
  
  wave(685, 200) #if FS= 10000, then W=900. If FS=500, then W= 190
  
  }else  
if(input$Input == "BTZ") {   
  x2<-ddt_BTZ$time
  y2<-ddt_BTZ$data.detrend
  
  datsm<-as.data.frame(supsmu(x2,y2))
  x2=datsm$x
  y2=datsm$y
  
  smooth2 <- function(x2, y2, w=60, ...) {   #smoothing data
    n <- length(y2)
    y.smooth<-data.detrend_BTZ
    y.med <- rollapply(zoo(y.smooth), 0.75*w+100, FUN=max, align="center") #2*w= 2 only highest peak, 1 rising peaks, 0.5 rising and decreasing peaks
    delta <- y.med - y.smooth[-c(1:w, n+1000:w)] 
    i.max <-which(delta == 0) + w
    list(y.hat=y.smooth, x=x[i.max], i=i.max)
  }
  
  pwave <- function(w,span){
    pwave <-smooth2(x2,y2,w=w,span=span)
    position_pwaves <- c() #create and empty vector to store values
    return(as.list(x2[pwave$i]))
  } 
  
  
  position_pwaves<-pwave(800, 200)
  
  wave2 <- function(w,span){ #PR
      wave2 <-smooth2(x2,y2,w=w,span=span)
      z2 <- NULL
      get_prepeaks <- function(y2){
        z2 <- NULL
        for(i in 2:length(y2)){
          sub_ddt <- ddt_BTZ[seq(as.numeric(y2[i-1]) , as.numeric(y2[i]), 1) , ]
          sub_zero <- which(abs(sub_ddt$data.detrend) <0.01)
          if(length(sub_zero) >0) z2[i-1] <- sub_ddt$time[tail(sub_zero,1)]}
        return(z2) }
      
      z3 <- NULL
      get_prepeaks2 <- function(y2){ #QRS
        z3 <- NULL
        for(i in 2:length(y2)){
          sub_ddt <- ddt_BTZ[seq(as.numeric(y2[i-1]) , as.numeric(y2[i]), 1) , ]
          sub_zero2 <- which(abs(sub_ddt$data.detrend) <0.005)
          if(length(sub_zero2) >0) z3[i-1] <- sub_ddt$time[head(sub_zero2,1)]}
        return(z3)}
      
    plot(x2, wave2$y.hat,  lwd=2, type="l", xlab= "Number of Observations", xlim=c(20000,30000), 
         main="P waves and Intervals",adj=0, bty="n",ylab="Voltage", ylim=c(-4,5), xaxt='n')
    points(x2[wave2$i], wave2$y.hat[wave2$i], col="Blue", pch=15, cex=1)
    temp_PR <- cbind(get_prepeaks(position_pwaves), 0)
    temp_QRS <- cbind(get_prepeaks2(position_pwaves), 0)
    points(temp_PR, col="Red", pch=15, cex=1)
    points(temp_QRS, col="red", pch=15, cex=1)
    
  } 
  
  wave2(600, 200) #if FS= 10000, then W=900. If FS=500, then W= 190
  
}else

if(input$Input == "CFZ") {   
  
  x2<-x3
  y2<-y3
  
  datsm<-as.data.frame(supsmu(x2,y2))
  x2=datsm$x
  y2=datsm$y
  
  
  smooth3 <- function(x2, y2, w=60, ...) {   #smoothing data
    n <- length(y2)
    y.smooth<-data.detrend_CFZ
    y.med <- rollapply(zoo(y.smooth), 1*w+100, FUN=max, align="center") #2*w= 2 only highest peak, 1 rising peaks, 0.5 rising and decreasing peaks
    delta <- y.med - y.smooth[-c(1:w, n+1000:w)] 
    i.max <-which(delta <= 0) + w
    list(y.hat=y.smooth, x=x[i.max], i=i.max)
  }
  
    pwave <- function(w,span){
      pwave <-smooth3(x2,y2,w=w,span=span)
      position_pwaves <- c() #create and empty vector to store values
      return(as.list(x2[pwave$i]))
    } 
    
    
    position_pwaves<-pwave(800, 200)
    
    wave3 <- function(w,span){
      wave3 <-smooth3(x2,y2,w=w,span=span)
      z2 <- NULL
      get_prepeaks <- function(y2){
        z2 <- NULL
        for(i in 2:length(y2)){
          sub_ddt <- ddt_CFZ[seq(as.numeric(y2[i-1]) , as.numeric(y2[i]), 1) , ]
          sub_zero <- which(abs(sub_ddt$data.detrend) <0.009)
          if(length(sub_zero) >0) z2[i-1] <- sub_ddt$time[tail(sub_zero,1)]}
        return(z2) }
      
      z3 <- NULL
      get_prepeaks2 <- function(y2){
        z3 <- NULL
        for(i in 2:length(y2)){
          sub_ddt <- ddt_CFZ[seq(as.numeric(y2[i-1]) , as.numeric(y2[i]), 1) , ]
          sub_zero2 <- which(abs(sub_ddt$data.detrend) <0.009)
          if(length(sub_zero2) >0) z3[i-1] <- sub_ddt$time[head(sub_zero2,1)]}
        return(z3)}
      
      
    plot(x3, wave3$y.hat,  lwd=2, type="l", col="black", xlab = "Time (Sec)",xlim=c(295000,305000),
         main="P waves and Intervals", bty="n",ylab="Voltage", ylim=c(-4,5), xaxt='n')
    points(x3[wave3$i], wave3$y.hat[wave3$i], col="Blue", pch=18, cex=2)
    temp_PR <- cbind(get_prepeaks(position_pwaves), 0)
    temp_QRS <- cbind(get_prepeaks2(position_pwaves), 0)
    points(temp_PR, col="Red", pch=15, cex=1)
    points(temp_QRS, col="Red", pch=15, cex=1)
    
    } 
    
  wave3(790, 200) #if FS= 10000, then W=900. If FS=500, then W= 190
  
} else

if(input$Input == "ATRA") {  
  
  x2<-x4
  y2<-y4
  
  datsm<-as.data.frame(supsmu(x2,y2))
  x2=datsm$x
  y2=datsm$y
  
  smooth4 <- function(x4, y4, w=60, ...) {   #smoothing data
    require(zoo)
    n <- length(y4)
    y.smooth<-data.detrend_ATRA
    y.med <- rollapply(zoo(y.smooth), 1*w+100, FUN=max, align="center") #2*w= 2 only highest peak, 1 rising peaks, 0.5 rising and decreasing peaks
    delta <- y.med - y.smooth[-c(1:w, n+1000:w)] 
    i.max <-which(delta == 0) + w
    list(y.hat=y.smooth, x=x[i.max], i=i.max)
  }
  
  pwave <- function(w,span){
    pwave <-smooth4(x2,y2,w=w,span=span)
    position_pwaves <- c() #create and empty vector to store values
    return(as.list(x2[pwave$i]))
  } 
  
  
  position_pwaves<-pwave(800, 200)
  
  wave4 <- function(w,span){
    wave4 <-smooth4(x2,y2,w=w,span=span)
    z2 <- NULL
    get_prepeaks <- function(y2){
      z2 <- NULL
      for(i in 2:length(y2)){
        sub_ddt <- ddt_ATRA[seq(as.numeric(y2[i-1]) , as.numeric(y2[i]), 1) , ]
        sub_zero <- which(abs(sub_ddt$data.detrend) <0.1)
        if(length(sub_zero) >0) z2[i-1] <- sub_ddt$time[tail(sub_zero,1)]}
      return(z2) }
    
    z3 <- NULL
    get_prepeaks2 <- function(y2){
      z3 <- NULL
      for(i in 2:length(y2)){
        sub_ddt <- ddt_ATRA[seq(as.numeric(y2[i-1]) , as.numeric(y2[i]), 1) , ]
        sub_zero2 <- which(abs(sub_ddt$data.detrend) <0.006)
        if(length(sub_zero2) >0) z3[i-1] <- sub_ddt$time[head(sub_zero2,1)]}
      return(z3)}
    
  
    plot(x4, wave4$y.hat,  lwd=2, type="l", col="black", xlab = "Time (Sec)",xlim=c(295000,305000),
         main="P waves and Intervals",adj=0, bty="n",ylab="Voltage", ylim=c(-4,5), xaxt='n')
    points(x4[wave4$i], wave4$y.hat[wave4$i], col="Blue", pch=18, cex=2)
    temp_PR <- cbind(get_prepeaks(position_pwaves), 0)
    temp_QRS <- cbind(get_prepeaks2(position_pwaves), 0)
    points(temp_PR, col="Red", pch=15, cex=1)
    points(temp_QRS, col="red", pch=15, cex=1)
    
  } 
  wave4(600, 200) #if FS= 10000, then W=900. If FS=500, then W= 190
  
}else

if(input$Input == "CFZATRA") { 
  
  x2<-x5
  y2<-y5
  
  datsm<-as.data.frame(supsmu(x2,y2))
  x2=datsm$x
  y2=datsm$y
  
  smooth5 <- function(x5, y5, w=60, ...) {   #smoothing data
    require(zoo)
    n <- length(y5)
    y.smooth<-data.detrend_CFZATRA
    y.med <- rollapply(zoo(y.smooth), 1*w+100, FUN=max, align="center") #2*w= 2 only highest peak, 1 rising peaks, 0.5 rising and decreasing peaks
    delta <- y.med - y.smooth[-c(1:w, n+1000:w)] 
    i.max <-which(delta <= 0) + w
    list(y.hat=y.smooth, x=x[i.max], i=i.max)
  }
  
  pwave <- function(w,span){
    pwave <-smooth5(x2,y2,w=w,span=span)
    position_pwaves <- c() #create and empty vector to store values
    return(as.list(x2[pwave$i]))
  } 
  
  
  position_pwaves<-pwave(800, 200)
  
  wave5 <- function(w,span){
    wave5 <-smooth5(x2,y2,w=w,span=span)
    z2 <- NULL
    get_prepeaks <- function(y2){
      z2 <- NULL
      for(i in 2:length(y2)){
        sub_ddt <- ddt_CFZATRA[seq(as.numeric(y2[i-1]) , as.numeric(y2[i]), 1) , ]
        sub_zero <- which(abs(sub_ddt$data.detrend) <0.006)
        if(length(sub_zero) >0) z2[i-1] <- sub_ddt$time[tail(sub_zero,1)]}
      return(z2) }
    
    z3 <- NULL
    get_prepeaks2 <- function(y2){
      z3 <- NULL
      for(i in 2:length(y2)){
        sub_ddt <- ddt_CFZATRA[seq(as.numeric(y2[i-1]) , as.numeric(y2[i]), 1) , ]
        sub_zero2 <- which(abs(sub_ddt$data.detrend) <0.009)
        if(length(sub_zero2) >0) z3[i-1] <- sub_ddt$time[head(sub_zero2,1)]}
      return(z3)}
    
    plot(x5, wave5$y.hat,  lwd=2, type="l", col="black", xlab = "Time (Sec)",xlim=c(290000,300000),
         main="P waves and Intervals",adj=0, bty="n",ylab="Voltage", ylim=c(-4,5), xaxt='n')
    points(x5[wave5$i], wave5$y.hat[wave5$i], col="Blue", pch=18, cex=2)
    temp_PR <- cbind(get_prepeaks(position_pwaves), 0)
    temp_QRS <- cbind(get_prepeaks2(position_pwaves), 0)
    points(temp_PR, col="Red", pch=15, cex=1)
    points(temp_QRS, col="Red", pch=15, cex=1)
    
  } 
  wave5(590, 200) #if FS= 10000, then W=900. If FS=500, then W= 190
}

}) 

return()
}

