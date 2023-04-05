## Script for Analysis of Mice ECG##
library(pracma)
library(zoo)
library(devtools)
#library(easyGgplot2)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(lattice)
library(plotrix)
library(grid)
library(ggplotify)
library(cowplot)
library(gridExtra)
library(pROC)
library(ggfortify)
library(changepoint)
library(strucchange) 
library(ggpmisc)
library(dlm)
library(TTR)
library(tsbox)
library(forecast)
library(signal)



###########################################################################
#setwd("~/Desktop/Thesis More Data/in vivo ECG/IN VIVO_FINAL")
setwd("/Volumes/NO NAME/LCMS murine hearts 03.12.21")

N.1<- "BTZ1.csv"
#######################MICE 1#####################################################
dat1<- read.csv(N.1, sep=";", header=FALSE) #change "," to ";" accordingly
dat1<- na.omit(dat1)
y<-as.numeric(dat1$V2) - mean(dat1$V2)
length<- as.numeric(length(dat1$V1))
x<-seq(from = 0, to = length/1000, by=0.001)
x<-x[-1]


dt<-as.matrix(cbind(x,y))
dt2<- apply(dt,1,median)

FS=10000

g0<-ggplot(dat1, aes(x=V1, y=V2))+ geom_line(lwd=0.5)+  scale_x_continuous(limits=c(44,45)) +
  scale_y_continuous(limits=c(-5,5))+ xlab("Time (Sec)") + ylab("Voltage (V)") +
  ggtitle("Raw ECG Input") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))

g0 


#####Low pass filter
#H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2
#lpf <- ((1-10^(-6))^2)/(1-10^(-1))^2
#Wn <- (1*2)/FS
N = 3 #order of 3 less processing
lowfilter<- butter(N, 0.1, type= "low")
ecg_l <-filtfilt(lowfilter, y )
plot(ecg_l, type="l", xlim=c(440000, 450000), ylim=c(-10,10)) 
#ecg_l <- (ecg_l/ max(abs(ecg_l)))
ecg_l <- (ecg_l-mean(ecg_l))/sd(ecg_l)
plot(ecg_l, type="l", xlim=c(440000, 450000))

#####High pass filter
highfilter<- butter(N, 0.01, type= "high")
ecg_h <-filtfilt(highfilter, y )
plot(ecg_h, type="l", xlim=c(440000, 450000), ylim=c(-10,10))
#ecg_h <- ecg_h/max(abs(ecg_h))
ecg_h <- (ecg_l-mean(ecg_h))/sd(ecg_h)
plot(ecg_h, type="l", xlim=c(440000, 450000))

ecg_x <- (ecg_l-mean(ecg_l))/sd(ecg_h)
plot(ecg_x, type="l", xlim=c(440000, 450000))

#####Bandpass filter
Wn <- (c(0.7,100)*2)/FS
ab<- butter(N,Wn)
ecg_h2 <- filtfilt(ab, ecg_x, y)
plot(ecg_h2, type="l", xlim=c(440001, 450000))

smooth1<-as.numeric(ecg_h)

##########################################Check if there is a hidden frequency harmonic
FFT<-fft(smooth1,  inverse=TRUE)
#plot(Re(FFT), type="l")


###########################################Apply detrend and Time Series
data.detrend <- detrend(smooth1, tt = 'linear')    

data.detrend <- ts(as.numeric(data.detrend),
                   start=c(first(data.detrend),1), frequency=FS) #original frequency was set to 7...for whatever reason, 100 looks really nice
str(data.detrend) #check

decomposed<- decompose(data.detrend)  # decompose time series to extract the trend
plot(decomposed, xlim=c(1,2))
data.detrend<-decomposed$x - decomposed$trend
data.detrend<-data.detrend - decomposed$seasonal

data.detrend<- na.omit(data.detrend)
plot(data.detrend, xlim=c(44,45))
#write.zoo(data.detrend, file = "data.detrend_CFZ.csv")

ddt<-as.data.frame(data.detrend) 
ddt<-na.omit(ddt)
ddt<-ddt %>% mutate(time= seq(from = 0, to = length(data.detrend[-1]), by=1))
#write.csv(ddt, file="ddt_CFZ.csv")

x<-ddt$time/FS
y<-as.numeric(ddt$x)


THR_SIG <- max(y)*1/3 # 0.25 of the max amplitude 
THR_NOISE = mean(y)*1/2; # 0.5 of the mean signal is considered to be noise



smoothing <- function(x, y, w=5, ...) {   #smoothing data
  require(zoo)
  n <- length(y)
  y.smooth<-data.detrend
  y.max <- rollapply(zoo(y.smooth), 2*w+5, max, align="center")
  y.min <- rollapply(zoo(y.smooth), 2*w+1, min, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+5-1:w)]
  i.max <- which(delta <= 0) + w 
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}


peaks <- function(w, span) {  
  peaks <- smoothing(x, y, w=w, span=span)#span=level of "smoothing"
  y.max <- max(y)
  y.min <- min(y)
  if (peaks$i < THR_SIG) {
  plot(x,y, type="l", col="white", xlab = "Time (Sec)", ylab="V", xlim=c(40,50)) # ylim=c(-4,6),
  lines(x, peaks$y.hat,  lwd=2)
  sapply(peaks$i, function(i) lines(c(x[i],x[i]), c(y.min, peaks$y.hat[i]),col="Red", lty=2))
  points(x[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=15, cex=1)
  position_peaks <- c() #create and empty vector to store values
  return(as.data.frame(cbind(x[peaks$i],peaks$y.hat[peaks$i])))
  
  } else { }
  plot(x,y, type="l", col="white", xlab = "Time (Sec)", ylab="V", xlim=c(40,50)) # ylim=c(-4,6),
  lines(x, peaks$y.hat,  lwd=2)
  points(x[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=15, cex=1)
  position_peaks <- c() #create and empty vector to store values
  return(as.list(x[peaks$i]))
  
  
}
position_peaks<-peaks(800,200)


smoothing2<- function(x, y, w=60, ...) {   #smoothing data
  require(zoo)
  n <- length(y)
  y.smooth<-data.detrend
  y.max <- rollapply(zoo(y.smooth), 1*w+100, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+100:w)]
  i.max <- which(delta <= 0) + w 
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}



pwave <- function(w,span){
  pwave <-smoothing2(x,y,w=w,span=span)
  position_pwaves <- c() #create and empty vector to store values
  return(as.list(x[pwave$i]))
} 


position_pwaves<-pwave(500,200)
############################################ Finding values the minimum prepeak value of X 

wave <- function(w,span){
  wave <-smoothing2(x,y,w=w,span=span)
  
  z2 <- NULL
  get_prepeaks <- function(y){
    
    z2 <- NULL
    for(i in 2:length(y)){
      sub_ddt <- ddt[seq(as.numeric(y[i-1]) , as.numeric(y[i]), 1) , ]
      sub_zero <- which(abs(sub_ddt$x) <0.008)
      if(length(sub_zero) >0) z2[i-1] <- sub_ddt$time[tail(sub_zero,1)]
    }
    return(z2)
    
  }
  
  z3 <- NULL
  get_prepeaks2 <- function(y){
    
    z3 <- NULL
    for(i in 2:length(y)){
      sub_ddt <- ddt[seq(as.numeric(y[i-1]) , as.numeric(y[i]), 1) , ]
      sub_zero2 <- which(abs(sub_ddt$x) <0.008)
      if(length(sub_zero2) >0) z3[i-1] <- sub_ddt$time[head(sub_zero2,1)]
    }
    return(z3)
    
  }
  
  
  plot(x,y, type="l", col="white", xlab = "Number of Observations",xlim=c(29,30), ylab="V", ylim=c(-4,6))
  lines(x, wave$y.hat,  lwd=2)
  #abline(h=0, col="blue", lwd=2, lty=2)
  points(x[wave$i], wave$y.hat[wave$i], col="Blue", pch=15, cex=1)
  temp_PR <- cbind(get_prepeaks(position_pwaves), 0)
  temp_QRS <- cbind(get_prepeaks2(position_pwaves), 0)
  points(temp_PR, col="Red", pch=15, cex=1)
  points(temp_QRS, col="Green", pch=15, cex=1)
  ##position_prepeaks <- c()
  #return(as.data.frame(cbind(temp_PR, temp_QRS)))
  
}
############################################ Identify the PQ intervals and Calculate HR
PQ_distance<-wave(500, 200)
position_peaks<- as.data.frame(position_peaks)
HR<- length(position_peaks)*10000/length(x)*60
HR<- as.data.frame(HR) 

initial<-1
last<-600

PQ_sample<-do.call(rbind, Map(data.frame, A=PQ_distance$V1[c(initial:last)], B=position_pwaves[c(initial:last)], C=PQ_distance$V3[c(initial:last)]))
colnames(PQ_sample)<-c("PQ", "P&Rwaves", "End of P and end of QRS")

dist_peaks_valls <- PQ_sample %>% 
  mutate (PQ_value= as.numeric(PQ - lag(PQ, default = first(PQ)))) %>%
  mutate(PR_value=as.numeric(`P&Rwaves`-(lag(PQ))+lag(PQ_value))) %>%
  mutate(QRS_value=as.numeric(`End of P and end of QRS`-(lag(PQ))))



startVal <- t(as.data.frame(rep(0, ncol(dist_peaks_valls))))  #Set start value to 0
colnames(startVal) <- colnames(dist_peaks_valls)
dist_peaks_valls_mod <- as.data.frame(rbind(startVal, dist_peaks_valls))
dist_peaks_valls_mod<-na.omit(dist_peaks_valls_mod)

###### DEFINE A SUBSET OF PEAKS TO TO WORK AND TO PLOT (SUBSET)
PQ_int<- as.matrix(dist_peaks_valls_mod$PQ_value[which(dist_peaks_valls_mod$PQ_value <=750 & dist_peaks_valls_mod$PQ_value >= 150)])
PR_int<- as.matrix(dist_peaks_valls_mod$PR_value[which(dist_peaks_valls_mod$PR_value <=750 & dist_peaks_valls_mod$PQ_value >= 150)])
QRS_int<- as.matrix(dist_peaks_valls_mod$QRS_value[which(dist_peaks_valls_mod$QRS_value >=100 & dist_peaks_valls_mod$QRS_value <= 250)])
Pwave_int<- as.matrix(dist_peaks_valls_mod$QRS_value[which(dist_peaks_valls_mod$QRS_value <100)])

PQ_int2<- apply(PQ_int, 2, function(x) c(median1=median(x), SD=sd(x)))/FS #gives the PQ interval in seconds (0.05 seconds=50 ms)
PR_int2<- apply(PR_int, 2, function(x) c(median1=median(x), SD=sd(x)))/FS #gives the PQ interval in seconds (0.05 seconds=50 ms)
QRS_int2<- apply(QRS_int, 2, function(x) c(median1=median(x), SD=sd(x)))/FS #gives the PQ interval in seconds (0.05 seconds=50 ms)
Pwave_int2<- apply(Pwave_int, 2, function(x) c(median1=mean(x), SD=sd(x)))/FS #gives the PQ interval in seconds (0.05 seconds=50 ms)


PQ_PR_Intervalls<-as.data.frame(cbind(PQ_int2, PR_int2, QRS_int2, Pwave_int2))
colnames(PQ_PR_Intervalls)<- c("PQ_Interval", "PR_Interval", "QRSp_Interval", "Pwave")
PQ_PR_Intervalls<-cbind(PQ_PR_Intervalls,HR)
PQ_PR_Intervalls
#setwd("/Volumes/Data/Documents HD (Samsung)/Documents/PhD/R training/Codes/EKG_measurement")
#write.csv(PQ_PR_Intervalls,file=paste(N.1, ".1h.csv"))


