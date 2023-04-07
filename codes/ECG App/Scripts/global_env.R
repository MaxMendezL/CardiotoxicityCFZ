#Example Data
Control <- read.csv("Control.csv", sep=";", header=FALSE) #change "," to ";" accordingly
BTZ <- read.csv("BTZ.csv", sep=";", header=FALSE)
CFZ <- read.csv("CFZ.csv", sep=";", header=FALSE)
ATRA <- read.csv("ATRA.csv", sep=";", header=FALSE)
CFZATRA <- read.csv("CFZATRA.csv", sep=";", header=FALSE)

data.detrend <- read.zoo("data.detrend_CTRL.csv")
data.detrend_BTZ <- read.zoo("data.detrend_BTZ.csv")
data.detrend_CFZ <- read.zoo("data.detrend_CFZ.csv")
data.detrend_ATRA <- read.zoo("data.detrend_ATRA.csv")
data.detrend_CFZATRA <- read.zoo("data.detrend_CFZATRA.csv")

ddt <- read.csv("ddt_CTRL.csv", sep=",", header=TRUE) #change "," to ";" accordingly
ddt_BTZ <- read.csv("ddt_BTZ.csv", sep=",", header=TRUE) #change "," to ";" accordingly
ddt_CFZ <- read.csv("ddt_CFZ.csv", sep=",", header=TRUE) #change "," to ";" accordingly
ddt_ATRA <- read.csv("ddt_ATRA.csv", sep=",", header=TRUE) #change "," to ";" accordingly
ddt_CFZATRA <- read.csv("ddt_CFZATRA.csv", sep=",", header=TRUE) #change "," to ";" accordingly

ECGvals_CTRL<- read.csv("ECGvals_CTRL.csv")
ECGvals_BTZ<- read.csv("ECGvals_BTZ.csv")
ECGvals_CFZ<- read.csv("ECGvals_CFZ.csv")
ECGvals_ATRA<- read.csv("ECGvals_ATRA.csv")
ECGvals_CFZATRA<- read.csv("ECGvals_CFZATRA.csv")

FS<-10000
x<-ddt$time
y<-ddt$data.detrend
x2<-ddt_BTZ$time
y2<-ddt_BTZ$data.detrend
x3<-ddt_CFZ$time
y3<-ddt_CFZ$data.detrend
x4<-ddt_ATRA$time
y4<-ddt_ATRA$data.detrend
x5<-ddt_CFZATRA$time
y5<-ddt_CFZATRA$data.detrend


