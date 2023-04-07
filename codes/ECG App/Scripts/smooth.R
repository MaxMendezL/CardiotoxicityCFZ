
smoothing <- function(x, y, w=5, ...) { #smoothing data
  require(zoo)
  n <- length(y)
  y.smooth<-data.detrend
  y.max <- rollapply(zoo(y.smooth), 2*w+5, max, align="center")
  y.min <- rollapply(zoo(y.smooth), 2*w+1, min, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+5-1:w)]
  i.max <- which(delta <= 0) + w 
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

smoothing2 <- function(x2, y2, w=5, ...) { #smoothing data
  require(zoo)
  n <- length(y2)
  y.smooth<-data.detrend_BTZ
  y.max <- rollapply(zoo(y.smooth), 2*w, max, align="center")
  y.min <- rollapply(zoo(y.smooth), 2*w+1, min, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+5-1:w)]
  i.max <- which(delta <= 0) + w 
  list(x2=x2[i.max], i=i.max, y.hat=y.smooth)
}

smoothing3 <- function(x3, y3, w=5, ...) { #smoothing data
  require(zoo)
  n <- length(y3)
  y.smooth<-data.detrend_CFZ
  y.max <- rollapply(zoo(y.smooth), 2*w+5, max, align="center")
  y.min <- rollapply(zoo(y.smooth), 2*w+1, min, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+5-1:w)]
  i.max <- which(delta <= 0) + w 
  list(x3=x3[i.max], i=i.max, y.hat=y.smooth)
}

smoothing4 <- function(x4, y4, w=5, ...) { #smoothing data
  require(zoo)
  n <- length(y4)
  y.smooth<-data.detrend_ATRA
  y.max <- rollapply(zoo(y.smooth), 2*w+5, max, align="center")
  y.min <- rollapply(zoo(y.smooth), 2*w+1, min, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+5-1:w)]
  i.max <- which(delta <= 0) + w 
  list(x4=x4[i.max], i=i.max, y.hat=y.smooth)
}

smoothing5 <- function(x5, y5, w=5, ...) { #smoothing data
  require(zoo)
  n <- length(y5)
  y.smooth<-data.detrend_CFZATRA
  y.max <- rollapply(zoo(y.smooth), 2*w+5, max, align="center")
  y.min <- rollapply(zoo(y.smooth), 2*w+1, min, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+5-1:w)]
  i.max <- which(delta <= 0) + w 
  list(x5=x5[i.max], i=i.max, y.hat=y.smooth)
}
