smoothing <- function(x, y, w = 5, width_extra = 5) {
  n <- length(y)
  
  if (length(x) != n) {
    stop("x and y must have the same length.", call. = FALSE)
  }
  
  y.smooth <- as.numeric(y)
  
  y.max <- zoo::rollapply(
    zoo::zoo(y.smooth),
    width = 2 * w + width_extra,
    FUN = max,
    align = "center",
    fill = NA
  )
  
  keep <- which(!is.na(y.max))
  delta <- y.max[keep] - y.smooth[keep]
  i.max <- keep[delta <= 0]
  
  list(
    x = x[i.max],
    i = i.max,
    y.hat = y.smooth
  )
}

smoothing2 <- function(x, y, w = 5, ...) {
  smoothing(x, y, w = w, width_extra = 0)
}

smoothing3 <- function(x, y, w = 5, ...) {
  smoothing(x, y, w = w, width_extra = 5)
}

smoothing4 <- function(x, y, w = 5, ...) {
  smoothing(x, y, w = w, width_extra = 5)
}

smoothing5 <- function(x, y, w = 5, ...) {
  smoothing(x, y, w = w, width_extra = 5)
}
