wavesUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    selectInput(
      ns("Input"),
      "P waves",
      choices = c("Control", "BTZ", "CFZ", "ATRA", "CFZATRA")
    ),
    plotOutput(ns("plot3"))
  )
}


customPlot3 <- function(input, output, session) {
  
  output$plot3 <- renderPlot({
    
    req(input$Input)
    
    # ---------------------------------
    # Condition-specific data + tuning
    # ---------------------------------
    cfg_list <- list(
      Control = list(
        ddt = ddt,
        detrend = data.detrend,
        xlim = c(295000, 305000),
        r_min_distance = 900,
        r_polarity = "absolute",
        r_quantile = 0.88,
        p_lookback = 500,
        p_guard = 100,
        p_polarity = "positive",
        p_frac = 0.25
      ),
      BTZ = list(
        ddt = ddt_BTZ,
        detrend = data.detrend_BTZ,
        xlim = c(20000, 30000),
        r_min_distance = 1100,
        r_polarity = "positive",
        r_quantile = 0.90,
        p_lookback = 500,
        p_guard = 100,
        p_polarity = "positive",
        p_frac = 0.25
      ),
      CFZ = list(
        ddt = ddt_CFZ,
        detrend = data.detrend_CFZ,
        xlim = c(295000, 305000),
        r_min_distance = 1300,
        r_polarity = "absolute",
        r_quantile = 0.88,
        p_lookback = 550,
        p_guard = 110,
        p_polarity = "positive",
        p_frac = 0.25
      ),
      ATRA = list(
        ddt = ddt_ATRA,
        detrend = data.detrend_ATRA,
        xlim = c(295000, 305000),
        r_min_distance = 1100,
        r_polarity = "absolute",
        r_quantile = 0.92,
        p_lookback = 500,
        p_guard = 100,
        p_polarity = "positive",
        p_frac = 0.25
      ),
      CFZATRA = list(
        ddt = ddt_CFZATRA,
        detrend = data.detrend_CFZATRA,
        xlim = c(290000, 300000),
        r_min_distance = 1300,
        r_polarity = "absolute",
        r_quantile = 0.88,
        p_lookback = 550,
        p_guard = 110,
        p_polarity = "positive",
        p_frac = 0.25
      )
    )
    
    cfg <- cfg_list[[input$Input]]
    dat <- cfg$ddt
    
    if (!all(c("time", "data.detrend") %in% names(dat))) {
      stop(sprintf("Dataset '%s' must contain columns 'time' and 'data.detrend'.", input$Input), call. = FALSE)
    }
    
    x <- dat$time
    y <- as.numeric(dat$data.detrend)
    
    # ---------------------------------
    # Light smoothing
    # ---------------------------------
    smooth_signal <- function(y, k = 5) {
      y_sm <- as.numeric(stats::filter(y, rep(1 / k, k), sides = 2))
      idx_na <- which(is.na(y_sm))
      if (length(idx_na) > 0) {
        y_sm[idx_na] <- y[idx_na]
      }
      y_sm
    }
    
    y_sm <- smooth_signal(y, k = 5)
    
    # ---------------------------------
    # Polarity-aware R/QRS detector
    # ---------------------------------
    detect_r_peaks <- function(y, min_height = NULL, min_distance = 800,
                               polarity = c("auto", "positive", "negative", "absolute")) {
      polarity <- match.arg(polarity)
      
      if (polarity == "auto") {
        pos_q <- stats::quantile(y, 0.98, na.rm = TRUE)
        neg_q <- abs(stats::quantile(y, 0.02, na.rm = TRUE))
        polarity <- if (neg_q > pos_q) "negative" else "positive"
      }
      
      y_work <- switch(
        polarity,
        positive = y,
        negative = -y,
        absolute = abs(y)
      )
      
      if (is.null(min_height)) {
        min_height <- stats::quantile(y_work, 0.90, na.rm = TRUE)
      }
      
      cand <- which(
        y_work > min_height &
          c(FALSE, diff(y_work) > 0) &
          c(diff(y_work) <= 0, FALSE)
      )
      
      if (length(cand) == 0) {
        return(integer(0))
      }
      
      keep <- cand[1]
      
      if (length(cand) > 1) {
        for (idx in cand[-1]) {
          if ((idx - tail(keep, 1)) >= min_distance) {
            keep <- c(keep, idx)
          } else {
            if (y_work[idx] > y_work[tail(keep, 1)]) {
              keep[length(keep)] <- idx
            }
          }
        }
      }
      
      keep
    }
    
    # ---------------------------------
    # P-wave detector: one candidate before each R
    # ---------------------------------
    detect_p_before_r <- function(y, r_idx, lookback = 350, guard = 80,
                                  p_polarity = c("positive", "absolute")) {
      p_polarity <- match.arg(p_polarity)
      
      p_idx <- integer(0)
      
      if (length(r_idx) == 0) {
        return(p_idx)
      }
      
      for (r in r_idx) {
        start <- max(1, r - lookback)
        end   <- max(start, r - guard)
        
        if (end <= start) next
        
        seg <- y[start:end]
        
        loc <- if (p_polarity == "positive") {
          which.max(seg)
        } else {
          which.max(abs(seg))
        }
        
        if (length(loc) == 1 && !is.na(loc)) {
          p_idx <- c(p_idx, start + loc - 1)
        }
      }
      
      unique(p_idx)
    }
    
    # ---------------------------------
    # Approximate P onset / offset
    # ---------------------------------
    get_p_bounds <- function(y, p_idx, frac = 0.25, max_back = 120, max_forward = 120) {
      onset <- integer(0)
      offset <- integer(0)
      
      if (length(p_idx) == 0) {
        return(list(onset = onset, offset = offset))
      }
      
      for (p in p_idx) {
        amp <- abs(y[p])
        thr <- frac * amp
        
        left <- p
        step <- 0
        while (left > 1 && step < max_back && abs(y[left]) > thr) {
          left <- left - 1
          step <- step + 1
        }
        
        right <- p
        step <- 0
        while (right < length(y) && step < max_forward && abs(y[right]) > thr) {
          right <- right + 1
          step <- step + 1
        }
        
        onset <- c(onset, left)
        offset <- c(offset, right)
      }
      
      list(onset = onset, offset = offset)
    }
    
    # ---------------------------------
    # Detect R/QRS peaks
    # ---------------------------------
    y_for_threshold <- switch(
      cfg$r_polarity,
      positive = y_sm,
      negative = -y_sm,
      absolute = abs(y_sm)
    )
    
    r_idx <- detect_r_peaks(
      y = y_sm,
      min_height = stats::quantile(y_for_threshold, cfg$r_quantile, na.rm = TRUE),
      min_distance = cfg$r_min_distance,
      polarity = cfg$r_polarity
    )
    
    # ---------------------------------
    # Detect P waves before each R
    # ---------------------------------
    p_idx <- detect_p_before_r(
      y = y_sm,
      r_idx = r_idx,
      lookback = cfg$p_lookback,
      guard = cfg$p_guard,
      p_polarity = cfg$p_polarity
    )
    
    # ---------------------------------
    # Approximate P bounds
    # ---------------------------------
    p_bounds <- get_p_bounds(
      y = y_sm,
      p_idx = p_idx,
      frac = cfg$p_frac,
      max_back = 120,
      max_forward = 120
    )
    
    # ---------------------------------
    # Plot
    # ---------------------------------
    plot(
      x, y,
      lwd = 2,
      type = "l",
      col = "black",
      xlab = "Time (Sec)",
      xlim = cfg$xlim,
      main = "P waves and Intervals",
      bty = "n",
      ylab = "Voltage",
      ylim = c(-4, 4),
      xaxt = "n"
    )
    
    # QRS / R anchors
    if (length(r_idx) > 0) {
      points(x[r_idx], y[r_idx], col = "darkgreen", pch = 17, cex = 1.3)
    }
    
    # P-wave peaks
    if (length(p_idx) > 0) {
      points(x[p_idx], y[p_idx], col = "blue", pch = 19, cex = 1.2)
    }
    
    # P-wave onset
    if (length(p_bounds$onset) > 0) {
      points(x[p_bounds$onset], y[p_bounds$onset], col = "red", pch = 15, cex = 1)
    }
    
    # P-wave offset
    if (length(p_bounds$offset) > 0) {
      points(x[p_bounds$offset], y[p_bounds$offset], col = "red", pch = 15, cex = 1)
    }
  })
  
  return()
}
