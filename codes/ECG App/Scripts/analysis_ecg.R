analyze_ecg_signal <- function(group_name) {
  
  cfg_list <- list(
    Control = list(
      ddt = ddt,
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
  
  cfg <- cfg_list[[group_name]]
  dat <- cfg$ddt
  
  x <- dat$time
  y <- as.numeric(dat$data.detrend)
  
  smooth_signal <- function(y, k = 5) {
    y_sm <- as.numeric(stats::filter(y, rep(1 / k, k), sides = 2))
    idx_na <- which(is.na(y_sm))
    if (length(idx_na) > 0) {
      y_sm[idx_na] <- y[idx_na]
    }
    y_sm
  }
  
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
    
    if (length(cand) == 0) return(integer(0))
    
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
  
  detect_p_before_r <- function(y, r_idx, lookback = 350, guard = 80,
                                p_polarity = c("positive", "absolute")) {
    p_polarity <- match.arg(p_polarity)
    
    p_idx <- integer(0)
    if (length(r_idx) == 0) return(p_idx)
    
    for (r in r_idx) {
      start <- max(1, r - lookback)
      end <- max(start, r - guard)
      if (end <= start) next
      
      seg <- y[start:end]
      loc <- if (p_polarity == "positive") which.max(seg) else which.max(abs(seg))
      
      if (length(loc) == 1 && !is.na(loc)) {
        p_idx <- c(p_idx, start + loc - 1)
      }
    }
    
    unique(p_idx)
  }
  
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
  
  y_sm <- smooth_signal(y, k = 5)
  
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
  
  p_idx <- detect_p_before_r(
    y = y_sm,
    r_idx = r_idx,
    lookback = cfg$p_lookback,
    guard = cfg$p_guard,
    p_polarity = cfg$p_polarity
  )
  
  p_bounds <- get_p_bounds(
    y = y_sm,
    p_idx = p_idx,
    frac = cfg$p_frac
  )
  
  # --- Metrics ---
  dt <- 1 / FS
  
  rr_sec <- if (length(r_idx) >= 2) diff(r_idx) * dt else numeric(0)
  hr_bpm <- if (length(rr_sec) > 0) 60 / rr_sec else numeric(0)
  
  pr_sec <- if (length(p_bounds$onset) == length(r_idx) && length(r_idx) > 0) {
    (r_idx - p_bounds$onset) * dt
  } else {
    numeric(0)
  }
  
  pq_sec <- if (length(p_idx) == length(r_idx) && length(r_idx) > 0) {
    (r_idx - p_idx) * dt
  } else {
    numeric(0)
  }
  
  pwave_sec <- if (length(p_bounds$onset) == length(p_bounds$offset) && length(p_bounds$onset) > 0) {
    (p_bounds$offset - p_bounds$onset) * dt
  } else {
    numeric(0)
  }
  
  # Placeholder QRS width estimate around anchor
  qrsp_sec <- if (length(r_idx) > 0) {
    rep(0.01, length(r_idx))
  } else {
    numeric(0)
  }
  
  make_row <- function(label, values) {
    data.frame(
      Value = label,
      Median = if (length(values) > 0) round(stats::median(values, na.rm = TRUE), 4) else NA_real_,
      STDev = if (length(values) > 1) round(stats::sd(values, na.rm = TRUE), 4) else 0,
      stringsAsFactors = FALSE
    )
  }
  
  summary_table <- rbind(
    make_row("PR (sec)", pr_sec),
    make_row("PQ (sec)", pq_sec),
    make_row("QRSp (sec)", qrsp_sec),
    make_row("Pwave (sec)", pwave_sec),
    make_row("Heart Rate (beats/min)", hr_bpm)
  )
  
  list(
    x = x,
    y = y,
    y_sm = y_sm,
    xlim = cfg$xlim,
    r_idx = r_idx,
    p_idx = p_idx,
    p_onset = p_bounds$onset,
    p_offset = p_bounds$offset,
    table = summary_table
  )
}
