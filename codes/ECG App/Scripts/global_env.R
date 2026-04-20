library(zoo)

# ---------------------------
# Helper: read csv safely
# ---------------------------
read_csv_safe <- function(file, sep = ",", header = TRUE, ...) {
  if (!file.exists(file)) {
    stop(sprintf("Required file not found: %s", file), call. = FALSE)
  }
  read.csv(file, sep = sep, header = header, ...)
}

# ---------------------------
# Helper: read zoo if present
# ---------------------------
read_zoo_optional <- function(file, ...) {
  if (!file.exists(file)) {
    return(NULL)
  }
  read.zoo(file, ...)
}

# ---------------------------
# Example raw input data
# ---------------------------
Control   <- read_csv_safe("Control.csv",   sep = ";", header = FALSE)
BTZ       <- read_csv_safe("BTZ.csv",       sep = ";", header = FALSE)
CFZ       <- read_csv_safe("CFZ.csv",       sep = ";", header = FALSE)
ATRA      <- read_csv_safe("ATRA.csv",      sep = ";", header = FALSE)
CFZATRA   <- read_csv_safe("CFZATRA.csv",   sep = ";", header = FALSE)

# ---------------------------
# Main detrended data frames
# ---------------------------
ddt         <- read_csv_safe("ddt_CTRL.csv",      sep = ",", header = TRUE)
ddt_BTZ     <- read_csv_safe("ddt_BTZ.csv",       sep = ",", header = TRUE)
ddt_CFZ     <- read_csv_safe("ddt_CFZ.csv",       sep = ",", header = TRUE)
ddt_ATRA    <- read_csv_safe("ddt_ATRA.csv",      sep = ",", header = TRUE)
ddt_CFZATRA <- read_csv_safe("ddt_CFZATRA.csv",   sep = ",", header = TRUE)

# ---------------------------
# Validate expected columns
# ---------------------------
validate_ddt <- function(dat, name) {
  required_cols <- c("time", "data.detrend")
  missing_cols <- setdiff(required_cols, names(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Dataset '%s' is missing required columns: %s",
        name,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

validate_ddt(ddt, "ddt_CTRL")
validate_ddt(ddt_BTZ, "ddt_BTZ")
validate_ddt(ddt_CFZ, "ddt_CFZ")
validate_ddt(ddt_ATRA, "ddt_ATRA")
validate_ddt(ddt_CFZATRA, "ddt_CFZATRA")

# ---------------------------
# Detrended signals
# Prefer standalone files if they exist,
# otherwise fall back to ddt_*$data.detrend
# ---------------------------
tmp_ctrl <- read_zoo_optional("data.detrend_CTRL.csv")
tmp_btz  <- read_zoo_optional("data.detrend_BTZ.csv")
tmp_cfz  <- read_zoo_optional("data.detrend_CFZ.csv")
tmp_atra <- read_zoo_optional("data.detrend_ATRA.csv")
tmp_cfza <- read_zoo_optional("data.detrend_CFZATRA.csv")

data.detrend         <- if (!is.null(tmp_ctrl)) as.numeric(tmp_ctrl) else as.numeric(ddt$data.detrend)
data.detrend_BTZ     <- if (!is.null(tmp_btz))  as.numeric(tmp_btz)  else as.numeric(ddt_BTZ$data.detrend)
data.detrend_CFZ     <- if (!is.null(tmp_cfz))  as.numeric(tmp_cfz)  else as.numeric(ddt_CFZ$data.detrend)
data.detrend_ATRA    <- if (!is.null(tmp_atra)) as.numeric(tmp_atra) else as.numeric(ddt_ATRA$data.detrend)
data.detrend_CFZATRA <- if (!is.null(tmp_cfza)) as.numeric(tmp_cfza) else as.numeric(ddt_CFZATRA$data.detrend)

# ---------------------------
# Precomputed ECG summary values
# ---------------------------
ECGvals_CTRL     <- read_csv_safe("ECGvals_CTRL.csv")
ECGvals_BTZ      <- read_csv_safe("ECGvals_BTZ.csv")
ECGvals_CFZ      <- read_csv_safe("ECGvals_CFZ.csv")
ECGvals_ATRA     <- read_csv_safe("ECGvals_ATRA.csv")
ECGvals_CFZATRA  <- read_csv_safe("ECGvals_CFZATRA.csv")

# ---------------------------
# Sampling frequency
# ---------------------------
FS <- 10000

# ---------------------------
# Backward-compatible vectors
# ---------------------------
x  <- ddt$time
y  <- ddt$data.detrend

x2 <- ddt_BTZ$time
y2 <- ddt_BTZ$data.detrend

x3 <- ddt_CFZ$time
y3 <- ddt_CFZ$data.detrend

x4 <- ddt_ATRA$time
y4 <- ddt_ATRA$data.detrend

x5 <- ddt_CFZATRA$time
y5 <- ddt_CFZATRA$data.detrend

# ---------------------------
# Unified object for cleaner future code
# ---------------------------
ecg_data <- list(
  Control = list(
    raw = Control,
    ddt = ddt,
    detrend = data.detrend,
    ecgvals = ECGvals_CTRL,
    x = x,
    y = y
  ),
  BTZ = list(
    raw = BTZ,
    ddt = ddt_BTZ,
    detrend = data.detrend_BTZ,
    ecgvals = ECGvals_BTZ,
    x = x2,
    y = y2
  ),
  CFZ = list(
    raw = CFZ,
    ddt = ddt_CFZ,
    detrend = data.detrend_CFZ,
    ecgvals = ECGvals_CFZ,
    x = x3,
    y = y3
  ),
  ATRA = list(
    raw = ATRA,
    ddt = ddt_ATRA,
    detrend = data.detrend_ATRA,
    ecgvals = ECGvals_ATRA,
    x = x4,
    y = y4
  ),
  CFZATRA = list(
    raw = CFZATRA,
    ddt = ddt_CFZATRA,
    detrend = data.detrend_CFZATRA,
    ecgvals = ECGvals_CFZATRA,
    x = x5,
    y = y5
  )
)
