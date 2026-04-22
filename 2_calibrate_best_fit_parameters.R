# ================================================================
# 2_calibrate_best_fit_parameters.R
#
# Calibration workflow to retrieve one best-fit parameter set for the
# Tree Water Imbalance and Storage Tracker (TWIST) module
#
# This script:
# 1) Loads TWIST functions
# 2) Defines calibration settings
# 3) Reads input and validation data
# 4) Prepares model inputs and validation targets
# 5) Builds a scaled calibration objective
# 6) Runs a multi-start optimization
# 7) Evaluates the best-fit parameter set
# 8) Writes two output tables
#
# The script is intentionally minimal and only returns:
# - one table with the best-fit TWIST parameter set
# - one table with calibration settings and performance metrics
#
# Note: All water-related variables (E, TWD, W) must share the same unit.
# ================================================================

# 1) Load functions ----
source("0_TWIST_functions.R")

# Load package used to write the metadata workbook
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop(
    "Package 'openxlsx' is required to write the metadata workbook. ",
    "Please install it first with install.packages('openxlsx')."
  )
}

# 2) Define calibration settings ----

# Input file paths
input_file <- "data/TWIST_input_data.rds"   
validation_file <- "data/TWIST_validation_data.rds"  

# Calibration years
calibration_years <- c(2018)

# TWIST parameter bounds
lower <- c(F_E = 0.0, F_TWD = 0.0, F_theta = 0.01)
upper <- c(F_E = 1.0, F_TWD = 1.0, F_theta = 1.0)

# Reference parameter set used to scale the hourly and daily RMSE terms
ref_par <- c(F_E = 0.6, F_TWD = 0.3, F_theta = 0.7)

# Multi-start optimization settings
n_starts <- 10
optim_control <- list(maxit = 400)
random_seed <- 123

# Objective weights
# Final objective = weighted hourly RMSE + weighted daily RMSE,
# with both components scaled by the reference parameter set.
w_hourly <- 1
w_daily  <- 1

# Parallel settings
use_parallel <- TRUE
n_cores <- max(1, parallel::detectCores() - 1)

# Water pool parameters
params_water_pool <- list(
  rho_sat = 1.07,
  rho_dry = 0.58
)

# Column names in the input data
col_names <- list(
  col_time   = "datetime",
  col_E      = "transpiration_l.m2",
  col_theta  = "theta_rel",
  col_m_wood = "m_dry_wood_kg.m2",
  col_W      = "W"
)

# Output directory
out_dir <- "output/twist_best_fit_calibration"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 3) Define helper functions ----

# Root mean square error
rmse_fast <- function(obs, sim) {
  ok <- is.finite(obs) & is.finite(sim)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((sim[ok] - obs[ok])^2))
}

# Mean bias
bias_fast <- function(obs, sim) {
  ok <- is.finite(obs) & is.finite(sim)
  if (!any(ok)) return(NA_real_)
  mean(sim[ok] - obs[ok])
}

# Coefficient of determination based on Pearson correlation
r2_fast <- function(obs, sim) {
  ok <- is.finite(obs) & is.finite(sim)
  if (sum(ok) < 2) return(NA_real_)
  cor(obs[ok], sim[ok])^2
}

# Daily mean helper used for aggregation of hourly TWD
# One value is returned per day.
daily_mean_fast <- function(x, group_factor) {
  tapply(
    x,
    group_factor,
    function(v) {
      if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
    }
  )
}

# Efficient calculation of q99_MDWD over the full period
# This is used to normalize TWD before model-data comparison.
compute_q99_mdwd_fast <- function(twd, dates, hrs) {
  idx_min <- hrs >= 4 & hrs <= 8
  idx_max <- hrs >= 10 & hrs <= 23

  min_by_day <- tapply(
    twd[idx_min],
    dates[idx_min],
    function(v) {
      if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
    }
  )

  max_by_day <- tapply(
    twd[idx_max],
    dates[idx_max],
    function(v) {
      if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    }
  )

  common_days <- intersect(names(min_by_day), names(max_by_day))
  if (length(common_days) == 0L) return(NA_real_)

  mdwd <- max_by_day[common_days] - min_by_day[common_days]
  if (!any(is.finite(mdwd))) return(NA_real_)

  quantile(mdwd, probs = 0.99, na.rm = TRUE, names = FALSE)
}

# Lightweight TWIST runner used inside the optimizer
# This reproduces the same model equations as TWIST_functions.R but avoids
# repeated data-frame overhead during calibration.
run_twist_fast <- function(E, theta_rel, W, F_E, F_TWD, F_theta, TWD_old = 0) {
  n <- length(E)

  TWD <- numeric(n)
  RWC <- numeric(n)

  for (i in seq_len(n)) {
    f_soil_val <- min(theta_rel[i] / F_theta, 1)
    uptake <- (F_E * E[i] + F_TWD * TWD_old) * f_soil_val
    TWD_new <- TWD_old + E[i] - uptake

    TWD[i] <- TWD_new
    RWC[i] <- max(0, (W[i] - TWD_new) / W[i])

    TWD_old <- TWD_new
  }

  list(TWD = TWD, RWC = RWC)
}

# Evaluate one parameter set and return the raw error components
# used to build the calibration objective.
evaluate_parameters_raw <- function(par) {
  F_E     <- par["F_E"]
  F_TWD   <- par["F_TWD"]
  F_theta <- par["F_theta"]

  run_res <- run_twist_fast(
    E = E_vec,
    theta_rel = theta_vec,
    W = W_vec,
    F_E = F_E,
    F_TWD = F_TWD,
    F_theta = F_theta,
    TWD_old = 0
  )

  q99_mdwd <- compute_q99_mdwd_fast(
    twd = run_res$TWD,
    dates = date_vec,
    hrs = hour_vec
  )

  if (!is.finite(q99_mdwd) || q99_mdwd <= 0) {
    return(list(
      hourly_rmse = NA_real_,
      daily_rmse = NA_real_,
      q99_mdwd = q99_mdwd,
      TWD_norm = NULL,
      RWC = NULL
    ))
  }

  twd_norm <- run_res$TWD / q99_mdwd

  hourly_rmse <- rmse_fast(obs_hourly_cal, twd_norm[cal_mask_obs])

  sim_daily_cal <- daily_mean_fast(twd_norm[cal_mask_obs], cal_day_factor)
  daily_rmse <- rmse_fast(obs_daily_cal, sim_daily_cal)

  list(
    hourly_rmse = hourly_rmse,
    daily_rmse = daily_rmse,
    q99_mdwd = q99_mdwd,
    TWD_norm = twd_norm,
    RWC = run_res$RWC
  )
}

# Evaluate one parameter set and return the scaled objective.
evaluate_parameters <- function(par) {
  raw_res <- evaluate_parameters_raw(par)

  if (!is.finite(raw_res$hourly_rmse) || !is.finite(raw_res$daily_rmse)) {
    return(list(
      objective = Inf,
      hourly_rmse = raw_res$hourly_rmse,
      daily_rmse = raw_res$daily_rmse,
      q99_mdwd = raw_res$q99_mdwd,
      TWD_norm = raw_res$TWD_norm,
      RWC = raw_res$RWC
    ))
  }

  objective <-
    w_hourly * (raw_res$hourly_rmse / scale_hourly) +
    w_daily  * (raw_res$daily_rmse  / scale_daily)

  list(
    objective = objective,
    hourly_rmse = raw_res$hourly_rmse,
    daily_rmse = raw_res$daily_rmse,
    q99_mdwd = raw_res$q99_mdwd,
    TWD_norm = raw_res$TWD_norm,
    RWC = raw_res$RWC
  )
}

# Objective function passed to optim()
objective_function <- function(par) {
  if (any(!is.finite(par))) return(Inf)
  if (par["F_E"]     < lower["F_E"]     || par["F_E"]     > upper["F_E"])     return(Inf)
  if (par["F_TWD"]   < lower["F_TWD"]   || par["F_TWD"]   > upper["F_TWD"])   return(Inf)
  if (par["F_theta"] < lower["F_theta"] || par["F_theta"] > upper["F_theta"]) return(Inf)

  evaluate_parameters(par)$objective
}

# Compute summary metrics for one subset of the time series.
compute_metrics_for_mask <- function(obs, sim, mask_obs, day_factor) {
  obs_hourly <- obs[mask_obs]
  sim_hourly <- sim[mask_obs]

  obs_daily <- daily_mean_fast(obs_hourly, day_factor)
  sim_daily <- daily_mean_fast(sim_hourly, day_factor)

  list(
    n_hourly    = sum(is.finite(obs_hourly) & is.finite(sim_hourly)),
    n_daily     = sum(is.finite(obs_daily) & is.finite(sim_daily)),
    rmse_hourly = rmse_fast(obs_hourly, sim_hourly),
    rmse_daily  = rmse_fast(obs_daily, sim_daily),
    bias_hourly = bias_fast(obs_hourly, sim_hourly),
    bias_daily  = bias_fast(obs_daily, sim_daily),
    r2_hourly   = r2_fast(obs_hourly, sim_hourly),
    r2_daily    = r2_fast(obs_daily, sim_daily)
  )
}

# Run one bounded optimization from one start point.
optimize_one_start <- function(i, start_points_mat) {
  p0 <- start_points_mat[i, ]

  fit <- tryCatch(
    optim(
      par = p0,
      fn = objective_function,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = optim_control
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  data.frame(
    start_id = i,
    convergence = fit$convergence,
    objective = fit$value,
    F_E = fit$par[1],
    F_TWD = fit$par[2],
    F_theta = fit$par[3],
    stringsAsFactors = FALSE
  )
}

# 4) Read and prepare data ----

# Read module input and validation data
# The exact file names and paths can be adapted later.
df_input <- readRDS(input_file)
validation_data <-  readRDS(validation_file)

# Keep only the required columns
required_input_cols <- c(
  col_names$col_time,
  col_names$col_E,
  col_names$col_theta,
  col_names$col_m_wood
)

missing_input_cols <- setdiff(required_input_cols, names(df_input))
if (length(missing_input_cols) > 0) {
  stop(
    "The following required input columns are missing: ",
    paste(missing_input_cols, collapse = ", ")
  )
}

df_input <- df_input[, required_input_cols]
validation_data <- validation_data[, c("datetime", "mean_TWD_norm")]

# Estimate water pool size W for each timestep
# This uses the same function as in example_run.R.
df_input[[col_names$col_W]] <- compute_water_pool(
  m_wood_dry = df_input[[col_names$col_m_wood]],
  params_water_pool = params_water_pool
)

# Extract time information once for efficient reuse
# These vectors are used repeatedly during calibration and evaluation.
time_vec <- df_input[[col_names$col_time]]
date_vec <- as.Date(time_vec)
hour_vec <- as.POSIXlt(time_vec)$hour
year_vec <- as.POSIXlt(time_vec)$year + 1900

# Extract model inputs as plain vectors for faster optimization
E_vec <- df_input[[col_names$col_E]]
theta_vec <- df_input[[col_names$col_theta]]
W_vec <- df_input[[col_names$col_W]]

# Align validation data to the TWIST model time axis
# Validation is assumed to provide normalized TWD observations.
validation_match <- match(time_vec, validation_data$datetime)
obs_full <- validation_data$mean_TWD_norm[validation_match]

# Define calibration and evaluation subsets
cal_mask <- year_vec %in% calibration_years
eval_mask <- !(year_vec %in% calibration_years)

# Keep only rows with observed validation data for scoring
cal_mask_obs <- cal_mask & is.finite(obs_full)
eval_mask_obs <- eval_mask & is.finite(obs_full)

if (sum(cal_mask_obs) < 2) {
  stop("Calibration period does not contain enough valid observations.")
}

# Precompute daily grouping for the calibration and evaluation periods
cal_day_factor <- factor(date_vec[cal_mask_obs])
eval_day_factor <- factor(date_vec[eval_mask_obs])

# Precompute observed validation series used in the objective function
obs_hourly_cal <- obs_full[cal_mask_obs]
obs_daily_cal <- daily_mean_fast(obs_hourly_cal, cal_day_factor)

# 5) Build objective scaling ----

# The objective combines hourly and daily RMSE.
# Both terms are scaled with the reference parameter set so that neither
# term dominates only because of its magnitude.
ref_eval <- evaluate_parameters_raw(ref_par)

if (!is.finite(ref_eval$hourly_rmse) || !is.finite(ref_eval$daily_rmse)) {
  stop("Reference parameter set did not produce valid calibration errors.")
}

scale_hourly <- ref_eval$hourly_rmse
scale_daily  <- ref_eval$daily_rmse

if (scale_hourly <= 0 || scale_daily <= 0) {
  stop("Objective scaling factors must be greater than zero.")
}

# 6) Run multi-start optimization ----

# Generate start points within the bounded parameter space.
# The reference parameter set is included as the first start.
set.seed(random_seed)

start_points <- rbind(
  ref_par,
  cbind(
    runif(n_starts - 1, lower["F_E"], upper["F_E"]),
    runif(n_starts - 1, lower["F_TWD"], upper["F_TWD"]),
    runif(n_starts - 1, lower["F_theta"], upper["F_theta"])
  )
)
colnames(start_points) <- c("F_E", "F_TWD", "F_theta")

cat("\nStarting multi-start optimization ...\n")
cat("Starts:", n_starts, "\n")
cat("Parallel:", use_parallel, "\n")
cat("Requested cores:", n_cores, "\n")

optimization_time <- system.time({
  if (use_parallel && n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c(
        "E_vec", "theta_vec", "W_vec",
        "date_vec", "hour_vec",
        "obs_hourly_cal", "obs_daily_cal",
        "cal_mask_obs", "cal_day_factor",
        "lower", "upper", "optim_control",
        "scale_hourly", "scale_daily",
        "w_hourly", "w_daily",
        "compute_q99_mdwd_fast",
        "run_twist_fast",
        "rmse_fast", "daily_mean_fast",
        "evaluate_parameters_raw",
        "evaluate_parameters",
        "objective_function",
        "optimize_one_start",
        "start_points"
      ),
      envir = environment()
    )

    optimization_list <- parallel::parLapply(
      cl,
      X = seq_len(nrow(start_points)),
      fun = optimize_one_start,
      start_points_mat = start_points
    )
  } else {
    optimization_list <- lapply(
      seq_len(nrow(start_points)),
      optimize_one_start,
      start_points_mat = start_points
    )
  }
})

optimization_list <- optimization_list[!vapply(optimization_list, is.null, logical(1))]
if (length(optimization_list) == 0L) {
  stop("Optimization did not return any valid result.")
}

optimization_summary <- do.call(rbind, optimization_list)
optimization_summary <- optimization_summary[order(optimization_summary$objective), ]
row.names(optimization_summary) <- NULL

# 7) Evaluate the best-fit parameter set ----

best_par <- as.numeric(optimization_summary[1, c("F_E", "F_TWD", "F_theta")])
names(best_par) <- c("F_E", "F_TWD", "F_theta")

best_eval <- evaluate_parameters(best_par)
if (is.null(best_eval$TWD_norm) || !is.finite(best_eval$objective)) {
  stop("Best-fit parameter set did not produce valid model output.")
}

metrics_cal <- compute_metrics_for_mask(
  obs = obs_full,
  sim = best_eval$TWD_norm,
  mask_obs = cal_mask_obs,
  day_factor = cal_day_factor
)

metrics_eval <- compute_metrics_for_mask(
  obs = obs_full,
  sim = best_eval$TWD_norm,
  mask_obs = eval_mask_obs,
  day_factor = eval_day_factor
)

# Create a compact parameter table that can be used directly in example_run.R.
best_fit_parameters <- data.frame(
  F_E = best_par["F_E"],
  F_TWD = best_par["F_TWD"],
  F_theta = best_par["F_theta"],
  stringsAsFactors = FALSE
)

# Create a compact metadata table describing the calibration setup
# and the resulting model performance.
best_fit_metadata <- data.frame(
  calibration_years = paste(calibration_years, collapse = ","),
  n_starts = n_starts,
  use_parallel = use_parallel,
  n_cores = if (use_parallel) n_cores else 1,
  objective_weight_hourly = w_hourly,
  objective_weight_daily = w_daily,
  lower_F_E = lower["F_E"],
  upper_F_E = upper["F_E"],
  lower_F_TWD = lower["F_TWD"],
  upper_F_TWD = upper["F_TWD"],
  lower_F_theta = lower["F_theta"],
  upper_F_theta = upper["F_theta"],
  reference_F_E = ref_par["F_E"],
  reference_F_TWD = ref_par["F_TWD"],
  reference_F_theta = ref_par["F_theta"],
  scale_hourly = scale_hourly,
  scale_daily = scale_daily,
  objective = best_eval$objective,
  q99_mdwd = best_eval$q99_mdwd,
  cal_n_hourly = metrics_cal$n_hourly,
  cal_n_daily = metrics_cal$n_daily,
  cal_rmse_hourly = metrics_cal$rmse_hourly,
  cal_rmse_daily = metrics_cal$rmse_daily,
  cal_bias_hourly = metrics_cal$bias_hourly,
  cal_bias_daily = metrics_cal$bias_daily,
  cal_r2_hourly = metrics_cal$r2_hourly,
  cal_r2_daily = metrics_cal$r2_daily,
  eval_n_hourly = metrics_eval$n_hourly,
  eval_n_daily = metrics_eval$n_daily,
  eval_rmse_hourly = metrics_eval$rmse_hourly,
  eval_rmse_daily = metrics_eval$rmse_daily,
  eval_bias_hourly = metrics_eval$bias_hourly,
  eval_bias_daily = metrics_eval$bias_daily,
  eval_r2_hourly = metrics_eval$r2_hourly,
  eval_r2_daily = metrics_eval$r2_daily,
  runtime_elapsed_sec = unname(optimization_time["elapsed"]),
  stringsAsFactors = FALSE
)

# Create a second table that explains all metadata columns.
# This sheet is intended for direct inspection in Excel.
metadata_column_description <- data.frame(
  column_name = c(
    "calibration_years",
    "n_starts",
    "use_parallel",
    "n_cores",
    "objective_weight_hourly",
    "objective_weight_daily",
    "lower_F_E",
    "upper_F_E",
    "lower_F_TWD",
    "upper_F_TWD",
    "lower_F_theta",
    "upper_F_theta",
    "reference_F_E",
    "reference_F_TWD",
    "reference_F_theta",
    "scale_hourly",
    "scale_daily",
    "objective",
    "q99_mdwd",
    "cal_n_hourly",
    "cal_n_daily",
    "cal_rmse_hourly",
    "cal_rmse_daily",
    "cal_bias_hourly",
    "cal_bias_daily",
    "cal_r2_hourly",
    "cal_r2_daily",
    "eval_n_hourly",
    "eval_n_daily",
    "eval_rmse_hourly",
    "eval_rmse_daily",
    "eval_bias_hourly",
    "eval_bias_daily",
    "eval_r2_hourly",
    "eval_r2_daily",
    "runtime_elapsed_sec"
  ),
  category = c(
    "calibration_setting",
    "calibration_setting",
    "calibration_setting",
    "calibration_setting",
    "objective_definition",
    "objective_definition",
    "parameter_bound",
    "parameter_bound",
    "parameter_bound",
    "parameter_bound",
    "parameter_bound",
    "parameter_bound",
    "reference_parameter",
    "reference_parameter",
    "reference_parameter",
    "objective_scaling",
    "objective_scaling",
    "best_fit_result",
    "best_fit_result",
    "calibration_performance",
    "calibration_performance",
    "calibration_performance",
    "calibration_performance",
    "calibration_performance",
    "calibration_performance",
    "calibration_performance",
    "calibration_performance",
    "evaluation_performance",
    "evaluation_performance",
    "evaluation_performance",
    "evaluation_performance",
    "evaluation_performance",
    "evaluation_performance",
    "evaluation_performance",
    "evaluation_performance",
    "runtime"
  ),
  unit_or_scale = c(
    "years",
    "count",
    "TRUE/FALSE",
    "count",
    "relative weight",
    "relative weight",
    "fraction",
    "fraction",
    "fraction",
    "fraction",
    "fraction",
    "fraction",
    "fraction",
    "fraction",
    "fraction",
    "normalized RMSE scale",
    "normalized RMSE scale",
    "dimensionless objective value",
    "same unit as TWD",
    "count",
    "count",
    "normalized TWD",
    "normalized TWD",
    "normalized TWD",
    "normalized TWD",
    "R²",
    "R²",
    "count",
    "count",
    "normalized TWD",
    "normalized TWD",
    "normalized TWD",
    "normalized TWD",
    "R²",
    "R²",
    "seconds"
  ),
  description = c(
    "Years used for parameter calibration.",
    "Number of optimization starts used in the multi-start calibration.",
    "Indicates whether the optimization was run in parallel.",
    "Number of CPU cores requested for the optimization.",
    "Weight assigned to the hourly RMSE term in the objective function.",
    "Weight assigned to the daily RMSE term in the objective function.",
    "Lower bound used for parameter F_E during optimization.",
    "Upper bound used for parameter F_E during optimization.",
    "Lower bound used for parameter F_TWD during optimization.",
    "Upper bound used for parameter F_TWD during optimization.",
    "Lower bound used for parameter F_theta during optimization.",
    "Upper bound used for parameter F_theta during optimization.",
    "Reference value of F_E used for objective scaling.",
    "Reference value of F_TWD used for objective scaling.",
    "Reference value of F_theta used for objective scaling.",
    "Hourly RMSE obtained for the reference parameter set and used to scale the hourly objective term.",
    "Daily RMSE obtained for the reference parameter set and used to scale the daily objective term.",
    "Final objective value of the best-fit parameter set.",
    "99th percentile of maximum daily water deficit range used to normalize simulated TWD.",
    "Number of valid hourly observations used to evaluate calibration performance.",
    "Number of valid daily aggregated observations used to evaluate calibration performance.",
    "Hourly RMSE between observed and simulated normalized TWD during the calibration period.",
    "Daily RMSE between observed and simulated normalized TWD during the calibration period.",
    "Hourly mean bias between observed and simulated normalized TWD during the calibration period.",
    "Daily mean bias between observed and simulated normalized TWD during the calibration period.",
    "Hourly coefficient of determination between observed and simulated normalized TWD during the calibration period.",
    "Daily coefficient of determination between observed and simulated normalized TWD during the calibration period.",
    "Number of valid hourly observations used to evaluate transfer performance outside the calibration period.",
    "Number of valid daily aggregated observations used to evaluate transfer performance outside the calibration period.",
    "Hourly RMSE between observed and simulated normalized TWD during the evaluation period.",
    "Daily RMSE between observed and simulated normalized TWD during the evaluation period.",
    "Hourly mean bias between observed and simulated normalized TWD during the evaluation period.",
    "Daily mean bias between observed and simulated normalized TWD during the evaluation period.",
    "Hourly coefficient of determination between observed and simulated normalized TWD during the evaluation period.",
    "Daily coefficient of determination between observed and simulated normalized TWD during the evaluation period.",
    "Elapsed wall-clock time of the optimization run."
  ),
  stringsAsFactors = FALSE
)

# 8) Write outputs ----

# Write the best-fit parameter table as CSV so that it can be reused
# directly in other scripts.
write.csv(
  best_fit_parameters,
  file = file.path(out_dir, "twist_best_fit_parameters.csv"),
  row.names = FALSE
)

# Write the metadata workbook as an Excel file with:
# 1) one sheet containing the calibration metadata
# 2) one sheet explaining all metadata columns
metadata_xlsx_file <- file.path(out_dir, "twist_best_fit_metadata.xlsx")

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "metadata")
openxlsx::writeData(wb, sheet = "metadata", x = best_fit_metadata)

openxlsx::addWorksheet(wb, "column_description")
openxlsx::writeData(wb, sheet = "column_description", x = metadata_column_description)

openxlsx::setColWidths(wb, sheet = "metadata", cols = 1:ncol(best_fit_metadata), widths = "auto")
openxlsx::setColWidths(
  wb,
  sheet = "column_description",
  cols = 1:ncol(metadata_column_description),
  widths = "auto"
)

openxlsx::saveWorkbook(wb, file = metadata_xlsx_file, overwrite = TRUE)

cat("Optimization finished.\n")
cat("Outputs written to:\n")
cat(normalizePath(out_dir), "\n")
