# ================================================================
# 3_analyze_parameter_uncertainty.R
#
# Uncertainty workflow for the Tree Water Imbalance and Storage Tracker (TWIST) module
#
# This script:
# 1) Loads TWIST functions
# 2) Defines uncertainty analysis settings
# 3) Reads input, validation, and best-fit parameter data
# 4) Prepares model inputs and validation targets
# 5) Builds a scaled objective identical to the calibration script
# 6) Samples the bounded parameter space around the best-fit solution
# 7) Evaluates model performance for all sampled parameter sets
# 8) Defines acceptable parameter sets relative to the best-fit solution
# 9) Writes two output tables
#
# The script is intentionally minimal and only returns:
# - one table with all evaluated parameter sets and their performance
# - one table summarising the resulting parameter uncertainty ranges
#
# Note: All water-related variables (E, TWD, W) must share the same unit.
# ================================================================

# 1) Load functions ----
source("0_TWIST_functions.R")

# Load package used to write the summary workbook
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop(
    "Package 'openxlsx' is required to write the uncertainty summary workbook. ",
    "Please install it first with install.packages('openxlsx')."
  )
}

# 2) Define uncertainty analysis settings ----

# Input file paths
input_file <- "data/TWIST_input_data.rds"   
validation_file <- "data/TWIST_validation_data.rds"  

best_fit_file <- file.path(
  "output",
  "twist_best_fit_calibration",
  "twist_best_fit_parameters.csv"
)


# Calibration years
# These must match the years used in the best-fit calibration script.
calibration_years <- c(2018)

# TWIST parameter bounds
lower <- c(F_E = 0.0, F_TWD = 0.0, F_theta = 0.01)
upper <- c(F_E = 1.0, F_TWD = 1.0, F_theta = 1.0)

# Reference parameter set used to scale the hourly and daily RMSE terms
ref_par <- c(F_E = 0.6, F_TWD = 0.3, F_theta = 0.7)

# Sampling settings
# n_random_samples defines how many random parameter sets are drawn in
# addition to the best-fit and reference parameter sets.
n_random_samples <- 8000
random_seed <- 123

# Definition of acceptable parameter sets
# A sampled set is considered acceptable if its objective is within
# acceptable_multiplier x best-fit objective.
acceptable_multiplier <- 1.05

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
out_dir <- "output/twist_parameter_uncertainty"
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

# Lightweight TWIST runner used inside the sampling routine
# This reproduces the same model equations as TWIST_functions.R but avoids
# repeated data-frame overhead during uncertainty analysis.
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
# used to build the uncertainty objective.
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
      TWD_norm = NULL
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
    TWD_norm = twd_norm
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
      TWD_norm = raw_res$TWD_norm
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
    TWD_norm = raw_res$TWD_norm
  )
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

# Evaluate one sampled parameter set and return one compact result row.
evaluate_parameter_set <- function(par, sample_id, sample_type) {
  eval_res <- evaluate_parameters(par)

  if (is.null(eval_res$TWD_norm) || !is.finite(eval_res$objective)) {
    return(NULL)
  }

  metrics_cal <- compute_metrics_for_mask(
    obs = obs_full,
    sim = eval_res$TWD_norm,
    mask_obs = cal_mask_obs,
    day_factor = cal_day_factor
  )

  metrics_eval <- compute_metrics_for_mask(
    obs = obs_full,
    sim = eval_res$TWD_norm,
    mask_obs = eval_mask_obs,
    day_factor = eval_day_factor
  )

  data.frame(
    sample_id = sample_id,
    sample_type = sample_type,
    F_E = par["F_E"],
    F_TWD = par["F_TWD"],
    F_theta = par["F_theta"],
    objective = eval_res$objective,
    q99_mdwd = eval_res$q99_mdwd,
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
    stringsAsFactors = FALSE
  )
}

# Draw random parameter sets within the bounded TWIST parameter space.
sample_parameter_sets <- function(n, lower, upper) {
  mat <- cbind(
    runif(n, lower["F_E"], upper["F_E"]),
    runif(n, lower["F_TWD"], upper["F_TWD"]),
    runif(n, lower["F_theta"], upper["F_theta"])
  )
  colnames(mat) <- c("F_E", "F_TWD", "F_theta")
  mat
}

# 4) Read and prepare data ----

# Read module input and validation data
# The exact file names and paths can be adapted later.
df_input <- readRDS(input_file)
validation_data <-  readRDS(validation_file)

# Read the best-fit parameter set produced by the `2_calibrate_best_fit_parameters.R` script.
best_fit_tbl <- read.csv(best_fit_file)
colnames(best_fit_tbl) <- gsub("\\.", "_", colnames(best_fit_tbl))

required_best_fit_cols <- c("F_E", "F_TWD", "F_theta")
missing_best_fit_cols <- setdiff(required_best_fit_cols, names(best_fit_tbl))
if (length(missing_best_fit_cols) > 0) {
  stop(
    "The best-fit parameter file is missing required columns: ",
    paste(missing_best_fit_cols, collapse = ", ")
  )
}

best_fit_par <- as.numeric(best_fit_tbl[1, required_best_fit_cols])
names(best_fit_par) <- required_best_fit_cols

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
# These vectors are used repeatedly during uncertainty analysis.
time_vec <- df_input[[col_names$col_time]]
date_vec <- as.Date(time_vec)
hour_vec <- as.POSIXlt(time_vec)$hour
year_vec <- as.POSIXlt(time_vec)$year + 1900

# Extract model inputs as plain vectors for faster repeated evaluation
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

# 6) Sample the parameter space ----

# Draw random samples from the bounded parameter space and prepend the
# best-fit and reference parameter sets for direct comparison.
set.seed(random_seed)
random_samples <- sample_parameter_sets(n_random_samples, lower, upper)

sample_matrix <- rbind(best_fit_par, ref_par, random_samples)
colnames(sample_matrix) <- c("F_E", "F_TWD", "F_theta")

sample_type <- c(
  "best_fit",
  "reference",
  rep("random", nrow(random_samples))
)

# 7) Run uncertainty analysis ----

cat("\nStarting parameter uncertainty analysis ...\n")
cat("Random samples:", n_random_samples, "\n")
cat("Total evaluations:", nrow(sample_matrix), "\n")
cat("Parallel:", use_parallel, "\n")
cat("Requested cores:", n_cores, "\n")

analysis_time <- system.time({
  if (use_parallel && n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c(
        "E_vec", "theta_vec", "W_vec",
        "date_vec", "hour_vec",
        "obs_full", "obs_hourly_cal", "obs_daily_cal",
        "cal_mask_obs", "eval_mask_obs",
        "cal_day_factor", "eval_day_factor",
        "scale_hourly", "scale_daily",
        "w_hourly", "w_daily",
        "compute_q99_mdwd_fast",
        "run_twist_fast",
        "rmse_fast", "bias_fast", "r2_fast", "daily_mean_fast",
        "evaluate_parameters_raw",
        "evaluate_parameters",
        "compute_metrics_for_mask",
        "evaluate_parameter_set",
        "sample_matrix",
        "sample_type"
      ),
      envir = environment()
    )

    results_list <- parallel::parLapply(
      cl,
      X = seq_len(nrow(sample_matrix)),
      fun = function(i) {
        evaluate_parameter_set(
          par = sample_matrix[i, ],
          sample_id = i,
          sample_type = sample_type[i]
        )
      }
    )
  } else {
    results_list <- lapply(
      seq_len(nrow(sample_matrix)),
      function(i) {
        evaluate_parameter_set(
          par = sample_matrix[i, ],
          sample_id = i,
          sample_type = sample_type[i]
        )
      }
    )
  }
})

results_list <- results_list[!vapply(results_list, is.null, logical(1))]
if (length(results_list) == 0L) {
  stop("Uncertainty analysis did not return any valid result.")
}

uncertainty_results <- do.call(rbind, results_list)
uncertainty_results <- uncertainty_results[order(uncertainty_results$objective), ]
row.names(uncertainty_results) <- NULL

# 8) Define acceptable parameter sets ----

# Acceptable solutions are defined relative to the objective of the
# calibrated best-fit parameter set.
best_fit_result <- uncertainty_results[
  uncertainty_results$sample_type == "best_fit",
  ,
  drop = FALSE
]

if (nrow(best_fit_result) != 1L) {
  stop("The best-fit parameter set could not be identified in the uncertainty results.")
}

best_fit_objective <- best_fit_result$objective[1]
acceptable_threshold <- acceptable_multiplier * best_fit_objective

acceptable_results <- uncertainty_results[
  uncertainty_results$objective <= acceptable_threshold,
  ,
  drop = FALSE
]

if (nrow(acceptable_results) == 0L) {
  stop("No acceptable parameter sets were identified.")
}

# Create a compact summary table of uncertainty ranges and key settings.
uncertainty_summary <- data.frame(
  calibration_years = paste(calibration_years, collapse = ","),
  n_random_samples = n_random_samples,
  n_total_evaluated = nrow(uncertainty_results),
  n_acceptable = nrow(acceptable_results),
  acceptable_multiplier = acceptable_multiplier,
  best_fit_objective = best_fit_objective,
  acceptable_threshold = acceptable_threshold,
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
  acceptable_F_E_min = min(acceptable_results$F_E, na.rm = TRUE),
  acceptable_F_E_median = median(acceptable_results$F_E, na.rm = TRUE),
  acceptable_F_E_max = max(acceptable_results$F_E, na.rm = TRUE),
  acceptable_F_TWD_min = min(acceptable_results$F_TWD, na.rm = TRUE),
  acceptable_F_TWD_median = median(acceptable_results$F_TWD, na.rm = TRUE),
  acceptable_F_TWD_max = max(acceptable_results$F_TWD, na.rm = TRUE),
  acceptable_F_theta_min = min(acceptable_results$F_theta, na.rm = TRUE),
  acceptable_F_theta_median = median(acceptable_results$F_theta, na.rm = TRUE),
  acceptable_F_theta_max = max(acceptable_results$F_theta, na.rm = TRUE),
  runtime_elapsed_sec = unname(analysis_time["elapsed"]),
  stringsAsFactors = FALSE
)

# Create a second table that explains all summary columns.
# This sheet is intended for direct inspection in Excel.
uncertainty_summary_column_description <- data.frame(
  column_name = c(
    "calibration_years",
    "n_random_samples",
    "n_total_evaluated",
    "n_acceptable",
    "acceptable_multiplier",
    "best_fit_objective",
    "acceptable_threshold",
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
    "acceptable_F_E_min",
    "acceptable_F_E_median",
    "acceptable_F_E_max",
    "acceptable_F_TWD_min",
    "acceptable_F_TWD_median",
    "acceptable_F_TWD_max",
    "acceptable_F_theta_min",
    "acceptable_F_theta_median",
    "acceptable_F_theta_max",
    "runtime_elapsed_sec"
  ),
  category = c(
    "analysis_setting",
    "analysis_setting",
    "analysis_result",
    "analysis_result",
    "acceptance_definition",
    "acceptance_definition",
    "acceptance_definition",
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
    "accepted_range",
    "accepted_range",
    "accepted_range",
    "accepted_range",
    "accepted_range",
    "accepted_range",
    "accepted_range",
    "accepted_range",
    "accepted_range",
    "runtime"
  ),
  unit_or_scale = c(
    "years",
    "count",
    "count",
    "count",
    "relative multiplier",
    "dimensionless objective value",
    "dimensionless objective value",
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
    "fraction",
    "fraction",
    "fraction",
    "fraction",
    "fraction",
    "fraction",
    "fraction",
    "fraction",
    "fraction",
    "seconds"
  ),
  description = c(
    "Years used for parameter calibration and therefore used to define the best-fit solution.",
    "Number of random parameter sets drawn in addition to the best-fit and reference parameter sets.",
    "Total number of parameter sets successfully evaluated during the uncertainty analysis.",
    "Number of parameter sets classified as acceptable based on the objective threshold.",
    "Multiplier applied to the best-fit objective to define the acceptance threshold.",
    "Objective value of the calibrated best-fit parameter set.",
    "Maximum objective value allowed for a parameter set to be classified as acceptable.",
    "Weight assigned to the hourly RMSE term in the objective function.",
    "Weight assigned to the daily RMSE term in the objective function.",
    "Lower bound used for parameter F_E during sampling.",
    "Upper bound used for parameter F_E during sampling.",
    "Lower bound used for parameter F_TWD during sampling.",
    "Upper bound used for parameter F_TWD during sampling.",
    "Lower bound used for parameter F_theta during sampling.",
    "Upper bound used for parameter F_theta during sampling.",
    "Reference value of F_E used for objective scaling.",
    "Reference value of F_TWD used for objective scaling.",
    "Reference value of F_theta used for objective scaling.",
    "Hourly RMSE obtained for the reference parameter set and used to scale the hourly objective term.",
    "Daily RMSE obtained for the reference parameter set and used to scale the daily objective term.",
    "Minimum accepted value of parameter F_E.",
    "Median accepted value of parameter F_E.",
    "Maximum accepted value of parameter F_E.",
    "Minimum accepted value of parameter F_TWD.",
    "Median accepted value of parameter F_TWD.",
    "Maximum accepted value of parameter F_TWD.",
    "Minimum accepted value of parameter F_theta.",
    "Median accepted value of parameter F_theta.",
    "Maximum accepted value of parameter F_theta.",
    "Elapsed wall-clock time of the uncertainty analysis."
  ),
  stringsAsFactors = FALSE
)

# 9) Write outputs ----

# Write the full uncertainty results table as CSV so that it can be
# reused directly in downstream analyses.
write.csv(
  uncertainty_results,
  file = file.path(out_dir, "twist_uncertainty_results.csv"),
  row.names = FALSE
)

# Write the compact uncertainty summary as an Excel file with:
# 1) one sheet containing the summary values
# 2) one sheet explaining all summary columns
summary_xlsx_file <- file.path(out_dir, "twist_uncertainty_summary.xlsx")

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "summary")
openxlsx::writeData(wb, sheet = "summary", x = uncertainty_summary)

openxlsx::addWorksheet(wb, "column_description")
openxlsx::writeData(
  wb,
  sheet = "column_description",
  x = uncertainty_summary_column_description
)

openxlsx::setColWidths(wb, sheet = "summary", cols = 1:ncol(uncertainty_summary), widths = "auto")
openxlsx::setColWidths(
  wb,
  sheet = "column_description",
  cols = 1:ncol(uncertainty_summary_column_description),
  widths = "auto"
)

openxlsx::saveWorkbook(wb, file = summary_xlsx_file, overwrite = TRUE)

cat("Uncertainty analysis finished.\n")
cat("Outputs written to:\n")
cat(normalizePath(out_dir), "\n")