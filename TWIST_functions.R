# ======================================================================
# TWIST_functions.R
#
# Core implementation of the Tree Water Imbalance and Storage Tracker (TWIST) module presented in:
#
#   Ziegler et al. (2025, in prep.)
#   "A simple framework for linking tree water deficit dynamics to drought risk across scales"
#
# This script implements Eqs. (1)–(5) of the manuscript:
#   Eq. (1) – Update of Tree Water Deficit (TWD)
#   Eq. (2) – Water uptake as a function of transpiration and TWD
#   Eq. (3) – Soil water limitation function
#   Eq. (4) – Relative Water Content (RWC)
#   Eq. (5) – Tree water pool estimation based on wood biomass and density
#
# An exemplary application reproducing manuscript results is provided in example_run.R
# ======================================================================

# ========== 1. Soil limitation function (Eq. 3) ==========
f_soil <- function(theta_rel, F_theta) {
  return(pmin(theta_rel / F_theta, 1))
}

# ========== 2. Water uptake (Eq. 2) ==========
compute_uptake <- function(E, TWD_old, theta_rel, params_TWD) {
  f_soil_val <- f_soil(theta_rel, params_TWD$F_theta)
  uptake <- (params_TWD$F_E * E + params_TWD$F_TWD * TWD_old) * f_soil_val
  return(uptake)
}

# ========== 3. Update Tree Water Deficit (Eq. 1) ==========
update_TWD <- function(E, TWD_old, theta_rel, params_TWD) {
  U <- compute_uptake(E, TWD_old, theta_rel, params_TWD)
  return(TWD_old + E - U)
}

# ========== 4. Tree water pool estimation (Eq. 5) ==========
compute_water_pool <- function(m_wood_dry, params_water_pool) {
  W <- (params_water_pool$rho_sat / params_water_pool$rho_dry - 1) * m_wood_dry         
  return(W)
}

# ========== 5. Relative Water Content (RWC) (Eq. 4) ==========
compute_rwc <- function(W, TWD) {
  return(pmax(0, (W - TWD) / W))  # bounded [0–1]
}

# ========== 6. TWD + RWC per timestep ==========

# E         ... Transpirational water loss. Unit Flexible: volume/mass of H2O, per tree or m² ground
# theta_rel ... Relative soil water content: 1 at field capacity, 0 at wilting point
# W         ... Tree water pool size [same unit as E]
# TWD_old   ... TWD at previous timestep [same unit as E]

run_TWD_model <- function(E, theta_rel, W, params_TWD, TWD_old = 0) {
  
  TWD_new <- update_TWD(E, TWD_old, theta_rel, params_TWD)
  
  RWC <- compute_rwc(W, TWD_new)
  
  return(list(TWD = TWD_new, RWC = RWC))
}

# ========== 7. Run model for timeseries input ==========

run_TWD_timeseries_df <- function(df, 
                                  params_TWD, 
                                  TWD_old = 0,               #Default: TWD_old 0 at first call
                                  col_names = list(
                                    col_time  = "datetime",
                                    col_E     = "transpiration",
                                    col_theta = "theta_rel",
                                    col_W     = "W"
                                  )) {      

  required_cols <- c(col_names$col_time,
                     col_names$col_E,
                     col_names$col_theta,
                     col_names$col_W)
  
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop("The following required columns are missing in df: ",
         paste(missing, collapse = ", "))
  }
  
  df_datetime <- df[[col_names$col_time]]
  df_E         <- df[[col_names$col_E]]
  df_theta_rel <- df[[col_names$col_theta]]
  df_W         <- df[[col_names$col_W]]
  
  n <- length(df_datetime)
  TWD_vec <- numeric(n)
  RWC_vec <- numeric(n)

  for (i in seq_len(n)) {
    result <- run_TWD_model(
      E             = df_E[i],
      theta_rel     = df_theta_rel[i],
      TWD_old       = TWD_old,
      W             = df_W[i],
      params_TWD    = params_TWD
      )
    
    TWD_vec[i] <- result$TWD
    RWC_vec[i] <- result$RWC

    TWD_old <- result$TWD
  }
  
  return(data.frame(
    datetime  = df_datetime,
    TWD       = TWD_vec,
    RWC       = RWC_vec
  ))
}


