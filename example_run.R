# ================================================================
# example_run.R
# 
# Example workflow to run the Tree Water Imbalance and Storage Tracker (TWIST) module
# presented in Ziegler et al. (2025, in prep.).
# 
# This script:
# 1) Loads TWIST functions
# 2) Reads input data
# 3) Defines module parameters 
# 4) Defines column names for input data
# 5) Calculates the water pool size
# 6) Runs the module
# 7) Visualizes results
# 
# For different applications, input data, parameter values, and column names
# can be adjusted.
#
# Note: All water-related variables (E, TWD, W) must share the same unit.
# In this example, transpiration is given in litres per m² ground and hour.
# ================================================================

# 1) Load functions ----
source("TWIST_functions.R")
library(dplyr)
library(ggplot2)

# 2) Read input data ---- 
# Contains: datetime, transpiration_l.m2, theta_rel, m_dry_wood_kg.m2
df_input <- readRDS("example_input_data.rds")

# 3) Define module parameters ----
# The parameters are defined for the Štítná Fagus sylvatica setup presented
# in the manuscript, but can be adjusted for other use cases.

# Define TWD parameters 
params_TWD <- list(
  F_E     = 0.6,   # Fraction of transpiration directly supplied by uptake
  F_TWD   = 0.3,   # Fraction of current TWD that can be refilled per timestep
  F_theta = 0.7    # Soil moisture threshold scaling uptake downregulation
)

# Define water pool parameters
params_water_pool <- list(
  rho_sat = 1.07,  # Fully saturated wood density [kg dm³]
  rho_dry = 0.58   # Oven-dried wood density [kg dm³]
)

# 4) Define column names for input data ----
# Can be adjusted for different input formats.
col_names <- list(
  col_time      = "datetime",              # POSIXct
  col_E         = "transpiration_l.m2",    # Transpirational water loss. Here: litres water per hour and m² ground, but other units are possible
  col_theta     = "theta_rel",             # Relative soil water content (0 = wilting; 1 = field capacity)
  col_m_wood    = "m_dry_wood_kg.m2",      # Oven-dry tree biomass contributing to tree water storage (kg per m² ground)
  col_W         = "W"                      # Water pool size (calculated in step 5); must be consistent with the unit of E; here: litres water per m² ground
)

# 5) Calculate the water pool size ----
# Estimate water pool size W for each timestep based on wood biomass and densities
df_input[[col_names$col_W]] <- compute_water_pool(
  m_wood_dry       = df_input[[col_names$col_m_wood]],
  params_water_pool = params_water_pool
)

# 6) Run model ----
output <- run_TWD_timeseries_df(
  df         = df_input,
  params_TWD = params_TWD,
  TWD_old    = 0,         # Initial TWD at the beginning of the simulation
  col_names  = col_names
)

# 7) Visualise output ----

# Tree Water Deficit
p_TWD <- ggplot(output, aes(datetime, TWD)) + 
  geom_line() + 
  labs(x = "Time", y = "Tree Water Deficit (same unit as E)") +
  theme_bw()

p_TWD

# Relative Water Content
p_RWC <- ggplot(output, aes(datetime, RWC)) + 
  geom_line() + 
  labs(x = "Time", y = "Relative water content (-)") +
  theme_bw()

p_RWC
