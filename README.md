# TWIST

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17979713.svg)](https://doi.org/10.5281/zenodo.17979713)

This repository provides an R implementation of the **Tree Water Imbalance and Storage Tracker (TWIST)** module presented in the manuscript:

**Ziegler et al. (2026, in prep.). _A simple framework for linking tree water deficit dynamics to drought risk across scales_.**

The repository contains the core TWIST functions, a reproducible example workflow, and two additional scripts for parameter calibration and parameter uncertainty analysis. Together, these files allow reproduction of the **TWD** and **RWC** time series shown in the accompanying publication, as well as derivation of the corresponding **best-fit parameter set** and the set of **acceptable parameter solutions**. At the same time, the code is intended as a compact template for applying TWIST to new datasets or coupling it to ecosystem models that provide the required inputs. 

---

## 1. Repository contents

| File / folder | Description |
| --- | --- |
| `0_TWIST_functions.R` | Core functions implementing the TWIST module, following Eqs. 1–5 in the manuscript. |
| `1_run_TWIST_example.R` | Minimal example workflow showing how to run TWIST with example input data and a best-fit parameter set. Reproduces the TWD and RWC time series presented in the publication. |
| `2_calibrate_best_fit_parameters.R` | Calibration workflow that retrieves a best-fit TWIST parameter set from input data and validation data. |
| `3_analyze_parameter_uncertainty.R` | Uncertainty workflow that evaluates parameter uncertainty around the calibrated best-fit solution and identifies acceptable parameter sets. |
| `data/` | Input data used by the example and calibration workflows. |
| `output/` | Output directory created by the calibration and uncertainty scripts. |
| `TWIST_module.Rproj` | R project file for opening the repository with the correct working directory and relative path structure in RStudio. |
| `README.md` | Documentation of repository structure, required inputs, workflow logic, and output files. |

---

## 2. Working directory and relative paths

All scripts use relative file paths. To ensure that these paths work correctly, open the repository via `TWIST_module.Rproj` and run the scripts from the repository root directory.

## 3. Software requirements

The repository was developed and tested in **R (version 4.4 or later)**. The scripts use only a small number of additional packages, but these must be installed before running the workflows.

Required packages:
- `ggplot2`
- `openxlsx`

The calibration and uncertainty scripts can optionally use parallel processing via base R functionality, so runtime will depend on the available hardware and the number of CPU cores used.

Packages can be installed with:

```r
install.packages(c("ggplot2", "openxlsx"))
```

## 4. Overview of the TWIST module

TWIST simulates two state variables:

1. **Tree Water Deficit (TWD)** – the cumulative internal water shortage caused by the imbalance between transpiration and water uptake.
2. **Relative Water Content (RWC)** – the fraction of the plant water pool currently filled.

The module uses only a small set of inputs and parameters.

**For TWD calculation:**
- Transpiration (`E`)
- Relative soil water content (`theta_rel`)
- Three parameters (`F_E`, `F_TWD`, `F_theta`)

**For RWC calculation:**
- TWD
- Tree water pool size (`W`), here estimated from dry wood biomass (`m_dry_wood_kg.m2`) and the ratio of saturated to oven-dry wood density (`rho_sat`, `rho_dry`)

All equations implemented in `0_TWIST_functions.R` correspond to those described in the manuscript. The implementation is intentionally compact so that TWIST can be used as a stand-alone module or embedded into larger model frameworks.

---

## 5. What can be reproduced with this repository

This repository supports three complementary use cases:

### Reproduction of the published example time series
Running `1_run_TWIST_example.R` reproduces the **TWD** and **RWC** time series shown for the example site in the accompanying publication, using the best-fit TWIST parameter set reported for that setup. 

### Reproduction of the calibrated parameter solution
Running `2_calibrate_best_fit_parameters.R` retrieves a best-fit parameter set from the provided input and validation data and writes:
- `twist_best_fit_parameters.csv`
- `twist_best_fit_metadata.xlsx` 

### Reproduction of the uncertainty analysis
Running `3_analyze_parameter_uncertainty.R` evaluates parameter uncertainty around the calibrated best-fit solution and writes:
- `twist_uncertainty_results.csv`
- `twist_uncertainty_summary.xlsx` 

These workflows are designed so that users can either reproduce the published TWIST application directly or adapt the scripts to their own data and modelling context. 

---

## 6. Input data requirements

### Example workflow input
The example run script expects a TWIST input dataset in `data/` containing the minimal variables required to run the module:

| Column | Meaning | Units / format |
| --- | --- | --- |
| `datetime` | Timestamp | `POSIXct`, hourly timestep |
| `transpiration_l.m2` | Transpirational water loss | litres water per hour and m² ground (other units are possible, but consistency with the water pool W unit is necessary) |
| `theta_rel` | Relative soil water content | dimensionless, 0 = wilting point, 1 = field capacity |
| `m_dry_wood_kg.m2` | Oven-dry tree biomass contributing to tree water storage | kg per m² ground |

### Calibration and uncertainty workflow inputs
The calibration and uncertainty scripts use:
- a TWIST input dataset containing the model drivers and wood biomass needed to run TWIST, and
- a validation dataset containing the observed TWD time series used for model-data comparison and parameter evaluation. 

### Validation data requirements and timestamp matching

For the calibration and uncertainty workflows, the validation dataset must contain the observed normalized TWD time series used for model-data comparison. In the current implementation, this dataset is expected to include a column named `mean_TWD_norm`.

In addition, the timestamps of the validation data must match the simulation timestamps exactly wherever observations are available. The calibration and uncertainty scripts align simulated and observed values by timestamp, so mismatched or inconsistent time formats will prevent correct comparison.

Before comparison, the simulated TWD is normalized over the full available time series using the same normalization approach described in the accompanying manuscript. This normalization is part of the calibration and uncertainty workflow and should be preserved when adapting the scripts to new datasets.


### Important note on units
The module can be used with flexible units and on different aggregation levels, for example per tree or per unit ground area. However, all water-related variables used within the TWIST calculations must consistently share the same unit, especially `W` and `E`, from which `TWD` inherits its unit. This is essential for both the example workflow and any adaptation to new datasets or coupled models. 

---

## 7. Model parameters

### TWD parameters

| Parameter | Meaning |
| --- | --- |
| `F_E` | Fraction of transpiration directly supplied by soil water uptake |
| `F_TWD` | Fraction of the current TWD that can be refilled per timestep |
| `F_theta` | Threshold for the soil water limitation onset |

### Water pool parameters

| Parameter | Meaning |
| --- | --- |
| `rho_sat` | Saturated wood density |
| `rho_dry` | Oven-dry wood density |

In `1_run_TWIST_example.R`, the TWD parameters are set to the best-fit solution used to reproduce the published example time series. In `2_calibrate_best_fit_parameters.R`, these parameters are estimated from the provided data. In `3_analyze_parameter_uncertainty.R`, uncertainty is assessed around the calibrated best-fit solution within predefined parameter bounds. 

---

## 8. Recommended workflow

The repository is structured so that users can enter at different levels, depending on their goal.

### 8.1 Run the published example
Use `1_run_TWIST_example.R` to see the minimal TWIST workflow:
1. Load the TWIST core functions
2. Read example input data
3. Define model parameters
4. Define input column names
5. Calculate tree water pool size
6. Run TWIST
7. Visualize TWD and RWC

This script is intended both as a direct reproduction of the published example time series and as a template for running TWIST on new data. 

### 8.2 Retrieve a best-fit parameter set
Use `2_calibrate_best_fit_parameters.R` to derive a best-fit TWIST parameter set from the provided input and validation data. Before model-data comparison, the simulated TWD is normalized over the full available time series using the same normalization approach described in the accompanying manuscript. Observed normalized TWD is expected in the validation file. The calibration objective is based on normalized TWD and combines hourly and daily RMSE terms. The script performs a multi-start calibration and saves a compact parameter table together with a metadata workbook documenting the calibration setup and model performance. Parameter estimation is performed for the defined calibration period, while model performance is reported separately for the calibration and evaluation periods.

Runtime note:
- The script uses a multi-start optimization approach and can use parallel processing.
- Runtime depends on the number of start values, the selected number of CPU cores, and the available hardware.


### 8.3 Analyze parameter uncertainty
Use `3_analyze_parameter_uncertainty.R` after calibration to assess uncertainty around the best-fit solution. The script samples parameter sets within the bounded TWIST parameter space, evaluates them using the same normalized TWD framework and objective function as in the calibration, and identifies acceptable solutions relative to the best-fit objective.

Prerequisite:
- `output/twist_best_fit_calibration/twist_best_fit_parameters.csv` must already exist, because it is read by the uncertainty script.

Runtime note:
- The script evaluates many parameter sets and can use parallel processing.
- Runtime depends strongly on the number of sampled parameter sets, the selected number of CPU cores, and the available hardware.

---

## 9. Output files

### 9.1 Calibration outputs
`2_calibrate_best_fit_parameters.R` writes its outputs to:

`output/twist_best_fit_calibration/`

Files:
- `twist_best_fit_parameters.csv`  
  Compact table containing the best-fit values of `F_E`, `F_TWD`, and `F_theta`.
- `twist_best_fit_metadata.xlsx`  
  Metadata workbook describing the calibration setup and model performance. The workbook contains one sheet with the metadata values and one sheet explaining all columns. 

### 9.2 Uncertainty outputs
`3_analyze_parameter_uncertainty.R` writes its outputs to:

`output/twist_parameter_uncertainty/`

Files:
- `twist_uncertainty_results.csv`  
  Full results table containing all evaluated parameter sets and their associated performance metrics.
- `twist_uncertainty_summary.xlsx`  
  Compact summary workbook describing the uncertainty analysis settings and the resulting ranges of acceptable parameter values. The workbook contains one sheet with the summary values and one sheet explaining all columns.
---

## 10. Interpretation of output variables

### TWIST model output
- **TWD** increases when transpiration exceeds uptake and declines when uptake exceeds transpiration.
- A **TWD of 0** corresponds to **RWC = 1** and represents full hydration relative to the assumed water-pool size.
- Higher **TWD** and lower **RWC** values indicate more severe internal water depletion. 

### Calibration output
The best-fit parameter set represents the parameter combination that minimizes the calibration objective function defined in `2_calibrate_best_fit_parameters.R`. The metadata workbook documents the optimization settings, objective scaling, parameter bounds, and calibration and evaluation performance metrics. 

### Uncertainty output
The uncertainty analysis evaluates parameter sets relative to the calibrated best-fit solution. The summary workbook reports the number of acceptable solutions, the acceptance threshold, and the resulting parameter ranges for `F_E`, `F_TWD`, and `F_theta`. The full results table can be used for further exploration, visualization, or sensitivity analysis.

---

## 11. Using TWIST with new data or in coupled models

The repository is not limited to the included example data. TWIST can also be coupled to ecosystem models, provided that the required inputs can be supplied:
- transpiration
- relative soil water content
- an estimate of wood biomass or tree water pool size

For such applications, users can:
- adapt the input file paths and column names in the scripts,
- replace the included data with model-generated inputs,
- recalibrate `F_E`, `F_TWD`, and `F_theta` for their own system, and
- use the example script as a minimal template for implementation. 

---

## 12. Reproducibility

Together, the scripts in this repository provide a compact and reproducible TWIST workflow:
- `0_TWIST_functions.R` contains the core model implementation
- `1_run_TWIST_example.R` reproduces the published example TWD and RWC time series
- `2_calibrate_best_fit_parameters.R` reproduces the best-fit parameter set
- `3_analyze_parameter_uncertainty.R` reproduces the uncertainty analysis and acceptable parameter solutions 

This combination is intended to support both transparent reproduction of the manuscript results and straightforward reuse of the TWIST module in other studies and modelling frameworks. 

---

## 13. License and citation

The code is distributed under an MIT License to encourage reuse and adaptation.

If using or modifying this code in your own work, please cite the accompanying paper:

**TO DO: include full manuscript reference once published**

For the code itself:

**Ziegler, Y., Ruehr, N. and Grote, R. (2026). _Tree Water Imbalance and Storage Tracker (TWIST) module_. Zenodo. https://doi.org/10.5281/zenodo.17979712**

---

## 14. Contact

For questions, extensions, or reporting issues:

**Yanick Ziegler**  
yanick.ziegler@kit.edu

Karlsruhe Institute of Technology, Institute of Meteorology and Climate Research – Atmospheric Environmental Research (KIT/IMK-IFU)
