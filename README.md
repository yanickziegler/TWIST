# **README**


This repository provides an R implementation of the Tree Water Imbalance and Storage Tracker (TWIST) module presented in the manuscript:

**Ziegler. et al. (2025, in prep.). *A simple framework for linking tree water deficit dynamics to drought risk across scales*.**

The material (code, example data, workflow) enables full reproducibility of the module and offers a minimal, easy-to-use template for applying it to other datasets and ecosystems.

---

## **1. Contents of the repository**

| File                       | Description                                                                                           |
| -------------------------- | ----------------------------------------------------------------------------------------------------- |
| **TWIST_functions.R**      | Core functions implementing the TWIST module, following Eqs. 1–5 in the manuscript.                   |
| **example_run.R**          | Example script that loads the input data, executes the module, and visualises TWD and RWC time series.|
| **example_input_data.rds** | Example dataset (Štítná site, 2018–2020) containing the minimal inputs required for the module.       |
| **README.md**              | Documentation of input structure, workflow, and interpretation of module output.                      |

---

## **2. Overview of the TWIST module**

TWIST simulates:

1. **Tree Water Deficit (TWD)** – the cumulative internal water shortage caused by the imbalance between transpiration and water uptake.
2. **Relative Water Content (RWC)** – the fraction of the plant’s water pool currently filled.

The module uses only:

**For TWD calculation:**
* Transpiration (`E`)
* Relative soil water content (`θ_rel`)
* Three parameters (`F_E`, `F_TWD`, `F_θ`)

**For RWC calculation:**
* TWD
* Tree water pool size (W), here estimated from dry wood biomass (`m_dry_wood_kg.m2`) and the ratio of saturated to oven-dry wood density (`rho_sat`, `rho_dry`)

and can therefore be integrated into a wide range of vegetation or land-surface models.

All equations implemented here correspond exactly to those described in the manuscript.

---

## **3. Input data requirements**

The example dataset (`example_input_data.rds`) contains:

| Column               | Meaning                                                       | Units                                              |
| ------------------   | ------------------------------------------------------------- | -------------------------------------------------- |
| `datetime`           | Timestamp                                                     | POSIXct, hourly timestep                           |
| `transpiration_l.m2` | Transpirational water loss                                    | liters water per hour and m² ground                |
| `theta_rel`          | Relative soil water content (0 = wilting; 1 = field capacity) | –                                                  |
| `m_dry_wood_kg.m2`   | Oven-dry tree biomass contributing to tree water storage      | kg per m² ground                                   |

### **Important:**
The module can be used with flexible units, and on a tree or on a stand-level.
However, all water-related variables (E, TWD, W) must concistently share the **same unit**.

---

## **4. Model parameters**

### **TWD parameters**

| Parameter | Meaning                                                          |
| --------- | ---------------------------------------------------------------- |
| `F_E`     | Fraction of transpiration directly supplied by soil water uptake |
| `F_TWD`   | Fraction of the current TWD that can be refilled per timestep    |
| `F_theta` | Threshold scaling soil water limitation                          |

### **Water pool parameters**

| Parameter | Meaning                          |
| --------- | -------------------------------- |
| `rho_sat` | Saturated wood density (kg dm⁻³) |
| `rho_dry` | Oven-dry wood density (kg dm⁻³)  |

---
## **5. Running the TWIST module**

An example of how to run TWIST is provided in `Example_run.R`. 

It follows the following logic:
1) Load TWIST functions
2) Read input data
3) Define module parameters 
4) Define column names for input data
5) Calculate the water pool size
6) Run module
7) Visualize results

For different applications, input data, parameter values, and column names can be adjusted.

---


## **6. Interpretation of output variables**

* TWD increases when transpiration exceeds uptake and declines when uptake exceeds transpiration.
* A TWD of 0 corresponds to a RWC of 1 and represents full hydration relative to the assumed water-pool size.
* Higher TWD and lower RWC values indicate more severe internal water depletion.

---

## **7. Reproducibility**

Running the two provided scripts (`TWIST_functions.R` and `Example_run.R`) will reproduce the TWD and RWC dynamics for the Štítná forest 2018-2020 shown in the manuscript.

The module is agnostic to the source of the input data and can be coupled to ecosystem, vegetation, and land-surface models as long as transpiration, relative soil water content, and wood biomass (or a different estimate for tree-water-pool size) are available or can be estimated.

---

## **8. License & citation**
Unless specified otherwise, the code is distributed under an MIT License to encourage reuse and adaptation.

If using or modifying this code in your own work, please cite:

INCLUDE REFERENCE ONCE THE PAPER IS PUBLISHED

---

## **9. Contact**

For questions, extensions, or reporting issues:

**Yanick Ziegler**, yanick.ziegler@kit.edu

Karlsruhe Institute of Technology, Institute of Meteorology and Climate Research - Atmospheric Environmental Research (KIT/IMKIFU)

