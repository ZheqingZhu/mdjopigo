# MD-JoPiGo: Multi-Dimensional Joint Patient Data Generator

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Statistical Computing](https://img.shields.io/badge/R-v4.2.2+-blue.svg)](https://www.r-project.org/)

**MD-JoPiGo** (Multidimensional Joint Patient Individual-data Generator and Optimizer) is an open-source framework designed to bridge the **"dimensionality gap"** in clinical trial reporting. It allows researchers to reconstruct high-dimensional individual patient data (IPD) from the one-dimensional (1D) survival summaries (Kaplan-Meier curves) typically found in published literature.

---

## 🌟 Overview

Randomized controlled trials (RCTs) often report efficacy through isolated 1D marginal summaries (e.g., KM curves for the overall population and specific subgroups). This prevents the observation of intersectional effects (e.g., the survival of a patient who is both "Elderly" and "ECOG 2-3"). 

**MD-JoPiGo** solves this by:
1.  **Estimating Joint Distributions**: Using a **Maximum Entropy (MaxEnt)** approach to find the most unbiased joint frequencies across clinical strata.
2.  **Individual-Level Synthesis**: Using **Simulated Annealing (SA)** to assign clinical labels to individuals, ensuring the resulting multidimensional profiles honor the observed 1D survival dynamics.
3.  **Topological Calibration**: Allowing for "Structural Priors" to correct biases caused by covariate collinearity (e.g., age-related frailty).

---

## 📊 Data Requirements

To use MD-JoPiGo, you need to provide two CSV files reconstructed from KM curves (using tools like KM-PoPiGo or Guyot's method).

### 1. `Overall.csv` (The ITT Population)
This file represents the total population of the treatment arm. It provides the "pool" of survival times and event statuses.

| Column | Type | Description |
| :--- | :--- | :--- |
| `time` | Numeric | The time-to-event or censoring time. |
| `status` | Binary | Event indicator: `1` for event (death/progression), `0` for censored. |
| `group` | String | The label for the overall population (e.g., `Overall` or `ITT`). |

### 2. `Layer.csv` (The Stratified Subgroups)
This file contains the IPD for all reported subgroups. The algorithm uses these curves as constraints to "anchor" the multidimensional optimization.

| Column | Type | Description |
| :--- | :--- | :--- |
| `time` | Numeric | Survival time. |
| `status` | Binary | Event indicator (`1` or `0`). |
| `group` | String | **Critical:** The label for the subgroup (e.g., `sex1`, `age_ge65`, `ECOG0_1`). These labels must match the names defined in the R script's feature dictionary. |

---

## 🚀 Quick Start

### 1. Configure Feature Mapping
Open the main R script and define your clinical dimensions in the `my_features` list. This tells the algorithm which labels in `Layer.csv` belong to the same feature (e.g., Sex, Age, ECOG).

```R
my_features <- list(
  Sex  = c("sex1", "sex2"),
  Age  = c("age_ge65", "age_ne65"),
  ECOG = c("ECGO0_1", "ECGO2_3")
)
```
## 2. Structural Calibration (Optional)
If your clinical variables are interdependent (e.g., Age influences ECOG score), provide a single joint proportion to guide the MaxEnt solver:

```R
my_constraints <- list(
  list(given = "age_ge65", target = "ECGO2_3", prob = 0.37)
)
```
## 3. Run the Pipeline
Execute the script to perform:

Phase 1: MaxEnt estimation of the joint distribution matrix.

Phase 2: SA optimization to generate the final multidimensional IPD table.

📈 Evaluation & Visualization
The output includes a high-fidelity "Digital Twin" cohort. You can validate the results using the integrated Nature-style RMST (Restricted Mean Survival Time) visualization module, which displays:

Shaded RMST areas for visual comparison.

Automated Hazard Ratio (HR) calculation via Cox Proportional Hazards.

Publication-ready aesthetics with Nature Publishing Group (NPG) color palettes.

📖 Citation
If you find this tool useful in your research, please cite our methodology:

Zheqing Zhu, et al. "Synthesizing multidimensional clinical profiles from published overall and stratified Kaplan–Meier images." (2026).

📄 License
This project is licensed under the MIT License. You are free to use, modify, and distribute the code for academic and commercial purposes, provided that original credit is given.
