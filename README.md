# MD-JoPiGo: Multi-Dimensional Joint Patient Data Generator

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Statistical Computing](https://img.shields.io/badge/R-v4.2.2+-blue.svg)](https://www.r-project.org/)

**MD-JoPiGo** (Multidimensional Joint Patient Individual-data Generator and Optimizer) is an open-source statistical framework designed to bridge the **"dimensionality gap"** in clinical trial reporting. It synthesizes high-dimensional individual patient data (IPD) profiles from one-dimensional (1D) Kaplan-Meier survival summaries typically found in published literature.

---

## 🌟 The Workflow Overview

Transforming published survival plots into multidimensional digital twin cohorts is a two-step pipeline:
1. **IPD Extraction**: Convert Kaplan-Meier images into 1D-IPD using **KM-PoPiGo**.
2. **Multidimensional Synthesis**: Reconstruct the joint distribution and assign individual multivariable profiles using **MD-JoPiGo**.

---

## 🛠️ Step 1: Data Extraction via KM-PoPiGo (Pre-requisite)

Before running the MD-JoPiGo algorithm, you must extract the raw numerical survival data (time and status) from your target Kaplan-Meier curve images. 



We highly recommend using our sister tool, **[KM-PoPiGo](https://kmpopigo.github.io/)**, a web-based genetic optimization framework for high-fidelity IPD reconstruction. 

**How to prepare your data:**
1. **Overall Population**: Upload the KM image representing the entire treatment arm to KM-PoPiGo. Extract the IPD and save it.
2. **Stratified Subgroups**: Upload the KM images for all reported clinical subgroups (e.g., Male, Female, Age <65, Age >=65). Extract the IPD for each subgroup curve.

---

## 📊 Step 2: Data Formatting (`Overall.csv` & `Layer.csv`)

Once you have extracted the 1D-IPD from KM-PoPiGo, organize them into two specific CSV files placed in your working directory.

### 1. `Overall.csv` (The ITT Population)
This file contains the extracted IPD from the **Overall** KM curve. It provides the "pool" of survival times and event statuses for the entire cohort.

| Column | Type | Description |
| :--- | :--- | :--- |
| `time` | Numeric | The time-to-event or censoring time extracted from KM-PoPiGo. |
| `status` | Binary | Event indicator: `1` for event (death/progression), `0` for censored. |
| `group` | String | A uniform label for the overall population (e.g., `Overall` or `ITT`). |

### 2. `Layer.csv` (The Stratified Subgroups)
This file combines the extracted IPD from **all the stratified subgroup curves** into a single spreadsheet. The algorithm uses these 1D marginals to "anchor" the multidimensional optimization.

| Column | Type | Description |
| :--- | :--- | :--- |
| `time` | Numeric | Survival time. |
| `status` | Binary | Event indicator (`1` or `0`). |
| `group` | String | **Critical:** The label identifying the specific subgroup (e.g., `sex1`, `age_ge65`, `ECOG0_1`). These labels must strictly match the names you define in the R script's feature dictionary. |

---

## 🚀 Step 3: Run MD-JoPiGo

### 1. Configure Feature Mapping
Open the main R script and define your clinical dimensions in the `my_features` list. This maps the subgroup labels from `Layer.csv` to their respective multidimensional features.

```R
my_features <- list(
  Sex  = c("sex1", "sex2"),
  Age  = c("age_ge65", "age_ne65"),
  ECOG = c("ECGO0_1", "ECGO2_3")
)
```
## 2. Structural Calibration (Optional but Recommended for Mediators)
If your clinical variables are causally interdependent (e.g., Age influences ECOG performance status), provide a known joint proportion (structural prior) to guide the Maximum Entropy solver and correct coefficient drift:
```R
my_constraints <- list(
  list(given = "age_ge65", target = "ECGO2_3", prob = 0.37)
)
```
## 3. Execute the Pipeline
Run the script to automatically perform:

Phase 1: Maximum Entropy estimation of the joint distribution matrix.

Phase 2: Simulated Annealing optimization to generate the final multidimensional IPD (MD_JoPiGo_Final_IPD.csv).

📖 Citation
If you find this tool useful in your research, please cite our methodology:

Zheqing Zhu, et al. "Synthesizing multidimensional clinical profiles from published overall and stratified Kaplan–Meier images." (2026).

📄 License
This project is licensed under the MIT License. You are free to use, modify, and distribute the code for academic and commercial purposes, provided that original credit is given.
