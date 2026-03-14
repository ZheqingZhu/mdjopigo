# MD-JoPiGo
Synthesizing multidimensional clinical profiles from published overall and stratified Kaplan–Meier images.

MD-JoPiGo (Multidimensional Joint Patient Individual-data Generator and Optimizer) is a statistical framework designed to bridge the "dimensionality gap" in clinical trial reporting. It allows researchers to reconstruct high-dimensional individual patient data (IPD) directly from aggregated one-dimensional survival summaries (KM curves).

🔬 Methodology
The framework employs a two-stage inversion process:

Stage 1: Maximum Entropy Estimation: Estimates the most unbiased joint probability distribution of clinical characteristics given 1D marginal constraints.

Stage 2: Simulated Annealing Optimization: Swaps individual clinical labels to minimize the discrepancy (Weighted ISE) between reconstructed and target survival trajectories.

🌟 Highlights
Structural Identifiability: Handles parallel, mediated, and collider causal topologies.

Minimal Priors: Corrects coefficient drift in interdependent clinical networks (e.g., age-related frailty) using single intersectional proportions.

Publication-Ready: Generates high-fidelity digital twin cohorts with Nature-style RMST and HR visualizations.
