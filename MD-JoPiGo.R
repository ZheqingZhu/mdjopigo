# =====================================================================
# MD-JoPiGo: Multi-Dimensional Joint Patient Data Generator
# =====================================================================
# This script integrates Maximum Entropy (MaxEnt) joint distribution 
# reconstruction and Simulated Annealing (SA) multidimensional optimization.
# =====================================================================

# ---------------------------------------------------------------------
# 0. INITIALIZATION (Install & Load Packages)
# ---------------------------------------------------------------------
if (!require("nloptr")) install.packages("nloptr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("survival")) install.packages("survival")

library(nloptr)
library(dplyr)
library(survival)

set.seed(2026) # Set seed for reproducibility

# =====================================================================
# 1. USER CONFIGURATION AREA (Users only need to modify this section)
# =====================================================================

# 1.1 Working Directory and File Paths
setwd("D:/paper/SCA/modle/opensource")      # Change to your folder path
overall_file <- "Overall.csv"        # Input 1: The overall 1D ITT IPD
layer_file   <- "Layer.csv"          # Input 2: The combined marginal 1D IPDs
output_file  <- "MD_JoPiGo_Final_IPD.csv" # Final Output IPD file

# 1.2 Dynamic Feature Dictionary
# Tell the algorithm which labels belong to the same feature dimension
my_features <- list(
  Sex  = c("sex1", "sex2"),
  Age  = c("age_ge65", "age_ne65"),
  ECOG = c("ECGO0_1", "ECGO2_3")
)

# 1.3 Structural Calibration Constraints (Prior Joint Probabilities)
# Example: P(Target | Given) = prob. Leave empty list() if no prior exists.
my_constraints <- list(
  #list(given = "age_ge65", target = "ECGO2_3", prob = 0.3733333),
  list(given = "age_ge65", target = "sex2",    prob = 0.3066667),
  list(given = "ECGO2_3", target = "sex2",    prob = 0.4102564)
)

# 1.4 Simulated Annealing Hyperparameters
sa_iterations   <- 50000
sa_initial_temp <- 2.0
sa_cooling_rate <- 0.9998


# =====================================================================
# 2. CORE ALGORITHM MODULES (Users DO NOT need to modify below)
# =====================================================================

# --- MODULE A: Maximum Entropy & Structural Calibration ---
run_max_entropy <- function(df_overall, df_layer, feature_mapping, extra_constraints) {
  cat("\n==================================================\n")
  cat("[PHASE 1] Structural Calibration & MaxEnt Reconstruction\n")
  cat("==================================================\n")
  
  total_n <- nrow(df_overall)
  subgroups_grid <- expand.grid(feature_mapping, stringsAsFactors = FALSE)
  n_subgroups <- nrow(subgroups_grid)
  subgroups_labels <- apply(subgroups_grid, 1, paste, collapse = " | ")
  
  marginals <- list()
  for (feat_name in names(feature_mapping)) {
    levels <- feature_mapping[[feat_name]]
    for (i in 1:(length(levels) - 1)) {
      lvl <- levels[i]
      count <- sum(df_layer$group == lvl, na.rm = TRUE)
      marginals[[lvl]] <- count / total_n
    }
  }
  
  n_constraints <- 1 + length(marginals) + length(extra_constraints)
  A <- matrix(0, nrow = n_constraints, ncol = n_subgroups)
  b <- numeric(n_constraints)
  
  A[1, ] <- 1; b[1] <- 1.0
  row_idx <- 2
  
  for (lvl in names(marginals)) {
    idx <- apply(subgroups_grid, 1, function(row) lvl %in% row)
    A[row_idx, idx] <- 1; b[row_idx] <- marginals[[lvl]]
    row_idx <- row_idx + 1
  }
  
  for (cond in extra_constraints) {
    idx_given <- apply(subgroups_grid, 1, function(r) cond$given %in% r)
    idx_both  <- apply(subgroups_grid, 1, function(r) (cond$given %in% r) && (cond$target %in% r))
    A[row_idx, idx_given] <- -cond$prob
    A[row_idx, idx_both]  <- A[row_idx, idx_both] + 1
    b[row_idx] <- 0
    cat(sprintf("-> Applied Prior: P(%s | %s) = %.2f\n", cond$target, cond$given, cond$prob))
    row_idx <- row_idx + 1
  }
  
  eval_f_penalty <- function(x) { 
    x_safe <- pmax(x, 1e-10)
    entropy <- sum(x_safe * log(x_safe)) 
    penalty <- 10000 * sum((A %*% x - b)^2)
    return(entropy + penalty)
  }
  
  x0 <- rep(1 / n_subgroups, n_subgroups)
  
  res <- slsqp(x0, fn = eval_f_penalty, lower = rep(0, n_subgroups), upper = rep(1, n_subgroups))
  
  if (res$convergence >= 0) {
    cat(sprintf("-> MaxEnt Penalty Value: %.6f\n", sum((A %*% res$par - b)^2)))
    
    exact_counts <- round(res$par * total_n)

    diff <- total_n - sum(exact_counts)
    if (diff != 0) {
      max_idx <- which.max(exact_counts)
      exact_counts[max_idx] <- exact_counts[max_idx] + diff
    }
    
    cat("-> Phase 1 Solved! Topology matrix generated via Lagrangian Relaxation.\n")
    return(list(df = data.frame(Subgroup = subgroups_labels, Exact_Count = exact_counts), grid = subgroups_grid))
  } else {
    stop("Optimization failed.")
  }
}

# --- MODULE B: High-Speed Simulated Annealing ---
fast_km_at_grid <- function(t, e, grid) {
  if (length(t) == 0) return(rep(1, length(grid)))
  ord <- order(t); t <- t[ord]; e <- e[ord]; n <- length(t)
  hazards <- e / (n:1); surv <- cumprod(1 - hazards)
  sf <- stepfun(t, c(1, surv), right = FALSE); return(sf(grid))
}

run_simulated_annealing <- function(df_overall, df_layer, maxent_data, feature_mapping, n_iters, temp, cooling_rate) {
  cat("\n==================================================\n")
  cat("[PHASE 2] Multi-Dimensional Simulated Annealing\n")
  cat("==================================================\n")
  
  df_overall_treat <- df_overall %>% arrange(time)
  overall_times <- df_overall_treat$time
  overall_events <- df_overall_treat$status
  N <- length(overall_times)
  time_grid <- seq(0, max(overall_times), length.out = 100)
  
  target_S <- list()
  for (feat_name in names(feature_mapping)) {
    for (lvl in feature_mapping[[feat_name]]) {
      df_sub <- df_layer %>% filter(group == lvl)
      fit <- survfit(Surv(time, status) ~ 1, data = df_sub)
      target_S[[lvl]] <- stepfun(fit$time, c(1, fit$surv), right = FALSE)(time_grid)
    }
  }
  
  K <- length(feature_mapping)
  labels_pool <- matrix(nrow = N, ncol = K)
  idx <- 1
  for (i in 1:nrow(maxent_data$df)) {
    if (maxent_data$df$Exact_Count[i] > 0) {
      int_vals <- sapply(1:K, function(k) match(as.character(maxent_data$grid[i, k]), feature_mapping[[k]]))
      for (c in 1:maxent_data$df$Exact_Count[i]) { labels_pool[idx, ] <- int_vals; idx <- idx + 1 }
    }
  }
  
  current_labels <- labels_pool[sample(N), ]
  
  calculate_loss <- function(labels) {
    loss <- 0
    for (k in 1:K) {
      for (lvl_idx in 1:length(feature_mapping[[k]])) {
        mask <- labels[, k] == lvl_idx
        if (sum(mask) == 0) loss <- loss + 1000 
        else loss <- loss + sum((fast_km_at_grid(overall_times[mask], overall_events[mask], time_grid) - target_S[[feature_mapping[[k]][lvl_idx]]])^2)
      }
    }
    return(loss)
  }
  
  current_loss <- calculate_loss(current_labels)
  best_labels <- current_labels; best_loss <- current_loss
  cat(sprintf("-> Initial Baseline Loss: %.4f\n", current_loss))
  
  for (i in 1:n_iters) {
    idx_swap <- sample(1:N, 2)
    new_labels <- current_labels
    new_labels[idx_swap[1], ] <- current_labels[idx_swap[2], ]
    new_labels[idx_swap[2], ] <- current_labels[idx_swap[1], ]
    
    new_loss <- calculate_loss(new_labels)
    if (new_loss < current_loss || runif(1) < exp(-(new_loss - current_loss) / temp)) {
      current_labels <- new_labels; current_loss <- new_loss
      if (new_loss < best_loss) { best_labels <- current_labels; best_loss <- new_loss }
    }
    temp <- temp * cooling_rate
    if (i %% 10000 == 0 || i == n_iters) cat(sprintf("   Iter %5d | Temp: %.4f | Best Loss: %.4f\n", i, temp, best_loss))
  }
  
  for (k in 1:K) df_overall_treat[[names(feature_mapping)[k]]] <- feature_mapping[[k]][best_labels[, k]]
  return(df_overall_treat %>% select(time, status, all_of(names(feature_mapping))))
}

# =====================================================================
# 3. EXECUTE PIPELINE
# =====================================================================
if (!file.exists(overall_file) || !file.exists(layer_file)) stop("Input CSV files not found! Please check file paths.")

df_overall <- read.csv(overall_file, stringsAsFactors = FALSE)
df_layer   <- read.csv(layer_file, stringsAsFactors = FALSE)

# Run Phase 1
maxent_result <- run_max_entropy(df_overall, df_layer, my_features, my_constraints)

# Run Phase 2
final_ipd <- run_simulated_annealing(df_overall, df_layer, maxent_result, my_features, sa_iterations, sa_initial_temp, sa_cooling_rate)

# Save Results
write.csv(final_ipd, output_file, row.names = FALSE)
cat("\n==================================================\n")
cat("🎉 SUCCESS! Multi-dimensional super IPD has been generated.\n")
cat(sprintf("💾 File saved to: %s\n", file.path(getwd(), output_file)))
cat("==================================================\n")
