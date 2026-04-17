# Step 9: Voxel-aware (VA) extension calibration
#
# Calibrates the per-DSB energy sensitivity parameter delta_f for each
# particle type (electron, proton, carbon) using the 2500-DSB calibration
# datasets from TOPAS-nBio simulations.
#
# Optimization via Pareto frontier: maximize within-nucleus CV while keeping
# population-mean error below 1%. Bootstrap uncertainty uses 1000 replicates.
#
# Inputs:  results/uvaom_v8_corrected_model.rds,
#          voxa_features_output_calibration/all_particles_calibration_energy_features.csv
# Outputs: results/voxa_voxel_aware_calibration.json,
#          results/voxa_voxel_aware_calibration.csv

library(tidyverse)
library(jsonlite)
set.seed(42)

# Create output directories if needed
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("figures")) dir.create("figures")


# Section 1: Load Voxa Om Base Model

cat("\n--- SECTION 1: LOADING VOxA OM BASE MODEL ---\n")
model <- readRDS("results/uvaom_v8_corrected_model.rds")
params <- model$parameters
FIXED_PARAMS <- model$fixed_params

cat(sprintf("Model version: %s\n", model$version))
cat(sprintf("R² (unweighted) = %.4f\n", model$fit_statistics$r2))
cat(sprintf("R² (weighted) = %.4f\n", model$fit_statistics$r2_weighted))
cat(sprintf("OER_max (theoretical, retention) = %.2f\n\n", model$OER_max_theoretical))


# Section 2: Extract Model Parameters

cat("\n--- SECTION 2: VOxA OM PARAMETERS ---\n")
# Oxygen kinetics
K_fix <- params$K_fix
K_repair <- params$K_repair

cat("Oxygen Kinetics (Michaelis-Menten):\n")
cat(sprintf("  K_fix    = %.4f%% O₂\n", K_fix))
cat(sprintf("  K_repair = %.4f%% O₂\n", K_repair))
cat(sprintf("  Ratio K_repair/K_fix = %.2f\n", K_repair / K_fix))
cat(sprintf("  Half-max O₂ ≈ %.2f mmHg\n\n", K_repair * 7.6))

# Case fractions (Sakata 2019 DSB combinatorics)
p1_low <- FIXED_PARAMS$p1_low
p2_low <- FIXED_PARAMS$p2_low
p3_low <- FIXED_PARAMS$p3_low
p1_high <- FIXED_PARAMS$p1_high

# Derive d and i fractions from p1_low = d²
d_frac <- sqrt(p1_low)
i_frac <- 1 - d_frac

cat("Sakata 2019 DSB Combinatorics:\n")
cat(sprintf("  d (single-hit direct)   = %.2f\n", d_frac))
cat(sprintf("  i (single-hit indirect) = %.2f\n", i_frac))
cat(sprintf("  p1_low  = d²  = %.4f (purely direct DSBs)\n", p1_low))
cat(sprintf("  p2_low  = 2di = %.4f (hybrid DSBs)\n", p2_low))
cat(sprintf("  p3_low  = i²  = %.4f (purely indirect DSBs)\n", p3_low))
cat(sprintf("  p1_high = %.2f (Nikjoo 2001, high-LET limit)\n\n", p1_high))

# Steepness parameters (hybrid light/heavy model)
s_dir_light <- params$s_dir_light
s_ind_light <- params$s_ind_light
s_dir_base <- params$s_dir_base
s_dir_scale <- params$s_dir_scale
s_ind_base <- params$s_ind_base
s_ind_scale <- params$s_ind_scale

cat("Steepness Parameters (Hybrid Model):\n")
cat(sprintf("  Light particles: s_dir = %.3f, s_ind = %.3f\n", s_dir_light, s_ind_light))
cat(sprintf("  Heavy (Z-scaled): s_dir(Z) = %.3f × (1 + %.3f × log(Z/2))\n", 
            s_dir_base, s_dir_scale))
cat(sprintf("                    s_ind(Z) = %.3f × (1 + %.3f × log(Z/2))\n\n",
            s_ind_base, s_ind_scale))

# Particle-specific x50 values from VOxA OM calibration
x50_dir_table <- c(
  electron = params$x50_dir_photon,  # Use photon as proxy for electron
  photon = params$x50_dir_photon,
  proton = params$x50_dir_proton,
  deuteron = params$x50_dir_deuteron,
  He = params$x50_dir_He,
  C = params$x50_dir_C,
  carbon = params$x50_dir_C,  # Alias
  Ne = params$x50_dir_Ne,
  Ar = params$x50_dir_Ar
)

x50_ind_table <- c(
  electron = params$x50_ind_photon,
  photon = params$x50_ind_photon,
  proton = params$x50_ind_proton,
  deuteron = params$x50_ind_deuteron,
  He = params$x50_ind_He,
  C = params$x50_ind_C,
  carbon = params$x50_ind_C,
  Ne = params$x50_ind_Ne,
  Ar = params$x50_ind_Ar
)

cat("Particle x50 Parameters (from VOxA OM):\n")
for (p in c("electron", "proton", "carbon")) {
  cat(sprintf("  %s: x50_dir = %.1f, x50_ind = %.1f\n", 
              p, x50_dir_table[p], x50_ind_table[p]))
}
# Theoretical OER_max
OER_max_theory <- model$OER_max_theoretical
OER_max_survival <- 1 + (OER_max_theory - 1) / 1.20

cat(sprintf("OER_max (theoretical):\n"))
cat(sprintf("  Retention: %.2f\n", OER_max_theory))
cat(sprintf("  Survival:  %.2f (conversion factor = 1.20)\n\n", OER_max_survival))

# VA calibration O2 level
O2_calibration <- 0.21  # Moderate hypoxia (% O2)

cat(sprintf("VA Calibration O2 level = %.2f%% (moderate hypoxia)\n\n", O2_calibration))

# Physical bounds for f_direct
f_min <- max(p1_low / 2, 0.02)  # Minimum direct fraction
f_max <- p1_high                 # Maximum direct fraction (Nikjoo 2001)

cat(sprintf("Physical bounds for f_direct: [%.2f, %.2f]\n\n", f_min, f_max))


# Section 3: Load 2500-Dsb Calibration Energy Features

cat("\n--- SECTION 3: LOADING 2500-DSB CALIBRATION ENERGY FEATURES ---\n")
# Try multiple possible paths for the energy features file
energy_features_paths <- c(
  "voxa_features_output_calibration/all_particles_calibration_energy_features.csv",
  "voxa_features_output_calibration/all_particles_energy_features.csv",
  "calibration_output/all_particles_calibration_energy_features.csv",
  "energy_features_output/all_particles_energy_features.csv"
)

energy_features_path <- NULL
for (path in energy_features_paths) {
  if (file.exists(path)) {
    energy_features_path <- path
    break
  }
}

if (is.null(energy_features_path)) {
  cat("ERROR: Energy features file not found!\n")
  cat("Searched paths:\n")
  for (path in energy_features_paths) {
    cat(sprintf("  - %s\n", path))
  }
  cat("\nPlease run extract_energy_features_voxa.py first:\n")
    stop("Energy features file not found")
}

# Load energy features
energy_data <- read_csv(energy_features_path, show_col_types = FALSE)

cat(sprintf("Energy features loaded from: %s\n", energy_features_path))
cat(sprintf("  Total DSBs: %d\n", nrow(energy_data)))
cat(sprintf("  Columns: %s\n", paste(names(energy_data)[1:8], collapse = ", ")))

# Check for required columns
required_cols <- c("particle", "E_local")
missing_cols <- setdiff(required_cols, names(energy_data))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

cat("\nParticles found:", paste(unique(energy_data$particle), collapse = ", "), "\n")
cat("\nDSBs per particle:\n")
dsb_counts <- energy_data %>% 
  count(particle) %>% 
  arrange(desc(n))
print(dsb_counts)
# E_local statistics by particle
cat("E_local statistics by particle:\n")
cat("─" %>% strrep(70), "\n")
energy_stats <- energy_data %>%
  group_by(particle) %>%
  summarise(
    n = n(),
    E_mean = mean(E_local),
    E_sd = sd(E_local),
    E_cv = sd(E_local) / mean(E_local) * 100,
    E_min = min(E_local),
    E_max = max(E_local),
    n_nonzero = sum(E_local > 0),
    .groups = "drop"
  )
print(energy_stats %>% mutate(across(where(is.numeric) & !matches("^n"), ~round(., 4))))
# Verify all DSBs have E_local > 0
n_zero_energy <- sum(energy_data$E_local == 0)
if (n_zero_energy > 0) {
  cat(sprintf("WARNING: %d DSBs have E_local = 0\n", n_zero_energy))
  } else {
  cat("✓ All DSBs have E_local > 0 (energy deposition present at each damage site)\n\n")
}


# Section 4: Calculate Particle-Specific Case Fractions

cat("\n--- SECTION 4: PARTICLE-SPECIFIC CASE FRACTIONS ---\n")
# Z-scaled steepness function
calc_steepness_Z <- function(Z, s_base, s_scale) {
  s_base * (1 + s_scale * log(max(Z, 2) / 2))
}

# Define particle physical properties for VA calibration
# These are the particles in our 2500-DSB calibration dataset
particles <- tibble(
  name = c("electron", "proton", "carbon"),
  Z = c(1, 1, 6),
  A_MeV = c(0.511, 938.3, 12 * 931.5),  # Rest mass energy
  E_kinetic = c(1, 10, 57 * 12),         # Kinetic energy in MeV (approximate)
  LET_keV_um = c(0.2, 4.6, 40.9),        # LET in keV/μm (from TOPAS simulations)
  particle_class = c("light", "light", "heavy")
) %>%
  mutate(
    # Relativistic calculations
    gamma = 1 + E_kinetic / A_MeV,
    beta = sqrt(1 - 1/gamma^2),
    
    # Effective charge (Barkas formula for ions)
    Z_eff = case_when(
      name == "carbon" ~ Z * (1 - exp(-125 * beta / Z^(2/3))),
      TRUE ~ as.numeric(Z)
    ),
    
    # Radiation quality parameter x = (Z_eff/β)²
    x = (Z_eff / beta)^2,
    
    # Particle-specific x50 values from VOxA OM
    x50_dir = x50_dir_table[name],
    x50_ind = x50_ind_table[name],
    
    # Get steepness based on particle class
    s_dir = case_when(
      particle_class == "light" ~ s_dir_light,
      TRUE ~ calc_steepness_Z(Z, s_dir_base, s_dir_scale)
    ),
    s_ind = case_when(
      particle_class == "light" ~ s_ind_light,
      TRUE ~ calc_steepness_Z(Z, s_ind_base, s_ind_scale)
    ),
    
    # LET transition functions (particle-specific)
    f_trans_dir = 1 / (1 + (x50_dir / pmax(x, 0.01))^s_dir),
    f_trans_ind = 1 / (1 + (x50_ind / pmax(x, 0.01))^s_ind),
    
    # Case fractions with LET transition
    # p1 increases with LET (more direct damage)
    p1 = p1_low + (p1_high - p1_low) * f_trans_dir,
    
    # p3 decreases with LET (less purely indirect damage)
    p3 = p3_low * (1 - f_trans_ind),
    
    # p2 fills the remainder (hybrid)
    p2_raw = 1 - p1 - p3,
    p2 = pmax(p2_raw, 0),  # Ensure non-negative
    
    # Renormalize to ensure sum = 1
    total = p1 + p2 + p3,
    p1 = p1 / total,
    p2 = p2 / total,
    p3 = p3 / total
  ) %>%
  select(-p2_raw, -total)

cat("Particle-Specific Parameters:\n")
cat("─" %>% strrep(80), "\n")
print(particles %>% 
        select(name, Z, LET_keV_um, x, x50_dir, s_dir, f_trans_dir, p1, p2, p3) %>%
        mutate(across(where(is.numeric), ~round(., 4))))
# Verify case fractions sum to 1
cat("Case fraction verification (should all be 1.0000):\n")
for (pname in particles$name) {
  pinfo <- particles %>% filter(name == pname)
  total <- pinfo$p1 + pinfo$p2 + pinfo$p3
  status <- ifelse(abs(total - 1.0) < 1e-10, "✓", "✗")
  cat(sprintf("  %s: p1 + p2 + p3 = %.6f %s\n", pname, total, status))
}
# Section 5: Core Physics Functions

cat("\n--- SECTION 5: CORE PHYSICS FUNCTIONS ---\n")
#' Calculate oxygen fixation probability (Michaelis-Menten)
#' p_ind = (O2 + K_fix) / (O2 + K_fix + K_repair)
calc_p_indirect <- function(O2) {
  (O2 + K_fix) / (O2 + K_fix + K_repair)
}

#' Calculate uniform model P_DSB (NORMALIZED to 21% O2)
#' P_DSB = P_raw(O2) / P_raw(21%)
#' This ensures P_DSB = 1.0 at normoxia
calc_P_DSB_uniform <- function(O2, p1, p2, p3) {
  p_ind <- calc_p_indirect(O2)
  p_ind_ref <- calc_p_indirect(21.0)  # Reference at normoxia
  
  # Raw retention probability
  P_raw <- p1 + p2 * p_ind + p3 * p_ind^2
  P_raw_ref <- p1 + p2 * p_ind_ref + p3 * p_ind_ref^2
  
  # Normalized (ensures P_DSB = 1.0 at 21% O2)
  P_raw / P_raw_ref
}

#' Calculate voxel-aware P_DSB (NORMALIZED, vectorized)
#' Each DSB gets its own f_direct based on local energy
#' f_direct(i) = p1_base + δf × E_zscore(i)
calc_P_DSB_voxel_aware <- function(O2, f_direct, p2_base, p3_base) {
  p_ind <- calc_p_indirect(O2)
  p_ind_ref <- calc_p_indirect(21.0)
  
  # Redistribute remaining probability mass between p2 and p3
  # Maintain ratio p2:p3 from base model
  remaining <- 1 - f_direct
  total_indirect_base <- p2_base + p3_base
  
  p2_local <- remaining * (p2_base / total_indirect_base)
  p3_local <- remaining * (p3_base / total_indirect_base)
  
  # Compute P_DSB (NORMALIZED)
  P_raw <- f_direct + p2_local * p_ind + p3_local * p_ind^2
  P_raw_ref <- f_direct + p2_local * p_ind_ref + p3_local * p_ind_ref^2
  
  P_raw / P_raw_ref
}

#' Calculate OER from P_DSB
calc_OER <- function(P_DSB) {
  1.0 / P_DSB
}

cat("Functions defined:\n")
# Verify normalization at 21% O2
cat("Normalization verification at 21% O2:\n")
for (pname in particles$name) {
  pinfo <- particles %>% filter(name == pname)
  P_DSB_21 <- calc_P_DSB_uniform(21.0, pinfo$p1, pinfo$p2, pinfo$p3)
  status <- ifelse(abs(P_DSB_21 - 1.0) < 1e-10, "✓", "✗")
  cat(sprintf("  %s: P_DSB(21%%) = %.10f %s\n", pname, P_DSB_21, status))
}
# Show OER predictions from uniform model
cat("Uniform model OER predictions (VOxA OM):\n")
O2_test <- c(21.0, 2.1, 0.21, 0.021, 0.001)
for (pname in particles$name) {
  pinfo <- particles %>% filter(name == pname)
  cat(sprintf("  %s (LET ≈ %.1f keV/μm):\n", toupper(pname), pinfo$LET_keV_um))
  for (O2 in O2_test) {
    P_DSB <- calc_P_DSB_uniform(O2, pinfo$p1, pinfo$p2, pinfo$p3)
    OER <- calc_OER(P_DSB)
    cat(sprintf("    O2=%6.3f%%: P_DSB=%.4f, OER=%.2f\n", O2, P_DSB, OER))
  }
}
# Section 6: Pareto Frontier Analysis

cat("\n--- SECTION 6: PARETO FRONTIER ANALYSIS ---\n")
#' Compute Pareto frontier for δf calibration
#' Trade-off between:
#'   - Maximizing heterogeneity (CV of P_DSB)
#'   - Preserving mean (small error vs uniform model)
compute_pareto_frontier <- function(edata, pinfo, n_points = 200) {
  n <- nrow(edata)
  p1 <- pinfo$p1
  p2 <- pinfo$p2
  p3 <- pinfo$p3
  
  # Standardize energy (z-scores)
  E_mean <- mean(edata$E_local)
  E_sd <- sd(edata$E_local)
  
  if (E_sd < 1e-10) {
    warning("E_local has zero variance - cannot compute Pareto frontier")
    return(NULL)
  }
  
  E_zscore <- (edata$E_local - E_mean) / E_sd
  
  # Target P_DSB from uniform model at calibration O2
  P_DSB_target <- calc_P_DSB_uniform(O2_calibration, p1, p2, p3)
  
  # Maximum δf that keeps f_direct within physical bounds for all DSBs
  score_range <- max(abs(E_zscore))
  if (score_range == 0) score_range <- 1
  
  delta_f_max <- min(
    (f_max - p1) / score_range,
    (p1 - f_min) / score_range
  ) * 0.95  # 5% safety margin
  
  # Ensure delta_f_max is positive
  delta_f_max <- max(delta_f_max, 0.001)
  
  # Test range of δf values
  delta_f_values <- seq(0, delta_f_max, length.out = n_points)
  
  results <- map_dfr(delta_f_values, function(delta_f) {
    # Calculate f_direct for each DSB (Equation 7)
    f_direct_raw <- p1 + delta_f * E_zscore
    f_direct <- pmin(pmax(f_direct_raw, f_min), f_max)
    
    # Calculate P_DSB for each DSB (Equation 10)
    P_DSB <- calc_P_DSB_voxel_aware(O2_calibration, f_direct, p2, p3)
    
    # Metrics
    tibble(
      delta_f = delta_f,
      f_direct_mean = mean(f_direct),
      f_direct_sd = sd(f_direct),
      f_direct_min = min(f_direct),
      f_direct_max = max(f_direct),
      P_DSB_mean = mean(P_DSB),
      P_DSB_sd = sd(P_DSB),
      P_DSB_cv = sd(P_DSB) / mean(P_DSB) * 100,
      mean_error_pct = abs(mean(P_DSB) - P_DSB_target) / P_DSB_target * 100,
      cor_E_P = cor(edata$E_local, P_DSB),
      n_clipped_low = sum(f_direct_raw < f_min),
      n_clipped_high = sum(f_direct_raw > f_max)
    )
  })
  
  results$P_DSB_target <- P_DSB_target
  results$particle <- pinfo$name
  results$n_dsbs <- n
  results$f_min_bound <- f_min
  results$f_max_bound <- f_max
  results$E_mean <- E_mean
  results$E_sd <- E_sd
  
  return(results)
}

cat("Pareto frontier optimization function defined.\n")
# Section 7: Select Optimal Operating Point

#' Select optimal δf from Pareto frontier
#' Maximize heterogeneity (CV) while keeping mean error < threshold
select_operating_point <- function(pareto_df, max_mean_error_pct = 1.0) {
  # Filter to points with acceptable mean error
  valid_points <- pareto_df %>%
    filter(mean_error_pct <= max_mean_error_pct)
  
  if (nrow(valid_points) == 0) {
    cat("    Warning: No points with error <", max_mean_error_pct, 
        "%. Relaxing constraint.\n")
    # Try progressively larger thresholds
    for (threshold in c(2.0, 5.0, 10.0)) {
      valid_points <- pareto_df %>%
        filter(mean_error_pct <= threshold)
      if (nrow(valid_points) > 0) {
        cat("    Using threshold:", threshold, "%\n")
        break
      }
    }
    
    if (nrow(valid_points) == 0) {
      # Use point with minimum error
      return(pareto_df %>%
               filter(mean_error_pct == min(mean_error_pct)) %>%
               slice(1))
    }
  }
  
  # Among valid points, choose one with maximum heterogeneity (P_DSB_cv)
  valid_points %>%
    filter(P_DSB_cv == max(P_DSB_cv)) %>%
    slice(1)
}


# Section 8: Bootstrap Uncertainty Quantification

cat("\n--- SECTION 8: BOOTSTRAP UNCERTAINTY QUANTIFICATION ---\n")
#' Bootstrap calibration for uncertainty quantification
#' Resample DSBs with replacement and re-optimize δf
bootstrap_calibrate <- function(edata, pinfo, n_bootstrap = 1000, 
                                max_mean_error_pct = 1.0,
                                show_progress = TRUE) {
  n <- nrow(edata)
  
  # Original calibration (full Pareto frontier)
  pareto_orig <- compute_pareto_frontier(edata, pinfo, n_points = 200)
  
  if (is.null(pareto_orig)) {
    warning("Could not compute Pareto frontier for ", pinfo$name)
    return(NULL)
  }
  
  selected_orig <- select_operating_point(pareto_orig, max_mean_error_pct)
  
  # Bootstrap iterations
  if (show_progress) {
    cat(sprintf("    Running %d bootstrap iterations...\n", n_bootstrap))
  }
  
  boot_results <- map_dfr(1:n_bootstrap, function(b) {
    # Resample DSBs with replacement
    idx <- sample(1:n, n, replace = TRUE)
    edata_boot <- edata[idx, ]
    
    # Compute Pareto frontier (fewer points for speed)
    pareto_boot <- compute_pareto_frontier(edata_boot, pinfo, n_points = 50)
    
    if (is.null(pareto_boot)) {
      return(tibble(delta_f = NA, P_DSB_sd = NA, P_DSB_cv = NA, 
                    mean_error_pct = NA, boot_id = b))
    }
    
    selected_boot <- select_operating_point(pareto_boot, max_mean_error_pct)
    
    selected_boot %>% 
      select(delta_f, P_DSB_sd, P_DSB_cv, mean_error_pct) %>%
      mutate(boot_id = b)
  }, .progress = FALSE)
  
  # Remove failed iterations
  boot_results <- boot_results %>% filter(!is.na(delta_f))
  
  if (nrow(boot_results) < n_bootstrap * 0.9) {
    warning(sprintf("Only %d/%d bootstrap iterations succeeded for %s",
                    nrow(boot_results), n_bootstrap, pinfo$name))
  }
  
  # Compute confidence intervals
  ci_summary <- boot_results %>%
    summarise(
      n_successful = n(),
      
      delta_f_mean = mean(delta_f),
      delta_f_se = sd(delta_f),
      delta_f_ci_low = quantile(delta_f, 0.025),
      delta_f_ci_high = quantile(delta_f, 0.975),
      
      P_DSB_sd_mean = mean(P_DSB_sd),
      P_DSB_sd_ci_low = quantile(P_DSB_sd, 0.025),
      P_DSB_sd_ci_high = quantile(P_DSB_sd, 0.975),
      
      P_DSB_cv_mean = mean(P_DSB_cv),
      P_DSB_cv_ci_low = quantile(P_DSB_cv, 0.025),
      P_DSB_cv_ci_high = quantile(P_DSB_cv, 0.975),
      
      mean_error_mean = mean(mean_error_pct),
      mean_error_ci_high = quantile(mean_error_pct, 0.975)
    )
  
  list(
    original = selected_orig,
    pareto = pareto_orig,
    bootstrap = boot_results,
    confidence_intervals = ci_summary,
    n_bootstrap = n_bootstrap,
    n_successful = nrow(boot_results)
  )
}

cat("Bootstrap calibration function defined.\n")
cat(sprintf("  Iterations: 1000\n"))
cat(sprintf("  Max mean error: 1.0%%\n\n"))


# Section 9: Main Calibration Loop

cat(sprintf("  SECTION 9: PARETO FRONTIER OPTIMIZATION (O2 = %.2f%%)\n", O2_calibration))
calibration_results <- list()
pareto_all <- list()
bootstrap_all <- list()

for (pname in particles$name) {
  pinfo <- particles %>% filter(name == pname)
  edata <- energy_data %>% filter(particle == pname)
  
  if (nrow(edata) == 0) {
    cat(sprintf("Skipping %s (no data)\n", pname))
    next
  }
  
  cat("─" %>% strrep(70), "\n")
  cat(sprintf("%s (n = %d DSBs)\n", toupper(pname), nrow(edata)))
  cat("─" %>% strrep(70), "\n")
  
  # Energy statistics
  cat(sprintf("  E_local: mean = %.4e, sd = %.4e, CV = %.1f%%\n",
              mean(edata$E_local), sd(edata$E_local),
              sd(edata$E_local) / mean(edata$E_local) * 100))
  cat(sprintf("  E_local range: [%.4e, %.4e]\n",
              min(edata$E_local), max(edata$E_local)))
  
  # Case fractions from VOxA OM
  cat(sprintf("  Case fractions (from OM): p1 = %.4f, p2 = %.4f, p3 = %.4f\n",
              pinfo$p1, pinfo$p2, pinfo$p3))
  
  # Uniform model target
  P_DSB_uniform <- calc_P_DSB_uniform(O2_calibration, pinfo$p1, pinfo$p2, pinfo$p3)
  OER_uniform <- calc_OER(P_DSB_uniform)
  cat(sprintf("  Uniform model at %.2f%% O2: P_DSB = %.4f, OER = %.2f\n", 
              O2_calibration, P_DSB_uniform, OER_uniform))
  
  # Bootstrap calibration (includes Pareto frontier)
  boot_result <- bootstrap_calibrate(edata, pinfo, n_bootstrap = 1000, 
                                     max_mean_error_pct = 1.0)
  
  if (is.null(boot_result)) {
        next
  }
  
  selected <- boot_result$original
  ci <- boot_result$confidence_intervals
  
  cat(sprintf("\n  OPTIMAL OPERATING POINT:\n"))
  cat(sprintf("    δf = %.5f [95%% CI: %.5f, %.5f]\n",
              selected$delta_f, ci$delta_f_ci_low, ci$delta_f_ci_high))
  cat(sprintf("    P_DSB CV = %.2f%% [95%% CI: %.2f%%, %.2f%%]\n",
              selected$P_DSB_cv, ci$P_DSB_cv_ci_low, ci$P_DSB_cv_ci_high))
  cat(sprintf("    Mean error = %.4f%% (target < 1%%)\n", selected$mean_error_pct))
  cat(sprintf("    f_direct range: [%.4f, %.4f]\n",
              selected$f_direct_min, selected$f_direct_max))
  cat(sprintf("    Cor(E_local, P_DSB) = %.4f\n\n", selected$cor_E_P))
  
  calibration_results[[pname]] <- selected
  pareto_all[[pname]] <- boot_result$pareto
  bootstrap_all[[pname]] <- boot_result
}


# Section 10: Validate At All O2 Levels

cat("\n--- SECTION 10: VALIDATION AT ALL O2 LEVELS ---\n")
O2_levels <- c(21.0, 2.1, 0.21, 0.021, 0.001)
O2_names <- c("Normoxia", "Mild hypoxia", "Moderate hypoxia", 
              "Severe hypoxia", "Anoxia")
validation_results <- list()

for (pname in particles$name) {
  pinfo <- particles %>% filter(name == pname)
  edata <- energy_data %>% filter(particle == pname)
  selected <- calibration_results[[pname]]
  
  if (is.null(selected) || nrow(edata) == 0) next
  
  # Compute f_direct for this particle using calibrated δf
  E_mean <- mean(edata$E_local)
  E_sd <- sd(edata$E_local)
  E_zscore <- (edata$E_local - E_mean) / E_sd
  
  f_direct <- pmin(pmax(pinfo$p1 + selected$delta_f * E_zscore, f_min), f_max)
  
  cat(sprintf("%s (δf = %.5f):\n", toupper(pname), selected$delta_f))
  cat(sprintf("  %8s %10s %12s %10s %10s %8s %10s\n",
              "O2 (%)", "Condition", "P_DSB±SD", "Uniform", "CV (%)", "OER", "Retained"))
  cat("  ", "─" %>% strrep(68), "\n", sep = "")
  
  # Test at each O2 level
  for (i in seq_along(O2_levels)) {
    O2 <- O2_levels[i]
    condition <- O2_names[i]
    
    P_DSB_voxel <- calc_P_DSB_voxel_aware(O2, f_direct, pinfo$p2, pinfo$p3)
    P_DSB_uniform <- calc_P_DSB_uniform(O2, pinfo$p1, pinfo$p2, pinfo$p3)
    OER_uniform <- calc_OER(P_DSB_uniform)
    
    validation_results[[paste0(pname, "_", O2)]] <- tibble(
      particle = pname,
      O2_percent = O2,
      condition = condition,
      P_DSB_uniform = P_DSB_uniform,
      OER_uniform = OER_uniform,
      P_DSB_voxel_mean = mean(P_DSB_voxel),
      P_DSB_voxel_sd = sd(P_DSB_voxel),
      P_DSB_voxel_cv = sd(P_DSB_voxel) / mean(P_DSB_voxel) * 100,
      mean_error_pct = abs(mean(P_DSB_voxel) - P_DSB_uniform) / P_DSB_uniform * 100,
      cor_E_P = cor(edata$E_local, P_DSB_voxel),
      n_dsbs = nrow(edata),
      expected_retained = round(mean(P_DSB_voxel) * nrow(edata))
    )
    
    cat(sprintf("  %8.3f %10s %6.4f±%.4f %10.4f %10.2f %8.2f %10d\n",
                O2, substr(condition, 1, 10),
                mean(P_DSB_voxel), sd(P_DSB_voxel), P_DSB_uniform,
                sd(P_DSB_voxel) / mean(P_DSB_voxel) * 100,
                OER_uniform,
                round(mean(P_DSB_voxel) * nrow(edata))))
  }
  }

validation_df <- bind_rows(validation_results)


# Section 11: Compile Final Calibration Summary

cat("\n--- SECTION 11: FINAL CALIBRATION SUMMARY ---\n")
final_calibration <- map_dfr(particles$name, function(pname) {
  pinfo <- particles %>% filter(name == pname)
  selected <- calibration_results[[pname]]
  boot <- bootstrap_all[[pname]]
  edata <- energy_data %>% filter(particle == pname)
  
  if (is.null(selected) || is.null(boot)) return(NULL)
  
  ci <- boot$confidence_intervals
  
  tibble(
    particle = pname,
    n_dsbs = nrow(edata),
    
    # Base case fractions (from VOxA OM)
    p1_base = pinfo$p1,
    p2_base = pinfo$p2,
    p3_base = pinfo$p3,
    
    # LET transition parameters
    LET_keV_um = pinfo$LET_keV_um,
    x50_dir = pinfo$x50_dir,
    x50_ind = pinfo$x50_ind,
    s_dir = pinfo$s_dir,
    s_ind = pinfo$s_ind,
    f_trans_dir = pinfo$f_trans_dir,
    f_trans_ind = pinfo$f_trans_ind,
    
    # Voxel-aware calibration results
    delta_f = selected$delta_f,
    delta_f_ci_low = ci$delta_f_ci_low,
    delta_f_ci_high = ci$delta_f_ci_high,
    
    # Physical bounds
    f_min = f_min,
    f_max = f_max,
    
    # Energy statistics
    E_mean = selected$E_mean,
    E_sd = selected$E_sd,
    
    # Operating metrics at calibration O2
    P_DSB_target = selected$P_DSB_target,
    P_DSB_sd = selected$P_DSB_sd,
    P_DSB_cv = selected$P_DSB_cv,
    P_DSB_cv_ci_low = ci$P_DSB_cv_ci_low,
    P_DSB_cv_ci_high = ci$P_DSB_cv_ci_high,
    mean_error_pct = selected$mean_error_pct,
    cor_E_P = selected$cor_E_P,
    
    # Effective f_direct range
    f_direct_min = selected$f_direct_min,
    f_direct_max = selected$f_direct_max,
    f_direct_range = selected$f_direct_max - selected$f_direct_min
  )
})

cat("Calibrated Voxel-Aware Parameters:\n")
cat("─" %>% strrep(80), "\n")
cat(sprintf("%-12s %10s %15s %25s %15s\n",
            "Particle", "N DSBs", "δf", "95% CI", "CV (%)"))
cat("─" %>% strrep(80), "\n")
for (i in 1:nrow(final_calibration)) {
  row <- final_calibration[i, ]
  cat(sprintf("%-12s %10d %15.5f [%10.5f, %10.5f] %10.2f\n",
              toupper(row$particle), row$n_dsbs, row$delta_f,
              row$delta_f_ci_low, row$delta_f_ci_high, row$P_DSB_cv))
}
cat("Effective f_direct Ranges:\n")
cat("─" %>% strrep(70), "\n")
print(final_calibration %>%
        select(particle, p1_base, f_direct_min, f_direct_max, f_direct_range) %>%
        mutate(range_over_base = f_direct_range / p1_base) %>%
        mutate(across(where(is.numeric), ~round(., 4))))
# Section 11A: Z-Indexed Δf Interpolation Framework
#
# WHY Z-INTERPOLATION, NOT CV(E_local)-BASED REPARAMETERISATION
# -------------------------------------------------------------
# A natural question is whether δf can be expressed as α × CV(E_local),
# making it fully derivable from each normoxic run without a separate
# calibration. Inspection of the calibration data rules this out:
#
#   Electron:  CV(E_local) = 20.8%,  δf = 0.00357  →  α = 0.0172
#   Proton:    CV(E_local) = 36.2%,  δf = 0.02588  →  α = 0.0715
#   Carbon:    CV(E_local) = 31.5%,  δf = 0.08753  →  α = 0.2778
#
# α varies 16× across three particles. The reason: CV(E_local) is
# non-monotonic with LET (proton > carbon because proton's narrow,
# stochastically placed track produces high relative energy variance
# despite lower LET). Meanwhile δf is monotonically increasing with Z.
# These two quantities are structurally decoupled.
#
# The physically correct generalisation is Z-INTERPOLATION — the same
# mechanism already validated for x50_dir and x50_ind (2.2% error for
# neon hold-out). δf is indexed by atomic number Z, interpolated
# log-linearly in Z-space. For any unseen particle (helium Z=2, neon Z=10,
# argon Z=18), δf is derived analytically with no additional simulations.
#
# For SOBP LET variants of the same particle (e.g. proton pSOBP vs dSOBP),
# a single δf per Z is used. The LET dependence of the direct fraction is
# carried entirely by p1(LET, Z) from the OM transition functions, which
# are already computed at runtime in compute_P_DSB. δf modulates the
# per-DSB spread around that p1 baseline; this spread is primarily a
# function of track structure geometry (Z), not the specific LET value
# within a particle's therapeutic range.
#
# IMPLEMENTATION
# --------------
# Z_eff assignments for interpolation:
#   electron → Z_interp = 0   (photon/electron, low-Z anchor)
#   proton   → Z_interp = 1
#   carbon   → Z_interp = 6
#
# For Z > 1, log-linear interpolation in log(Z)-space between the two
# nearest calibrated points. For Z = 0 (electron/photon), the electron
# value is used directly as the lower anchor.
#

cat("\n--- SECTION 11A: Z-INDEXED δf INTERPOLATION FRAMEWORK ---\n")
# Build Z-indexed δf table from calibration results
# electron is assigned Z_interp = 0 as the lower anchor
Z_interp_map <- c(electron = 0, proton = 1, carbon = 6)

Z_delta_f_table <- final_calibration %>%
  mutate(Z_interp = Z_interp_map[particle]) %>%
  select(particle, Z_interp, delta_f, delta_f_ci_low, delta_f_ci_high, P_DSB_cv) %>%
  arrange(Z_interp)

cat("Z-indexed δf calibration table:\n")
cat("─" %>% strrep(70), "\n")
print(Z_delta_f_table %>% mutate(across(where(is.numeric) & !matches("Z"), ~round(., 5))))
#' Interpolate δf by atomic number Z using log-linear interpolation in Z-space.
#'
#' For Z > 0: log-linear interpolation between the two bracketing calibrated
#' particles in log(Z)-space (same method as OM x50 Z-interpolation).
#' For Z = 0 (electron/photon): returns the electron anchor value directly.
#' For Z > max(calibrated Z): linear extrapolation from the top two points.
#'
#' @param Z_target  Atomic number of the target particle (numeric).
#' @param Z_tbl     Z-indexed δf table (from Z_delta_f_table).
#' @return          Named list: delta_f, delta_f_ci_low, delta_f_ci_high, method.
interp_delta_f_by_Z <- function(Z_target, Z_tbl = Z_delta_f_table) {
  
  Z_vals    <- Z_tbl$Z_interp
  df_vals   <- Z_tbl$delta_f
  df_lo     <- Z_tbl$delta_f_ci_low
  df_hi     <- Z_tbl$delta_f_ci_high
  particles <- Z_tbl$particle
  
  # --- Exact match ---
  if (Z_target %in% Z_vals) {
    idx <- which(Z_vals == Z_target)
    return(list(
      delta_f        = df_vals[idx],
      delta_f_ci_low = df_lo[idx],
      delta_f_ci_high = df_hi[idx],
      method         = "exact",
      particle_source = particles[idx]
    ))
  }
  
  # --- Z = 0 (electron/photon anchor) ---
  if (Z_target == 0) {
    idx <- which(Z_vals == 0)
    return(list(
      delta_f        = df_vals[idx],
      delta_f_ci_low = df_lo[idx],
      delta_f_ci_high = df_hi[idx],
      method         = "electron_anchor"
    ))
  }
  
  # --- Interpolation/extrapolation for Z > 0 ---
  # Use only Z > 0 calibrated points for log-space interpolation
  pos_mask  <- Z_vals > 0
  Z_pos     <- Z_vals[pos_mask]
  df_pos    <- df_vals[pos_mask]
  df_lo_pos <- df_lo[pos_mask]
  df_hi_pos <- df_hi[pos_mask]
  
  log_Z_target <- log(Z_target)
  log_Z_pos    <- log(Z_pos)
  
  if (Z_target < min(Z_pos)) {
    # Below minimum Z > 0 calibrated point: interpolate between electron (Z=0)
    # treated as log(0.5) for numerical continuity and proton (Z=1)
    idx_e  <- which(Z_vals == 0)
    idx_p  <- which(Z_pos == min(Z_pos))
    log_Z_lo <- log(0.5)   # pseudo-log for electron anchor
    log_Z_hi <- log_Z_pos[idx_p]
    frac <- (log_Z_target - log_Z_lo) / (log_Z_hi - log_Z_lo)
    frac <- max(0, min(1, frac))  # clamp
    
    interp_log <- function(a, b) exp(log(a) + frac * (log(b) - log(a)))
    
    return(list(
      delta_f        = interp_log(df_vals[idx_e], df_pos[idx_p]),
      delta_f_ci_low = interp_log(df_lo[idx_e], df_lo_pos[idx_p]),
      delta_f_ci_high = interp_log(df_hi[idx_e], df_hi_pos[idx_p]),
      method         = "log_interp_low",
      Z_bracket      = c(0, Z_pos[idx_p])
    ))
  }
  
  if (Z_target > max(Z_pos)) {
    # Above maximum calibrated Z: extrapolate from top two positive-Z points
    ord    <- order(Z_pos)
    n_pos  <- length(Z_pos)
    idx_lo <- ord[n_pos - 1]
    idx_hi <- ord[n_pos]
    log_Z_lo <- log_Z_pos[idx_lo]
    log_Z_hi <- log_Z_pos[idx_hi]
    frac <- (log_Z_target - log_Z_lo) / (log_Z_hi - log_Z_lo)
    
    interp_log <- function(a, b) exp(log(a) + frac * (log(b) - log(a)))
    
    return(list(
      delta_f        = interp_log(df_pos[idx_lo], df_pos[idx_hi]),
      delta_f_ci_low = interp_log(df_lo_pos[idx_lo], df_lo_pos[idx_hi]),
      delta_f_ci_high = interp_log(df_hi_pos[idx_lo], df_hi_pos[idx_hi]),
      method         = "log_extrap",
      Z_bracket      = c(Z_pos[idx_lo], Z_pos[idx_hi])
    ))
  }
  
  # --- Standard bracket interpolation ---
  idx_lo <- max(which(Z_pos <= Z_target))
  idx_hi <- min(which(Z_pos >= Z_target))
  
  log_Z_lo <- log_Z_pos[idx_lo]
  log_Z_hi <- log_Z_pos[idx_hi]
  frac <- (log_Z_target - log_Z_lo) / (log_Z_hi - log_Z_lo)
  
  interp_log <- function(a, b) exp(log(a) + frac * (log(b) - log(a)))
  
  list(
    delta_f        = interp_log(df_pos[idx_lo], df_pos[idx_hi]),
    delta_f_ci_low = interp_log(df_lo_pos[idx_lo], df_lo_pos[idx_hi]),
    delta_f_ci_high = interp_log(df_hi_pos[idx_lo], df_hi_pos[idx_hi]),
    method         = "log_interp",
    Z_bracket      = c(Z_pos[idx_lo], Z_pos[idx_hi])
  )
}

# ── Validate on calibrated particles (should recover exact values) ──────────
cat("Interpolation validation on calibrated particles (should be exact):\n")
for (pname in c("electron", "proton", "carbon")) {
  Z_p    <- Z_interp_map[pname]
  result <- interp_delta_f_by_Z(Z_p)
  cal_df <- final_calibration %>% filter(particle == pname) %>% pull(delta_f)
  err    <- abs(result$delta_f - cal_df) / cal_df * 100
  cat(sprintf("  %-10s Z=%d: δf = %.6f (calibrated = %.6f, error = %.4f%%) [%s]\n",
              pname, Z_p, result$delta_f, cal_df, err, result$method))
}
# ── Predict unseen particles ──────────────────────────────────────────────
unseen_particles <- tibble(
  name = c("helium", "nitrogen", "neon",  "argon"),
  Z    = c(2,        7,          10,       18)
)

cat("δf predictions for unseen particles (Z-interpolation):\n")
cat("─" %>% strrep(70), "\n")
cat(sprintf("  %-12s  %4s  %10s  %25s  %-20s\n",
            "Particle", "Z", "δf", "95% CI", "Method"))
cat("─" %>% strrep(70), "\n")

interp_predictions <- map_dfr(1:nrow(unseen_particles), function(i) {
  p   <- unseen_particles[i, ]
  res <- interp_delta_f_by_Z(p$Z)
  cat(sprintf("  %-12s  %4d  %10.6f  [%.6f, %.6f]  %s\n",
              p$name, p$Z, res$delta_f, res$delta_f_ci_low,
              res$delta_f_ci_high, res$method))
  tibble(
    particle        = p$name,
    Z_interp        = p$Z,
    delta_f         = res$delta_f,
    delta_f_ci_low  = res$delta_f_ci_low,
    delta_f_ci_high = res$delta_f_ci_high,
    method          = res$method
  )
})
# ── Monotonicity check across Z ──────────────────────────────────────────
cat("Monotonicity check (δf should increase with Z):\n")
all_Z <- tibble(
  particle = c("electron", "proton", "helium", "carbon", "neon", "argon"),
  Z        = c(0, 1, 2, 6, 10, 18)
) %>%
  mutate(delta_f = map_dbl(Z, ~interp_delta_f_by_Z(.x)$delta_f))

print(all_Z %>% mutate(delta_f = round(delta_f, 6)))
is_monotone <- all(diff(all_Z$delta_f) > 0)
cat(sprintf("\n  Monotone increasing: %s\n\n",
            ifelse(is_monotone, "✓ YES", "✗ NO — check interpolation")))

# Export Z-δf table for downstream use
write_csv(Z_delta_f_table, "results/voxa_Z_delta_f_table.csv")
write_csv(interp_predictions, "results/voxa_Z_delta_f_interpolated.csv")
# Section 12: Save Results

cat("\n--- SECTION 12: SAVING RESULTS ---\n")
# CSV exports
write_csv(final_calibration, "results/voxa_voxel_aware_calibration.csv")
write_csv(validation_df, "results/voxa_voxel_aware_validation.csv")
# Pareto frontiers
pareto_combined <- bind_rows(pareto_all)
write_csv(pareto_combined, "results/voxa_voxel_aware_pareto_frontiers.csv")
# JSON export for Python integration
json_export <- list(
  model_version = "VOxA v1.0 (Voxel-Aware Oxygen Model)",
  calibration_type = "Voxel-Aware (VA)",
  base_model_version = model$version,
  base_model_r2 = model$fit_statistics$r2,
  base_model_r2_weighted = model$fit_statistics$r2_weighted,
  calibration_date = as.character(Sys.time()),
  
  base_parameters = list(
    K_fix = K_fix,
    K_repair = K_repair,
    d_fraction = d_frac,
    i_fraction = i_frac,
    p1_low = p1_low,
    p2_low = p2_low,
    p3_low = p3_low,
    p1_high = p1_high,
    s_dir_light = s_dir_light,
    s_ind_light = s_ind_light,
    s_dir_base = s_dir_base,
    s_dir_scale = s_dir_scale,
    s_ind_base = s_ind_base,
    s_ind_scale = s_ind_scale,
    OER_max_retention = OER_max_theory,
    OER_max_survival = OER_max_survival,
    conversion_factor = 1.20
  ),
  
  va_calibration_settings = list(
    O2_calibration_pct = O2_calibration,
    f_min = f_min,
    f_max = f_max,
    description = "Pareto frontier optimization with bootstrap uncertainty",
    n_bootstrap = 1000,
    max_mean_error_pct = 1.0
  ),
  
  particle_calibration = final_calibration %>%
    split(.$particle) %>%
    map(function(x) {
      list(
        n_dsbs = x$n_dsbs,
        LET_keV_um = x$LET_keV_um,
        p1_base = x$p1_base,
        p2_base = x$p2_base,
        p3_base = x$p3_base,
        x50_dir = x$x50_dir,
        x50_ind = x$x50_ind,
        s_dir = x$s_dir,
        s_ind = x$s_ind,
        delta_f = x$delta_f,
        delta_f_ci_low = x$delta_f_ci_low,
        delta_f_ci_high = x$delta_f_ci_high,
        E_mean = x$E_mean,
        E_sd = x$E_sd,
        P_DSB_cv = x$P_DSB_cv,
        P_DSB_cv_ci_low = x$P_DSB_cv_ci_low,
        P_DSB_cv_ci_high = x$P_DSB_cv_ci_high,
        mean_error_pct = x$mean_error_pct,
        f_direct_min = x$f_direct_min,
        f_direct_max = x$f_direct_max
      )
    }),
  
  validation = validation_df %>%
    split(paste(.$particle, .$O2_percent, sep = "_")) %>%
    map(~as.list(.x)),

  # Z-indexed δf interpolation framework
  # Use interp_delta_f_by_Z(Z) for any unseen particle.
  # Z_interp assignments: electron=0, proton=1, carbon=6.
  Z_delta_f_table = Z_delta_f_table %>%
    split(.$particle) %>%
    map(function(x) list(
      particle        = x$particle,
      Z_interp        = x$Z_interp,
      delta_f         = x$delta_f,
      delta_f_ci_low  = x$delta_f_ci_low,
      delta_f_ci_high = x$delta_f_ci_high,
      P_DSB_cv        = x$P_DSB_cv
    )),

  Z_delta_f_interpolated_particles = interp_predictions %>%
    split(.$particle) %>%
    map(function(x) list(
      particle        = x$particle,
      Z_interp        = x$Z_interp,
      delta_f         = x$delta_f,
      delta_f_ci_low  = x$delta_f_ci_low,
      delta_f_ci_high = x$delta_f_ci_high,
      method          = x$method
    )),

  Z_interpolation_notes = list(
    description = paste(
      "δf is indexed by atomic number Z and interpolated log-linearly in",
      "log(Z)-space between the two nearest calibrated particles.",
      "Z_interp assignments: electron=0 (anchor), proton=1, carbon=6.",
      "For SOBP LET variants of the same particle (pSOBP vs dSOBP),",
      "use the same δf; LET dependence is captured by p1(LET,Z) from OM."
    ),
    Z_interp_assignments = list(electron = 0L, proton = 1L, carbon = 6L,
                                helium = 2L, nitrogen = 7L,
                                neon = 10L, argon = 18L),
    why_not_CV_E_local = paste(
      "CV(E_local) implies alpha values of 0.017 (electron), 0.072 (proton),",
      "0.278 (carbon) — a 16x range. Not a universal constant.",
      "CV(E_local) is non-monotonic with LET; δf is monotonic with Z.",
      "Z-interpolation is the physically consistent choice."
    )
  )
)

write_json(json_export, "results/voxa_voxel_aware_calibration.json",
           pretty = TRUE, auto_unbox = TRUE)
# Final Summary

cat(sprintf("║   Base Model: VOxA OM (R² = %.4f)                                 ║\n", 
            model$fit_statistics$r2))
cat(sprintf("║   K_fix = %.4f%% O₂  |  K_repair = %.4f%% O₂                     ║\n", 
            K_fix, K_repair))
cat(sprintf("║   OER_max: %.2f (retention) | %.2f (survival)                      ║\n", 
            OER_max_theory, OER_max_survival))
cat(sprintf("║   Case fractions (low LET): p1=%.2f, p2=%.2f, p3=%.2f            ║\n",
            p1_low, p2_low, p3_low))
for (pname in particles$name) {
  cal <- final_calibration %>% filter(particle == pname)
  if (nrow(cal) == 0) next
  
  cat(sprintf("║   %-10s: δf = %.5f, CV = %.2f%% [%.2f%%, %.2f%%]          ║\n",
              toupper(pname), cal$delta_f, cal$P_DSB_cv,
              cal$P_DSB_cv_ci_low, cal$P_DSB_cv_ci_high))
}

cat("\nStep 9 complete.\n\n")