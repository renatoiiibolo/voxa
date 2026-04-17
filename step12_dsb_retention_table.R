# Step 12: DSB retention table
#
# Computes expected retained DSB counts per nucleus across five oxygen
# conditions using the VA model and the 400-DSB validation datasets.
# Starting from ~400 initial DSBs (≈2 Gy carbon pSOBP), reports counts
# for electron, proton, and carbon with 95% Monte Carlo confidence intervals.
#
# Inputs:  results/voxa_voxel_aware_calibration.json,
#          voxa_features_output_validation/all_particles_validation_energy_features.csv
# Outputs: results/voxa_dsb_retention_table.csv,
#          results/voxa_dsb_retention_detailed.csv,
#          figures/fig47_dsb_retention.png

library(tidyverse)
library(jsonlite)
library(gridExtra)
library(grid)

set.seed(42)

# Create output directories
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("figures")) dir.create("figures")


mediterranean_blue <- "#1E5B8C"
pompeii_red <- "#C75B5B"
terracotta <- "#E07B3C"
lemon_yellow <- "#D4A84B"
olive_green <- "#2D7D46"
warm_stone <- "#8B7355"
sandy_beige <- "#F5E6D3"
ancient_stone <- "#4A3728"

particle_colors <- c(
  "electron" = mediterranean_blue,
  "proton" = lemon_yellow,
  "carbon" = pompeii_red
)


# Section 1: Load Voxa Calibration

cat("\n--- SECTION 1: LOADING VOxA CALIBRATION (Step 9) ---\n")
calibration_paths <- c(
  "results/voxa_voxel_aware_calibration.json",
  "voxa_voxel_aware_calibration.json"
)

calibration_file <- NULL
for (path in calibration_paths) {
  if (file.exists(path)) {
    calibration_file <- path
    break
  }
}

if (is.null(calibration_file)) {
  stop("ERROR: VOxA calibration file not found!\n",
       "Please run Step 9 first: Rscript step9_voxa_voxel_aware_calibration.R")
}

cat(sprintf("Loading calibration from: %s\n", calibration_file))
calibration <- fromJSON(calibration_file)

# Extract parameters
VOXA_PARAMS <- list(
  model_version = calibration$model_version,
  base_model_version = calibration$base_model_version,
  K_fix = calibration$base_parameters$K_fix,
  K_repair = calibration$base_parameters$K_repair,
  O2_normoxia = 21.0,
  O2_calibration = calibration$va_calibration_settings$O2_calibration_pct,
  f_min = calibration$va_calibration_settings$f_min,
  f_max = calibration$va_calibration_settings$f_max,
  p1_low = calibration$base_parameters$p1_low,
  p2_low = calibration$base_parameters$p2_low,
  p3_low = calibration$base_parameters$p3_low,
  p1_high = calibration$base_parameters$p1_high,
  OER_max_retention = calibration$base_parameters$OER_max_retention,
  OER_max_survival = calibration$base_parameters$OER_max_survival,
  conversion_factor = calibration$base_parameters$conversion_factor,
  R_squared = calibration$base_model_r2,
  R_squared_weighted = calibration$base_model_r2_weighted,
  particles = list()
)

# Load particle-specific parameters
PARTICLE_KEYS <- names(calibration$particle_calibration)
for (p in PARTICLE_KEYS) {
  pc <- calibration$particle_calibration[[p]]
  VOXA_PARAMS$particles[[p]] <- list(
    n_dsbs_calibration = pc$n_dsbs,
    LET_keV_um = pc$LET_keV_um,
    p1 = pc$p1_base,
    p2 = pc$p2_base,
    p3 = pc$p3_base,
    delta_f = pc$delta_f,
    delta_f_ci_low = pc$delta_f_ci_low,
    delta_f_ci_high = pc$delta_f_ci_high,
    P_DSB_cv = pc$P_DSB_cv,
    E_mean = pc$E_mean,
    E_sd = pc$E_sd
  )
}

cat(sprintf("\nModel version: %s\n", VOXA_PARAMS$model_version))
cat(sprintf("Base model R²: %.4f (weighted: %.4f)\n", 
            VOXA_PARAMS$R_squared, VOXA_PARAMS$R_squared_weighted))

cat("\nOxygen Kinetics:\n")
cat(sprintf("  K_fix    = %.4f%% O₂\n", VOXA_PARAMS$K_fix))
cat(sprintf("  K_repair = %.4f%% O₂\n", VOXA_PARAMS$K_repair))
cat(sprintf("  OER_max  = %.2f (retention) / %.2f (survival)\n\n",
            VOXA_PARAMS$OER_max_retention, VOXA_PARAMS$OER_max_survival))

cat("Case Fractions (Low LET):\n")
cat(sprintf("  p1_low  = %.2f (purely direct DSBs)\n", VOXA_PARAMS$p1_low))
cat(sprintf("  p2_low  = %.2f (hybrid DSBs)\n", VOXA_PARAMS$p2_low))
cat(sprintf("  p3_low  = %.2f (purely indirect DSBs)\n", VOXA_PARAMS$p3_low))
cat(sprintf("  p1_high = %.2f (direct at high LET)\n\n", VOXA_PARAMS$p1_high))

cat("VA Calibration:\n")
for (p in PARTICLE_KEYS) {
  params <- VOXA_PARAMS$particles[[p]]
  cat(sprintf("  %s: δf = %.5f, CV = %.2f%%, LET ≈ %.1f keV/μm\n",
              toupper(p), params$delta_f, params$P_DSB_cv, params$LET_keV_um))
}
# ── Z-indexed δf interpolation (from Step 9 Section 11A) ─────────────────────
# Load Z_delta_f_table from JSON if present; otherwise build from particle_calibration.
# Z_interp assignments: electron=0, proton=1, carbon=6.
Z_interp_map_s12 <- c(electron = 0L, proton = 1L, carbon = 6L,
                       helium = 2L, nitrogen = 7L, neon = 10L, argon = 18L)

if (!is.null(calibration$Z_delta_f_table)) {
  Z_delta_f_tbl <- bind_rows(lapply(calibration$Z_delta_f_table, as_tibble)) %>%
    arrange(Z_interp)
} else {
  # Fallback: build from particle_calibration with hardcoded Z_interp
  Z_delta_f_tbl <- bind_rows(lapply(PARTICLE_KEYS, function(p) {
    pc <- calibration$particle_calibration[[p]]
    tibble(particle = p, Z_interp = Z_interp_map_s12[p],
           delta_f = pc$delta_f, delta_f_ci_low = pc$delta_f_ci_low,
           delta_f_ci_high = pc$delta_f_ci_high)
  })) %>% arrange(Z_interp)
}

#' Log-linear Z-interpolation for δf.  Same logic as Section 11A of Step 9.
interp_delta_f_Z_s12 <- function(Z_target) {
  Z_vals <- Z_delta_f_tbl$Z_interp
  df_vals <- Z_delta_f_tbl$delta_f
  df_lo   <- Z_delta_f_tbl$delta_f_ci_low
  df_hi   <- Z_delta_f_tbl$delta_f_ci_high
  if (Z_target %in% Z_vals) {
    idx <- which(Z_vals == Z_target)
    return(list(delta_f = df_vals[idx], delta_f_ci_low = df_lo[idx],
                delta_f_ci_high = df_hi[idx]))
  }
  # Use Z > 0 for log-space; electron is an anchor at Z_pseudo = 0.5
  pos_mask  <- Z_vals > 0
  Z_pos     <- Z_vals[pos_mask]; df_pos <- df_vals[pos_mask]
  df_lo_pos <- df_lo[pos_mask];  df_hi_pos <- df_hi[pos_mask]
  log_Zt    <- log(max(Z_target, 1e-6))
  log_Zp    <- log(Z_pos)
  if (Z_target < min(Z_pos)) {
    idx_e  <- which(Z_vals == 0); idx_p <- 1
    frac <- (log_Zt - log(0.5)) / (log_Zp[idx_p] - log(0.5))
    frac <- max(0, min(1, frac))
    il   <- function(a, b) exp(log(a) + frac * (log(b) - log(a)))
    return(list(delta_f = il(df_vals[idx_e], df_pos[idx_p]),
                delta_f_ci_low  = il(df_lo[idx_e], df_lo_pos[idx_p]),
                delta_f_ci_high = il(df_hi[idx_e], df_hi_pos[idx_p])))
  }
  idx_lo <- max(which(Z_pos <= Z_target))
  idx_hi <- min(which(Z_pos >= Z_target))
  frac   <- (log_Zt - log_Zp[idx_lo]) / (log_Zp[idx_hi] - log_Zp[idx_lo])
  il     <- function(a, b) exp(log(a) + frac * (log(b) - log(a)))
  list(delta_f = il(df_pos[idx_lo], df_pos[idx_hi]),
       delta_f_ci_low  = il(df_lo_pos[idx_lo], df_lo_pos[idx_hi]),
       delta_f_ci_high = il(df_hi_pos[idx_lo], df_hi_pos[idx_hi]))
}

# Add helium to VOXA_PARAMS using Z-interpolated δf.
# p1/p2/p3 for helium at a given LET would be computed at runtime from OM
# transition functions; here we use the pSOBP LET (10 keV/µm) as default.
# These are updated once helium TOPAS-nBio simulations are available.
he_df_interp <- interp_delta_f_Z_s12(2L)
VOXA_PARAMS$particles[["helium"]] <- list(
  n_dsbs_calibration = 0L,          # No direct calibration yet
  LET_keV_um         = 10.0,        # Helium pSOBP (placeholder; update for dSOBP)
  p1                 = NA_real_,    # Must be computed from OM at runtime for exact LET
  p2                 = NA_real_,
  p3                 = NA_real_,
  delta_f            = he_df_interp$delta_f,
  delta_f_ci_low     = he_df_interp$delta_f_ci_low,
  delta_f_ci_high    = he_df_interp$delta_f_ci_high,
  P_DSB_cv           = NA_real_,    # Not yet validated; use Z-interpolation estimate
  E_mean             = NA_real_,
  E_sd               = NA_real_,
  interpolated       = TRUE
)

cat(sprintf("  HELIUM (Z-interpolated): δf = %.6f [%.6f, %.6f]\n",
            VOXA_PARAMS$particles$helium$delta_f,
            VOXA_PARAMS$particles$helium$delta_f_ci_low,
            VOXA_PARAMS$particles$helium$delta_f_ci_high))
# Section 2: Load Validation Data (400-Dsb)

cat("\n--- SECTION 2: LOADING VALIDATION DATA (400-DSB) ---\n")
val_data_paths <- c(
  "voxa_features_output_validation/all_particles_validation_energy_features.csv",
  "voxa_features_output_validation/all_particles_energy_features.csv",
  "validation_output/all_particles_validation_energy_features.csv"
)

val_data_file <- NULL
for (path in val_data_paths) {
  if (file.exists(path)) {
    val_data_file <- path
    break
  }
}

if (is.null(val_data_file)) {
  stop("ERROR: Validation energy features file not found!\n",
       "Please run: python extract_energy_features_voxa.py --mode validation")
}

cat(sprintf("Loading validation data from: %s\n", val_data_file))
val_data <- read_csv(val_data_file, show_col_types = FALSE)
cat(sprintf("  Total validation DSBs: %d\n\n", nrow(val_data)))

# Summary by particle
cat("Validation dataset summary:\n")
val_summary <- val_data %>%
  group_by(particle) %>%
  summarise(
    n_dsbs = n(),
    E_local_mean = mean(E_local),
    E_local_cv = sd(E_local) / mean(E_local) * 100,
    .groups = "drop"
  )
print(val_summary)
# Store DSB counts
DSB_counts <- setNames(val_summary$n_dsbs, val_summary$particle)


# Section 3: Voxa Model Functions

cat("\n--- SECTION 3: VOxA MODEL FUNCTIONS ---\n")
#' Get particle parameters safely
get_particle_params <- function(particle_name) {
  particle_name <- as.character(particle_name)
  
  if (particle_name %in% names(VOXA_PARAMS$particles)) {
    return(VOXA_PARAMS$particles[[particle_name]])
  }
  
  for (key in names(VOXA_PARAMS$particles)) {
    if (tolower(key) == tolower(particle_name)) {
      return(VOXA_PARAMS$particles[[key]])
    }
  }
  
  stop(sprintf("Particle '%s' not found", particle_name))
}

#' Calculate oxygen fixation probability (Michaelis-Menten)
calc_p_indirect <- function(O2) {
  (O2 + VOXA_PARAMS$K_fix) / (O2 + VOXA_PARAMS$K_fix + VOXA_PARAMS$K_repair)
}

#' Calculate uniform model P_DSB (normalized to 21% O2)
calc_P_DSB_uniform <- function(O2, particle) {
  particle <- as.character(particle)
  params <- get_particle_params(particle)
  
  p1 <- params$p1
  p2 <- params$p2
  p3 <- params$p3
  
  p_ind <- calc_p_indirect(O2)
  p_ind_ref <- calc_p_indirect(VOXA_PARAMS$O2_normoxia)
  
  P_raw <- p1 + p2 * p_ind + p3 * p_ind^2
  P_raw_ref <- p1 + p2 * p_ind_ref + p3 * p_ind_ref^2
  
  P_raw / P_raw_ref
}

#' Calculate OER from P_DSB
calc_OER <- function(O2, particle) {
  1.0 / calc_P_DSB_uniform(O2, particle)
}

#' Compute voxel-aware P_DSB for each DSB
compute_P_DSB_voxel <- function(E_local, particle, O2) {
  particle <- as.character(particle)
  params <- get_particle_params(particle)
  
  E_mean <- mean(E_local)
  E_std <- sd(E_local)
  
  if (E_std < 1e-10) {
    E_zscore <- rep(0, length(E_local))
  } else {
    E_zscore <- (E_local - E_mean) / E_std
  }
  
  f_direct_raw <- params$p1 + params$delta_f * E_zscore
  f_direct <- pmax(pmin(f_direct_raw, VOXA_PARAMS$f_max), VOXA_PARAMS$f_min)
  
  remaining <- 1.0 - f_direct
  total_indirect <- params$p2 + params$p3
  
  if (total_indirect > 1e-10) {
    p2_local <- remaining * (params$p2 / total_indirect)
    p3_local <- remaining * (params$p3 / total_indirect)
  } else {
    p2_local <- rep(0, length(f_direct))
    p3_local <- remaining
  }
  
  p_ind <- calc_p_indirect(O2)
  p_ind_ref <- calc_p_indirect(VOXA_PARAMS$O2_normoxia)
  
  P_raw <- f_direct + p2_local * p_ind + p3_local * p_ind^2
  P_raw_ref <- f_direct + p2_local * p_ind_ref + p3_local * p_ind_ref^2
  
  P_DSB <- ifelse(P_raw_ref > 0, P_raw / P_raw_ref, 1.0)
  pmax(pmin(P_DSB, 1.0), 0.0)
}

cat("Model functions defined.\n\n")


# Section 4: Define Conditions

cat("\n--- SECTION 4: OXYGEN CONDITIONS ---\n")
conditions <- tibble(
  condition = c("Normoxia", "Mild hypoxia", "Moderate hypoxia", "Severe hypoxia", "Anoxia"),
  O2 = c(21.0, 2.1, 0.21, 0.021, 0.001),
  description = c("Air (21%)", "~2%", "~0.2%", "~0.02%", "<0.01%")
)

cat("Oxygen Conditions:\n")
print(conditions)
# Particle information
particle_info <- tibble(
  particle = c("electron", "proton", "carbon"),
  particle_label = c("Electron", "Proton", "Carbon"),
  LET_keV_um = sapply(c("electron", "proton", "carbon"), function(p) {
    get_particle_params(p)$LET_keV_um
  })
)

cat("Particle LET Values (from calibration):\n")
print(particle_info)
# Section 5: Calculate Dsb Retention (Va Model)

cat("\n--- SECTION 5: DSB RETENTION CALCULATIONS (VA MODEL) ---\n")
# Standard order
ordered_particles <- c("electron", "proton", "carbon")

# Calculate retention for each particle and condition
results_list <- list()

for (particle in ordered_particles) {
  cat(sprintf("Processing %s...\n", toupper(particle)))
  
  # Get validation data for this particle
  E_local <- val_data %>% 
    filter(particle == !!particle) %>% 
    pull(E_local)
  
  n_dsbs <- length(E_local)
  params <- get_particle_params(particle)
  
  for (i in 1:nrow(conditions)) {
    O2 <- conditions$O2[i]
    cond_name <- conditions$condition[i]
    
    # Compute P_DSB for each DSB (VA model)
    P_DSB_voxel <- compute_P_DSB_voxel(E_local, particle, O2)
    
    # Uniform model for comparison
    P_DSB_uniform <- calc_P_DSB_uniform(O2, particle)
    
    # Statistics
    P_DSB_mean <- mean(P_DSB_voxel)
    P_DSB_std <- sd(P_DSB_voxel)
    P_DSB_cv <- ifelse(P_DSB_mean > 1e-10, P_DSB_std / P_DSB_mean * 100, 0)
    
    # Expected retained DSBs
    expected_retained <- round(P_DSB_mean * n_dsbs)
    expected_retained_uniform <- round(P_DSB_uniform * n_dsbs)
    
    # Monte Carlo simulation for confidence interval
    n_mc <- 1000
    retained_counts <- numeric(n_mc)
    for (j in 1:n_mc) {
      retained_counts[j] <- sum(runif(n_dsbs) < P_DSB_voxel)
    }
    ci_low <- quantile(retained_counts, 0.025)
    ci_high <- quantile(retained_counts, 0.975)
    
    # OER
    OER <- calc_OER(O2, particle)
    
    results_list[[length(results_list) + 1]] <- tibble(
      particle = particle,
      particle_label = particle_info$particle_label[particle_info$particle == particle],
      LET_keV_um = params$LET_keV_um,
      condition = cond_name,
      O2_pct = O2,
      n_dsbs_total = n_dsbs,
      
      # VA model results
      P_DSB_va_mean = P_DSB_mean,
      P_DSB_va_std = P_DSB_std,
      P_DSB_va_cv = P_DSB_cv,
      retained_va = expected_retained,
      retained_va_ci_low = ci_low,
      retained_va_ci_high = ci_high,
      retained_va_pct = P_DSB_mean * 100,
      
      # Uniform model comparison
      P_DSB_uniform = P_DSB_uniform,
      retained_uniform = expected_retained_uniform,
      retained_uniform_pct = P_DSB_uniform * 100,
      
      # OER
      OER = OER,
      
      # Delta_f used
      delta_f = params$delta_f
    )
  }
}

results <- bind_rows(results_list)

# Set factor levels for ordering
results$condition <- factor(results$condition, levels = conditions$condition)
results$particle <- factor(results$particle, levels = ordered_particles)

# Section 6: Display Results

cat("\n--- SECTION 6: DSB RETENTION TABLES ---\n")
# OER values table
cat("OER values by particle and condition:\n")
cat("┌───────────────────┬──────────┬──────────┬──────────┬──────────┐\n")
cat("│ Condition         │ O₂ (%)   │ Electron │ Proton   │ Carbon   │\n")
cat("├───────────────────┼──────────┼──────────┼──────────┼──────────┤\n")

for (cond in levels(results$condition)) {
  r <- results %>% filter(condition == cond)
  cat(sprintf("│ %-17s │ %8.3f │ %8.2f │ %8.2f │ %8.2f │\n",
              cond, r$O2_pct[1],
              r$OER[r$particle == "electron"],
              r$OER[r$particle == "proton"],
              r$OER[r$particle == "carbon"]))
}
cat("└─────────────���─────┴──────────┴──────────┴──────────┴──────────┘\n\n")

# P_DSB retention percentage table - FIXED header
cat("P_DSB retention (% of normoxia, VA model):\n")
cat("┌───────────────────┬──────────┬─────────────────┬─────────────────┬─────────────────┐\n")
cat(sprintf("│ Condition         │ O₂ (%%)   │ Electron        │ Proton          │ Carbon          │\n"))
cat(sprintf("│                   │          │ (n=%d)          │ (n=%d)          │ (n=%d)          │\n",
            DSB_counts["electron"], DSB_counts["proton"], DSB_counts["carbon"]))
cat("├───────────────────┼──────────┼─────────────────┼─────────────────┼─────────────────┤\n")

for (cond in levels(results$condition)) {
  r <- results %>% filter(condition == cond)
  
  e_val <- r %>% filter(particle == "electron")
  p_val <- r %>% filter(particle == "proton")
  c_val <- r %>% filter(particle == "carbon")
  
  cat(sprintf("│ %-17s │ %8.3f │ %5.1f%% ±%4.1f%%  │ %5.1f%% ±%4.1f%%  │ %5.1f%% ±%4.1f%%  │\n",
              cond, r$O2_pct[1],
              e_val$retained_va_pct, e_val$P_DSB_va_cv,
              p_val$retained_va_pct, p_val$P_DSB_va_cv,
              c_val$retained_va_pct, c_val$P_DSB_va_cv))
}
cat("└───────────────────┴──────────┴─────────────────┴─────────────────┴─────────────────┘\n\n")

# Expected retained DSB counts - FIXED
cat("Expected retained DSBs (VA model with 95% CI):\n")
cat("┌───────────────────┬──────────┬─────────────────┬─────────────────┬─────────────────┐\n")
cat("│ Condition         │ O₂ (%)   │ Electron        │ Proton          │ Carbon          │\n")
cat("├───────────────────┼──────────┼─────────────────┼─────────────────┼─────────────────┤\n")

for (cond in levels(results$condition)) {
  r <- results %>% filter(condition == cond)
  
  e_val <- r %>% filter(particle == "electron")
  p_val <- r %>% filter(particle == "proton")
  c_val <- r %>% filter(particle == "carbon")
  
  # Use %.0f instead of %d, or convert to integer
  cat(sprintf("│ %-17s │ %8.3f │ %3.0f [%3.0f, %3.0f]  │ %3.0f [%3.0f, %3.0f]  │ %3.0f [%3.0f, %3.0f]  │\n",
              cond, r$O2_pct[1],
              e_val$retained_va, e_val$retained_va_ci_low, e_val$retained_va_ci_high,
              p_val$retained_va, p_val$retained_va_ci_low, p_val$retained_va_ci_high,
              c_val$retained_va, c_val$retained_va_ci_low, c_val$retained_va_ci_high))
}
cat("└───────────────────┴──────────┴─────────────────┴─────────────────┴─────────────────┘\n\n")


# Section 7: Case Fractions

cat("\n--- SECTION 7: CASE FRACTIONS BY PARTICLE ---\n")
for (particle in ordered_particles) {
  params <- get_particle_params(particle)
  label <- particle_info$particle_label[particle_info$particle == particle]
  
  cat(sprintf("%s (LET ≈ %.1f keV/μm):\n", label, params$LET_keV_um))
  cat(sprintf("  p1 (direct)   = %.4f (%.1f%%)\n", params$p1, 100*params$p1))
  cat(sprintf("  p2 (hybrid)   = %.4f (%.1f%%)\n", params$p2, 100*params$p2))
  cat(sprintf("  p3 (indirect) = %.4f (%.1f%%)\n", params$p3, 100*params$p3))
  cat(sprintf("  δf            = %.5f\n\n", params$delta_f))
}


# Section 8: Oxygen Fixation Probability

cat("\n--- SECTION 8: OXYGEN FIXATION PROBABILITY (p_ind) ---\n")
p_ind_table <- conditions %>%
  mutate(
    p_ind = sapply(O2, calc_p_indirect),
    p_ind_pct = 100 * p_ind
  )

cat("p_ind = (O₂ + K_fix) / (O₂ + K_fix + K_repair)\n")
cat(sprintf("K_fix = %.4f%%, K_repair = %.4f%%\n\n", VOXA_PARAMS$K_fix, VOXA_PARAMS$K_repair))

cat("┌─��─────────────────┬──────────┬──────────┬──────────┐\n")
cat("│ Condition         │ O₂ (%)   │ p_ind    │ p_ind(%) │\n")
cat("├───────────────────┼──────────┼──────────┼──────────┤\n")
for (i in 1:nrow(p_ind_table)) {
  r <- p_ind_table[i,]
  cat(sprintf("│ %-17s │ %8.3f │ %8.4f │ %7.1f%% │\n",
              r$condition, r$O2, r$p_ind, r$p_ind_pct))
}
cat("└───────────────────┴──────────┴──────────┴──────────┘\n\n")


# Section 9: Generate Figure 47

cat("\n--- SECTION 9: GENERATING FIGURE 47 ---\n")
# Create local color mapping
local_particle_colors <- setNames(
  c(mediterranean_blue, lemon_yellow, pompeii_red),
  ordered_particles
)

# Panel A: Retained DSBs by condition
plot_data_a <- results %>%
  mutate(particle_label = factor(particle_label, 
                                 levels = c("Electron", "Proton", "Carbon")))

panel_a <- ggplot(plot_data_a, aes(x = condition, y = retained_va, fill = particle)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = retained_va_ci_low, ymax = retained_va_ci_high),
                position = position_dodge(width = 0.8), width = 0.2, linewidth = 0.5) +
  scale_fill_manual(values = local_particle_colors,
                    labels = c("Electron", "Proton", "Carbon"),
                    name = "Particle") +
  labs(
    title = "A. Expected Retained DSBs by Oxygen Condition",
    subtitle = "VOxA VA model predictions for 400-DSB validation datasets",
    x = "Oxygen Condition",
    y = "Expected Retained DSBs"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = warm_stone),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# Panel B: P_DSB retention percentage
panel_b <- ggplot(plot_data_a, aes(x = condition, y = retained_va_pct, 
                                   color = particle, group = particle)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = retained_va_ci_low / n_dsbs_total * 100,
                  ymax = retained_va_ci_high / n_dsbs_total * 100,
                  fill = particle), alpha = 0.2, color = NA) +
  scale_color_manual(values = local_particle_colors,
                     labels = c("Electron", "Proton", "Carbon"),
                     name = "Particle") +
  scale_fill_manual(values = local_particle_colors, guide = "none") +
  labs(
    title = "B. P_DSB Retention (% of Normoxia)",
    subtitle = "Shaded region = 95% CI from Monte Carlo",
    x = "Oxygen Condition",
    y = "Retention (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = warm_stone),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 105))

# Panel C: OER comparison
panel_c <- ggplot(plot_data_a %>% filter(condition != "Normoxia"), 
                  aes(x = condition, y = OER, fill = particle)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", OER)),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  scale_fill_manual(values = local_particle_colors,
                    labels = c("Electron", "Proton", "Carbon"),
                    name = "Particle") +
  labs(
    title = "C. OER by Particle and Condition",
    subtitle = "OER = 1 / P_DSB (higher = more radiosensitization needed)",
    x = "Oxygen Condition",
    y = "OER"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = warm_stone),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0, max(results$OER) * 1.15))

# Combine panels
combined <- grid.arrange(
  panel_a, panel_b, panel_c,
  ncol = 3,
  top = textGrob("VOxA Model DSB Retention Predictions (400-DSB Validation Datasets)",
                 gp = gpar(fontface = "bold", fontsize = 14))
)

ggsave("figures/fig47_dsb_retention.png", combined,
       width = 15, height = 6, dpi = 300, bg = "white")
# Section 10: Latex Table Output

cat("\n--- SECTION 10: LATEX TABLE OUTPUT ---\n")
cat("% LaTeX table for manuscript - DSB Retention (VOxA VA Model)\n")
cat("\\begin{table}[H]\n")
cat("\\centering\n")
cat(sprintf("\\caption{DSB retention predictions by VOxA VA model ($K_{\\text{fix}}$ = %.4f\\%%, $K_{\\text{repair}}$ = %.4f\\%%, $R^2$ = %.3f). Values show expected retained DSBs [95\\%% CI] from 400-DSB validation datasets.}\n",
            VOXA_PARAMS$K_fix, VOXA_PARAMS$K_repair, VOXA_PARAMS$R_squared))
cat("\\label{tab:voxa_dsb_retention}\n")
cat("\\begin{tabular}{@{}lcccc@{}}\n")
cat("\\toprule\n")
cat("\\textbf{Condition} & \\textbf{O$_2$ (\\%)} & \\textbf{Electron} & \\textbf{Proton} & \\textbf{Carbon} \\\\\n")
cat(" & & (%.1f keV/$\\mu$m) & (%.1f keV/$\\mu$m) & (%.1f keV/$\\mu$m) \\\\\n",
    VOXA_PARAMS$particles$electron$LET_keV_um,
    VOXA_PARAMS$particles$proton$LET_keV_um,
    VOXA_PARAMS$particles$carbon$LET_keV_um)
cat("\\midrule\n")

for (cond in levels(results$condition)) {
  r <- results %>% filter(condition == cond)
  
  e_val <- r %>% filter(particle == "electron")
  p_val <- r %>% filter(particle == "proton")
  c_val <- r %>% filter(particle == "carbon")
  
  cat(sprintf("%s & %.3f & %.0f [%.0f, %.0f] & %.0f [%.0f, %.0f] & %.0f [%.0f, %.0f] \\\\\n",
              cond, r$O2_pct[1],
              e_val$retained_va, e_val$retained_va_ci_low, e_val$retained_va_ci_high,
              p_val$retained_va, p_val$retained_va_ci_low, p_val$retained_va_ci_high,
              c_val$retained_va, c_val$retained_va_ci_low, c_val$retained_va_ci_high))
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n\n")

# Compact table with percentages
cat("% Compact LaTeX table (percentages)\n")
cat("\\begin{table}[H]\n")
cat("\\centering\n")
cat("\\caption{DSB retention percentage by VOxA model (relative to normoxia)}\n")
cat("\\label{tab:voxa_retention_pct}\n")
cat("\\begin{tabular}{@{}lccccc@{}}\n")
cat("\\toprule\n")
cat("\\textbf{Condition} & \\textbf{O$_2$ (\\%)} & \\textbf{$p_{\\text{ind}}$} & \\textbf{Electron} & \\textbf{Proton} & \\textbf{Carbon} \\\\\n")
cat("\\midrule\n")

for (cond in levels(results$condition)) {
  r <- results %>% filter(condition == cond)
  p_ind <- calc_p_indirect(r$O2_pct[1])
  
  cat(sprintf("%s & %.3f & %.3f & %.1f\\%% & %.1f\\%% & %.1f\\%% \\\\\n",
              cond, r$O2_pct[1], p_ind,
              r$retained_va_pct[r$particle == "electron"],
              r$retained_va_pct[r$particle == "proton"],
              r$retained_va_pct[r$particle == "carbon"]))
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n\n")


# Section 11: Save Results

cat("\n--- SECTION 11: SAVING RESULTS ---\n")
# Detailed results
write_csv(results, "results/voxa_dsb_retention_detailed.csv")
# Summary table
results_summary <- results %>%
  select(particle, particle_label, LET_keV_um, condition, O2_pct,
         P_DSB_va_mean, P_DSB_va_cv, retained_va, retained_va_pct, OER)
write_csv(results_summary, "results/voxa_dsb_retention_table.csv")
# p_ind table
write_csv(p_ind_table, "results/voxa_p_ind_table.csv")
# JSON export
json_export <- list(
  model_version = VOXA_PARAMS$model_version,
  generation_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  parameters = list(
    K_fix = VOXA_PARAMS$K_fix,
    K_repair = VOXA_PARAMS$K_repair,
    R_squared = VOXA_PARAMS$R_squared,
    OER_max_retention = VOXA_PARAMS$OER_max_retention
  ),
  validation_datasets = list(
    source = val_data_file,
    electron = list(n_dsbs = DSB_counts["electron"]),
    proton = list(n_dsbs = DSB_counts["proton"]),
    carbon = list(n_dsbs = DSB_counts["carbon"])
  ),
  results = split(results, results$condition) %>%
    lapply(function(df) {
      split(df, df$particle) %>%
        lapply(function(x) as.list(x[1,]))
    })
)

write_json(json_export, "results/voxa_dsb_retention_table.json", 
           pretty = TRUE, auto_unbox = TRUE)
# Final Summary

p_ind_anoxia <- calc_p_indirect(0.001)
p_ind_normoxia <- calc_p_indirect(21.0)

cat(sprintf("║   Model: VOxA v1.0                                                  ║\n"))
cat(sprintf("║   Base R² = %.4f                                                    ║\n", VOXA_PARAMS$R_squared))
cat(sprintf("║   K_fix = %.4f%%  |  K_repair = %.4f%%                             ║\n", 
            VOXA_PARAMS$K_fix, VOXA_PARAMS$K_repair))
cat(sprintf("║   OER_max = %.2f (retention) | %.2f (survival)                      ║\n",
            VOXA_PARAMS$OER_max_retention, VOXA_PARAMS$OER_max_survival))
cat(sprintf("║   p_ind(anoxia)   = %.4f                                           ║\n", p_ind_anoxia))
cat(sprintf("║   p_ind(normoxia) = %.4f                                           ║\n", p_ind_normoxia))
for (p in ordered_particles) {
  params <- get_particle_params(p)
  n <- DSB_counts[p]
  anoxia_ret <- results %>% 
    filter(particle == p, condition == "Anoxia") %>% 
    pull(retained_va)
  cat(sprintf("║     %s: n=%d, anoxia retained=%d (%.1f%%)                   ║\n",
              toupper(substr(p, 1, 1)), n, anoxia_ret, anoxia_ret/n*100))
}
cat("\nStep 12 complete.\n\n")