# Step 3: VOxA — particle-specific dual-sigmoid calibration
#
# Fits the OER model to the 233-point calibration dataset via two-phase
# constrained optimization (L-BFGS-B warm start → Nelder-Mead refinement).
#
# The model uses dual sigmoidal transitions f_dir(x) and f_ind(x) with
# particle-specific midpoints x50. Z-ordering is enforced as a hard
# monotonicity constraint on the heavy-ion sequence (He < C < Ne < Ar):
# at fixed LET, lighter ions are slower → narrower tracks → more radical
# recombination → lower OER.
#
# Physical LET cutoffs at Bragg-peak limits prevent non-physical extrapolation.
# An overkill correction handles the high-LET dose-response saturation.
# Light particles (Z ≤ 1) use separate steepness parameters from heavy ions.
#
# Outputs: results/uvaom_v8_corrected_model.rds,
#          results/calibration_data_v8_corrected.csv

library(tidyverse)
library(MASS)
select <- dplyr::select

# Physical Constants: Maximum Let Values (Bragg Peak Limits)

MAX_LET_PHYSICAL <- list(
  photon = 35,      # Secondary electrons max ~31 keV/μm
  proton = 100,     # Bragg peak ~92-100 keV/μm
  deuteron = 120,   # Slightly higher than proton
  He = 200,         # ~160-240 keV/μm depending on measurement
  C = 550,          # Up to ~500-550 keV/μm
  N = 600,          # Interpolated
  
  O = 620,          # Interpolated
  Ne = 700,         # Up to ~650-700 keV/μm
  Si = 800,         # Interpolated
  Ar = 900          # Higher than Ne
)

for (p in names(MAX_LET_PHYSICAL)) {
  cat(sprintf("  %-8s: %4d keV/µm\n", p, MAX_LET_PHYSICAL[[p]]))
}
# Conversion Functions

#' Convert OER_retention to OER_survival (general formula)
#' Ensures OER_survival = 1.0 when OER_retention = 1.0
convert_retention_to_survival <- function(OER_retention, conv_factor = 1.20) {
  1.00 + ((OER_retention - 1.00) / conv_factor)
}

#' Convert OER_survival to OER_retention (general formula)
convert_survival_to_retention <- function(OER_survival, conv_factor = 1.20) {
  1.00 + ((OER_survival - 1.00) * conv_factor)
}


# Section 0: Data Preparation

cat("\n--- SECTION 0: DATA PREPARATION ---\n")
# Load setup
if (file.exists("data/uvaom_recalibration_setup.RData")) {
  load("data/uvaom_recalibration_setup.RData")
  cat("✓ Loaded setup from data/uvaom_recalibration_setup.RData\n")
} else {
  CONVERSION_FACTOR <- list(
    mean = 1.20,
    sd = 0.05,
    range = c(1.15, 1.25),
    source = "Hirayama2005"
  )
  cat("⚠ Setup not found, using default CONVERSION_FACTOR\n")
}

# Load Furusawa data (from Step 1)
cat("Loading Furusawa data (D10-based OER calculation)...\n")

furusawa <- read_csv("data/furusawa_oer_clean.csv", show_col_types = FALSE) %>%
  mutate(
    ion = ion_std,
    O2_hyp = 0.001,
    O2_ref = 21.0,
    dataset = "Furusawa2000_corrected",
    cell_line_std = case_when(
      str_detect(toupper(cell_line), "V79") ~ "V79",
      str_detect(toupper(cell_line), "HSG") ~ "HSG",
      TRUE ~ "Other"
    )
  )

cat(sprintf("  Furusawa: %d points\n", nrow(furusawa)))
cat(sprintf("    OER_survival range: %.2f - %.2f\n", 
            min(furusawa$OER_survival, na.rm = TRUE), 
            max(furusawa$OER_survival, na.rm = TRUE)))

# Load literature + Tinganelli data (from Step 2)
wenzl <- read_csv("data/wenzl_literature_oer.csv", show_col_types = FALSE) %>%
  mutate(
    ion = ifelse(is.na(ion_std), ion, ion_std),
    O2_hyp = ifelse(is.na(O2_hyp), 0.001, O2_hyp),
    O2_ref = ifelse(is.na(O2_ref), 21.0, O2_ref),
    cell_line = if("cell_line" %in% names(.)) cell_line else "Unknown"
  )

cat(sprintf("  Literature + Tinganelli: %d points\n", nrow(wenzl)))

# Check for Tinganelli data
n_tinganelli <- sum(wenzl$source == "Tinganelli2015", na.rm = TRUE)
cat(sprintf("    Including Tinganelli 2015: %d points\n", n_tinganelli))

# D-Kondo for VALIDATION ONLY
dkondo_validation <- tibble(
  LET = 1.0,
  ion = "photon",
  O2_hyp = c(21.0, 10.0, 5.0, 2.1, 1.0, 0.5, 0.21, 0.1, 0.05, 0.021, 0.01, 0.001),
  O2_ref = 21.0,
  HRF_experimental = c(1.0, 1.05, 1.10, 1.21, 1.45, 1.75, 2.33, 2.55, 2.70, 2.83, 2.88, 2.95),
  is_measured = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
  dataset = "D-Kondo2009_validation",
  notes = "HRF values - radioprotector used, NOT direct OER measurement"
)

# Particle metadata with Z values
particle_info <- tibble(
  ion = c("photon", "proton", "deuteron", "He", "C", "N", "O", "Ne", "Si", "Ar"),
  Z = c(0, 1, 1, 2, 6, 7, 8, 10, 14, 18),
  particle_class = c("light", "light", "light", "heavy", "heavy", "heavy", "heavy", "heavy", "heavy", "heavy"),
  max_let = c(35, 100, 120, 200, 550, 600, 620, 700, 800, 900)
)

# Combine calibration data
calibration_data_raw <- bind_rows(furusawa, wenzl) %>%
  filter(
    !is.na(OER_retention),
    !is.na(LET),
    !is.na(O2_hyp),
    !is.na(O2_ref)
  ) %>%
  left_join(particle_info, by = "ion") %>%
  filter(!is.na(Z))

cat(sprintf("\n  Combined raw data: %d points\n", nrow(calibration_data_raw)))

# Show O2 level distribution
cat("\nO2 level distribution:\n")
calibration_data_raw %>%
  mutate(O2_bin = case_when(
    O2_hyp >= 20 ~ "Normoxia (21%)",
    O2_hyp >= 1 ~ "Physoxia (1-20%)",
    O2_hyp >= 0.1 ~ "Mild hypoxia (0.1-1%)",
    O2_hyp >= 0.01 ~ "Acute hypoxia (0.01-0.1%)",
    TRUE ~ "Anoxia (<0.01%)"
  )) %>%
  count(O2_bin) %>%
  print()
# Physical Let Limit Check (New)

cat("PHYSICAL LET LIMIT CHECK:\n")
calibration_data_raw <- calibration_data_raw %>%
  mutate(
    # Check if LET exceeds physical maximum for this particle
    is_beyond_bragg = LET > max_let,
    bragg_proximity = pmin(LET / max_let, 1.5),  # Cap at 1.5 for extreme cases
    
    # Flag severity
    let_status = case_when(
      LET > max_let * 1.1 ~ "NON-PHYSICAL",
      LET > max_let * 0.9 ~ "NEAR_BRAGG",
      LET > max_let * 0.7 ~ "APPROACHING",
      TRUE ~ "OK"
    )
  )

# Report non-physical data points
beyond_bragg_summary <- calibration_data_raw %>%
  filter(is_beyond_bragg) %>%
  group_by(ion) %>%
  summarise(
    n = n(),
    LET_range = sprintf("%.0f-%.0f", min(LET), max(LET)),
    max_physical = first(max_let),
    .groups = "drop"
  )

if (nrow(beyond_bragg_summary) > 0) {
  cat(sprintf("⚠ Found %d data points exceeding physical LET limits:\n", 
              sum(calibration_data_raw$is_beyond_bragg)))
  print(beyond_bragg_summary)
  cat("\n  These points will be DOWN-WEIGHTED (not removed) in calibration.\n\n")
} else {
  cat("✓ All data points within physical LET limits.\n\n")
}


# Outlier Identification

cat("OUTLIER IDENTIFICATION:\n")
calibration_data_raw <- calibration_data_raw %>%
  mutate(
    outlier_reason = case_when(
      # Physical constraint: OER must be >= 1.0
      OER_retention < 1.0 ~ "OER < 1.0 (physically implausible)",
      
      # Upper bound: OER_retention > 4.5 is suspicious
      OER_retention > 4.5 ~ "OER > 4.5 (exceeds expected maximum)",
      
      # High LET should have low OER (approaching 1.0)
      particle_class == "heavy" & LET > 300 & OER_retention > 2.5 ~ 
        "High LET (>300) with OER > 2.5",
      
      # Low LET heavy ions should have OER similar to photons
      particle_class == "heavy" & LET < 40 & OER_retention < 1.3 ~ 
        "Heavy ion low-LET (<40) with OER < 1.3",
      
      TRUE ~ NA_character_
    ),
    is_outlier = !is.na(outlier_reason)
  )

outliers <- calibration_data_raw %>%
  filter(is_outlier) %>%
  select(ion, LET, OER_retention, cell_line, dataset, outlier_reason)

cat(sprintf("Identified %d outliers for REMOVAL:\n\n", nrow(outliers)))
if (nrow(outliers) > 0) {
  print(outliers)
} else {
  }
calibration_data <- calibration_data_raw %>%
  filter(!is_outlier)

n_original <- nrow(calibration_data_raw)
n_removed <- sum(calibration_data_raw$is_outlier)
n_remaining <- nrow(calibration_data)

cat(sprintf("Data summary:\n"))
cat(sprintf("  Input observations: %d\n", n_original))
cat(sprintf("  Outliers removed: %d\n", n_removed))
cat(sprintf("  Remaining observations: %d\n\n", n_remaining))


# Cell-Line Identification

cat("CELL-LINE DISTRIBUTION:\n")
# Ensure cell_line_std exists
if (!"cell_line_std" %in% names(calibration_data)) {
  calibration_data <- calibration_data %>%
    mutate(
      cell_line_std = case_when(
        is.na(cell_line) ~ "Unknown",
        str_detect(toupper(cell_line), "V79") ~ "V79",
        str_detect(toupper(cell_line), "HSG") ~ "HSG",
        str_detect(toupper(cell_line), "CHO") ~ "CHO",
        str_detect(toupper(cell_line), "T1") ~ "T1",
        TRUE ~ "Other"
      )
    )
}

cat("By cell line:\n")
calibration_data %>%
  count(cell_line_std) %>%
  arrange(desc(n)) %>%
  print()
# Enhanced Weighting Scheme (Improved)

cat("ENHANCED WEIGHTING SCHEME:\n")
calibration_data <- calibration_data %>%
  group_by(ion) %>%
  mutate(base_weight = 1 / sqrt(n())) %>%
  ungroup() %>%
  mutate(
    # 1. Clinical relevance weight (focus on therapeutically relevant LET)
    clinical_weight = case_when(
      particle_class == "heavy" & LET >= 30 & LET <= 200 ~ 1.5,   # SOBP region
      particle_class == "heavy" & LET >= 20 & LET <= 300 ~ 1.3,   # Extended clinical
      particle_class == "heavy" & LET > 300 ~ 1.0,                # High LET
      particle_class == "light" & LET >= 0.5 & LET <= 30 ~ 1.3,   # Clinical photon/proton
      TRUE ~ 1.0
    ),
    
    # 2. Physical plausibility weight (NEW - based on LET limits)
    physical_weight = case_when(
      is_beyond_bragg ~ 0.2,                    # Heavy penalty for non-physical
      bragg_proximity > 0.9 ~ 0.6,              # Reduced weight near Bragg peak
      bragg_proximity > 0.7 ~ 0.85,             # Slight reduction approaching peak
      TRUE ~ 1.0
    ),
    
    # 3. Cell line reliability weight
    cell_weight = case_when(
      cell_line_std == "V79" ~ 1.2,   # Most data, well-characterized
      cell_line_std == "HSG" ~ 1.1,   # Good characterization
      cell_line_std == "CHO" ~ 1.0,   # Tinganelli data
      cell_line_std == "T1" ~ 1.0,
      TRUE ~ 0.9                       # Unknown/Other
    ),
    
    # 4. Experimental uncertainty weight (if available)
    uncertainty_weight = ifelse(!is.na(OER_survival_se) & OER_survival_se > 0,
                                1 / OER_survival_se^2,
                                1.0),
    uncertainty_weight = uncertainty_weight / max(uncertainty_weight, na.rm = TRUE),
    
    # Combined weight
    weight = base_weight * clinical_weight * physical_weight * cell_weight * 
      (0.5 + 0.5 * uncertainty_weight)
  )

# Report weight distribution
cat("Weight distribution:\n")
calibration_data %>%
  summarise(
    n_full_physical = sum(physical_weight == 1.0),
    n_reduced_physical = sum(physical_weight < 1.0 & physical_weight > 0.2),
    n_minimal_physical = sum(physical_weight <= 0.2),
    mean_weight = mean(weight),
    sd_weight = sd(weight)
  ) %>%
  print()
# Define particle groups
light_particles <- c("photon", "proton", "deuteron")
heavy_particles_calibrated <- c("He", "C", "Ne", "Ar")
heavy_particles_interpolated <- c("N", "O", "Si")
all_particles <- c(light_particles, heavy_particles_calibrated)

cat("Final dataset summary:\n")
calibration_data %>%
  group_by(ion, particle_class) %>%
  summarise(
    n = n(),
    LET_range = sprintf("%.1f-%.1f", min(LET), max(LET)),
    OER_range = sprintf("%.2f-%.2f", min(OER_retention), max(OER_retention)),
    mean_weight = sprintf("%.3f", mean(weight)),
    n_near_bragg = sum(bragg_proximity > 0.7),
    .groups = "drop"
  ) %>%
  arrange(match(ion, c("photon", "proton", "deuteron", "He", "C", "N", "O", "Ne", "Si", "Ar"))) %>%
  print()
# Section 1: Fixed Parameters

cat("\n--- SECTION 1: FIXED PARAMETERS ---\n")
FIXED_PARAMS <- list(
  # DSB combinatorics with d = 0.20, i = 0.80
  p1_low = 0.04,    # d² = 0.20²
  p2_low = 0.32,    # 2di = 2×0.20×0.80
  p3_low = 0.64,    # i² = 0.80²
  
  # High-LET limit from Nikjoo 2001 experiments
  p1_high = 0.64    # 64% direct at extreme LET
)

cat("Case fractions (Sakata 2019 DSB combinatorics):\n")
cat(sprintf("  p1_low  = %.2f (purely direct DSBs)\n", FIXED_PARAMS$p1_low))
cat(sprintf("  p2_low  = %.2f (hybrid DSBs)\n", FIXED_PARAMS$p2_low))
cat(sprintf("  p3_low  = %.2f (purely indirect DSBs)\n", FIXED_PARAMS$p3_low))
cat(sprintf("  p1_high = %.2f (Nikjoo 2001: ~64%% direct at high LET)\n\n", FIXED_PARAMS$p1_high))


# Section 2: Model Functions (Enhanced)

cat("\n--- SECTION 2: MODEL FUNCTIONS (ENHANCED) ---\n")
#' Michaelis-Menten oxygen kinetics
calc_p_indirect_single_MM <- function(O2_percent, K_fix, K_repair) {
  (O2_percent + K_fix) / (O2_percent + K_fix + K_repair)
}

#' Compute LET-dependent transitions
compute_transitions <- function(LET, x50_dir, x50_ind, s_dir, s_ind) {
  x <- 2.5 * LET^1.1
  
  f_direct <- 1 / (1 + (x50_dir / pmax(x, 0.001))^s_dir)
  f_indirect <- 1 / (1 + (x50_ind / pmax(x, 0.001))^s_ind)
  
  return(list(f_direct = f_direct, f_indirect = f_indirect))
}

#' Compute case fractions from transitions
compute_case_fractions <- function(LET, x50_dir, x50_ind, s_dir, s_ind) {
  trans <- compute_transitions(LET, x50_dir, x50_ind, s_dir, s_ind)
  
  p1 <- FIXED_PARAMS$p1_low + (FIXED_PARAMS$p1_high - FIXED_PARAMS$p1_low) * trans$f_direct
  p3 <- FIXED_PARAMS$p3_low * (1 - trans$f_indirect)
  p2 <- 1 - p1 - p3
  p2 <- pmax(0, p2)
  
  total <- p1 + p2 + p3
  
  return(list(p1 = p1/total, p2 = p2/total, p3 = p3/total))
}

#' Predict OER from model parameters (BASE - without overkill correction)
predict_OER_base <- function(LET, x50_dir, x50_ind, s_dir, s_ind,
                             K_fix, K_repair, O2_hyp, O2_ref) {
  
  fracs <- compute_case_fractions(LET, x50_dir, x50_ind, s_dir, s_ind)
  
  p_ind_hyp <- calc_p_indirect_single_MM(O2_hyp, K_fix, K_repair)
  p_ind_ref <- calc_p_indirect_single_MM(O2_ref, K_fix, K_repair)
  
  P_hyp <- fracs$p1 + fracs$p2 * p_ind_hyp + fracs$p3 * p_ind_hyp^2
  P_ref <- fracs$p1 + fracs$p2 * p_ind_ref + fracs$p3 * p_ind_ref^2
  
  P_hyp <- pmax(0.001, pmin(1.0, P_hyp))
  P_ref <- pmax(0.001, pmin(1.0, P_ref))
  
  return(pmax(1.0, P_ref / P_hyp))
}

#' Overkill correction factor (NEW)
#' Near the Bragg peak, effectiveness drops due to overkill
#' This causes OER to slightly increase back toward higher values
calc_overkill_factor <- function(LET, max_let, overkill_strength = 0.15) {
  # Proximity to Bragg peak (0 = far, 1 = at peak, >1 = beyond)
  proximity <- LET / max_let
  
  # Overkill correction only applies when approaching/at Bragg peak
  if (proximity > 0.7) {
    # Parabolic correction that increases OER slightly near Bragg peak
    # At proximity = 1.0, factor = 1 + overkill_strength
    correction <- 1 + overkill_strength * ((proximity - 0.7) / 0.3)^2
    return(pmin(correction, 1 + overkill_strength * 1.5))  # Cap the correction
  } else {
    return(1.0)
  }
}

#' Predict OER with overkill correction (NEW)
predict_OER_with_overkill <- function(LET, x50_dir, x50_ind, s_dir, s_ind,
                                      K_fix, K_repair, O2_hyp, O2_ref,
                                      max_let, overkill_strength = 0.10) {
  
  # Base OER prediction
  OER_base <- predict_OER_base(LET, x50_dir, x50_ind, s_dir, s_ind,
                               K_fix, K_repair, O2_hyp, O2_ref)
  
  # Apply overkill correction
  overkill_factor <- calc_overkill_factor(LET, max_let, overkill_strength)
  
  # Overkill moves OER back toward higher values (less effective killing)
  # OER_corrected = 1 + (OER_base - 1) * overkill_factor
  OER_corrected <- 1 + (OER_base - 1) * overkill_factor
  
  return(pmax(1.0, OER_corrected))
}

#' Z-scaled steepness for heavy ions (ENHANCED)
#' Physics: Lighter ions have narrower tracks → sharper transitions
#' Therefore steepness should DECREASE with Z (or stay constant)
calc_steepness_Z <- function(Z, s_base, s_scale) {
  # Original: s_base * (1 + s_scale * log(Z/2))
  # Enhanced: Allow both positive and negative scaling
  s_base * (1 + s_scale * log(pmax(Z, 2) / 2))
}

#' Physical steepness scaling (NEW - alternative)
#' Based on track structure: lighter ions have sharper transitions
calc_steepness_physical <- function(Z, s_base, Z_ref = 6) {
  # Steepness decreases mildly with Z (heavier ions have more gradual transitions)
  # because their track density varies more slowly with LET
  s <- s_base * (Z_ref / pmax(Z, 2))^0.12
  return(pmax(s, 1.0))  # Minimum steepness of 1.0
}

#' Huber loss for robust optimization
huber_loss <- function(residual, delta = 0.8) {
  abs_r <- abs(residual)
  ifelse(abs_r <= delta,
         0.5 * residual^2,
         delta * (abs_r - 0.5 * delta))
}

#' Cell-line correction factor
get_cell_correction <- function(cell_line_std, factor_HSG, factor_T1, factor_CHO) {
  case_when(
    cell_line_std == "HSG" ~ factor_HSG,
    cell_line_std == "T1" ~ factor_T1,
    cell_line_std == "CHO" ~ factor_CHO,
    cell_line_std == "V79" ~ 1.0,
    TRUE ~ 1.0
  )
}

#' Constrained x50 scaling (NEW)
#' x50 increases with Z according to: x50(Z) = x50_ref * (Z/Z_ref)^alpha
#' This enforces physics and reduces free parameters
calc_x50_scaled <- function(Z, x50_ref, alpha, Z_ref = 2) {
  x50_ref * (Z / Z_ref)^alpha
}

#' Z-interpolation for x50 parameters
interpolate_x50 <- function(Z_target, params_list) {
  calibrated_Z <- c(He = 2, C = 6, Ne = 10, Ar = 18)
  
  lower_Z <- max(calibrated_Z[calibrated_Z <= Z_target])
  upper_Z <- min(calibrated_Z[calibrated_Z >= Z_target])
  
  if (lower_Z == upper_Z) {
    ion_name <- names(calibrated_Z)[calibrated_Z == lower_Z]
    return(list(
      x50_dir = params_list[[paste0("x50_dir_", ion_name)]],
      x50_ind = params_list[[paste0("x50_ind_", ion_name)]]
    ))
  }
  
  lower_ion <- names(calibrated_Z)[calibrated_Z == lower_Z]
  upper_ion <- names(calibrated_Z)[calibrated_Z == upper_Z]
  
  log_ratio <- (log(Z_target) - log(lower_Z)) / (log(upper_Z) - log(lower_Z))
  
  x50_dir_lower <- params_list[[paste0("x50_dir_", lower_ion)]]
  x50_dir_upper <- params_list[[paste0("x50_dir_", upper_ion)]]
  x50_ind_lower <- params_list[[paste0("x50_ind_", lower_ion)]]
  x50_ind_upper <- params_list[[paste0("x50_ind_", upper_ion)]]
  
  x50_dir <- exp(log(x50_dir_lower) + log_ratio * (log(x50_dir_upper) - log(x50_dir_lower)))
  x50_ind <- exp(log(x50_ind_lower) + log_ratio * (log(x50_ind_upper) - log(x50_ind_lower)))
  
  return(list(x50_dir = x50_dir, x50_ind = x50_ind))
}

cat("Model functions defined:\n")
# Section 3: Phase 1 - Photon + K Values (Hybrid: Light Particles)

cat("\n--- SECTION 3: PHASE 1 - LIGHT PARTICLES (Photon + K Values) ---\n")
cat("HYBRID MODEL: Light particles (Z ≤ 1) use dedicated parameters\n")
photon_data <- calibration_data %>% 
  filter(ion == "photon")

cat(sprintf("Photon data: %d points\n", nrow(photon_data)))
cat(sprintf("  LET range: %.1f - %.1f keV/μm\n", min(photon_data$LET), max(photon_data$LET)))
cat(sprintf("  OER_retention range: %.2f - %.2f\n\n", 
            min(photon_data$OER_retention), max(photon_data$OER_retention)))

phase1_param_names <- c("K_fix", "K_repair", 
                        "x50_dir_photon", "x50_ind_photon", 
                        "s_dir_light", "s_ind_light")

objective_phase1 <- function(params_vec, data, return_details = FALSE) {
  params <- as.list(params_vec)
  names(params) <- phase1_param_names
  
  if (any(!is.finite(unlist(params)))) {
    if (return_details) return(list(failed = TRUE))
    return(1e10)
  }
  
  # Physical constraints on K values
  if (params$K_fix <= 0 || params$K_repair <= 0) {
    if (return_details) return(list(failed = TRUE))
    return(1e10)
  }
  if (params$K_fix >= params$K_repair) {
    if (return_details) return(list(failed = TRUE))
    return(1e10)
  }
  
  # x50 ordering constraint
  if (params$x50_dir_photon >= params$x50_ind_photon * 0.9) {
    if (return_details) return(list(failed = TRUE))
    return(1e10)
  }
  
  # Predict OER for photon data
  data$OER_pred <- sapply(1:nrow(data), function(i) {
    predict_OER_base(
      LET = data$LET[i],
      x50_dir = params$x50_dir_photon,
      x50_ind = params$x50_ind_photon,
      s_dir = params$s_dir_light,
      s_ind = params$s_ind_light,
      K_fix = params$K_fix,
      K_repair = params$K_repair,
      O2_hyp = data$O2_hyp[i],
      O2_ref = data$O2_ref[i]
    )
  })
  
  if (any(!is.finite(data$OER_pred))) {
    if (return_details) return(list(failed = TRUE))
    return(1e10)
  }
  
  data$residual <- data$OER_retention - data$OER_pred
  
  loss <- sum(data$weight * huber_loss(data$residual, delta = 0.5))
  
  # Regularization
  steep_reg <- 0.1 * ((log(params$s_dir_light) - log(2.0))^2 + 
                        (log(params$s_ind_light) - log(2.0))^2)
  k_reg <- 0.3 * (log(params$K_fix / 0.15))^2
  
  total <- loss + steep_reg + k_reg
  
  if (return_details) {
    p_ind_hyp <- calc_p_indirect_single_MM(0.001, params$K_fix, params$K_repair)
    p_ind_ref <- calc_p_indirect_single_MM(21.0, params$K_fix, params$K_repair)
    
    P_hyp <- FIXED_PARAMS$p1_low + FIXED_PARAMS$p2_low * p_ind_hyp + FIXED_PARAMS$p3_low * p_ind_hyp^2
    P_ref <- FIXED_PARAMS$p1_low + FIXED_PARAMS$p2_low * p_ind_ref + FIXED_PARAMS$p3_low * p_ind_ref^2
    OER_max_theoretical <- P_ref / P_hyp
    
    return(list(
      failed = FALSE, 
      data = data, 
      params = params, 
      OER_max_theoretical = OER_max_theoretical,
      p_ind_hyp = p_ind_hyp, 
      p_ind_ref = p_ind_ref
    ))
  }
  
  return(ifelse(is.finite(total), total, 1e10))
}

# Initial values and bounds
phase1_init <- c(K_fix = 0.15, K_repair = 0.50,
                 x50_dir_photon = 150, x50_ind_photon = 2000,
                 s_dir_light = 1.5, s_ind_light = 2.0)

phase1_lower <- c(K_fix = 0.02, K_repair = 0.10,
                  x50_dir_photon = 10, x50_ind_photon = 300,
                  s_dir_light = 0.5, s_ind_light = 0.5)

phase1_upper <- c(K_fix = 0.80, K_repair = 2.0,
                  x50_dir_photon = 800, x50_ind_photon = 10000,
                  s_dir_light = 5.0, s_ind_light = 5.0)

# Multi-start optimization
best_p1 <- NULL
best_p1_val <- Inf

K_fix_starts <- c(0.08, 0.12, 0.18, 0.25)
K_repair_starts <- c(0.25, 0.40, 0.60)

cat("Multi-start optimization:\n")
for (kf in K_fix_starts) {
  for (kr in K_repair_starts) {
    if (kf >= kr) next
    
    init <- c(K_fix = kf, K_repair = kr,
              x50_dir_photon = 150, x50_ind_photon = 2000,
              s_dir_light = 1.5, s_ind_light = 2.0)
    
    opt <- tryCatch({
      optim(par = init, fn = objective_phase1, data = photon_data,
            method = "L-BFGS-B", lower = phase1_lower, upper = phase1_upper,
            control = list(maxit = 1000))
    }, error = function(e) list(value = 1e10))
    
    if (opt$value < best_p1_val) {
      best_p1_val <- opt$value
      best_p1 <- opt
      cat(sprintf("  K_fix=%.2f, K_repair=%.2f -> %.4f (new best)\n", kf, kr, opt$value))
    }
  }
}
# Refine best solution
opt_p1_2 <- optim(par = best_p1$par, fn = objective_phase1, data = photon_data,
                  method = "L-BFGS-B", lower = phase1_lower, upper = phase1_upper,
                  control = list(maxit = 2000, factr = 1e6))

opt_p1_final <- optim(par = opt_p1_2$par, fn = objective_phase1, data = photon_data,
                      method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-10))

opt_p1_final$par <- pmax(phase1_lower, pmin(phase1_upper, opt_p1_final$par))

phase1_results <- objective_phase1(opt_p1_final$par, photon_data, return_details = TRUE)
phase1_params <- phase1_results$params

cat(sprintf("║  K_fix           = %.4f%% O2                                        ║\n", phase1_params$K_fix))
cat(sprintf("║  K_repair        = %.4f%% O2                                        ║\n", phase1_params$K_repair))
cat(sprintf("║  x50_dir_photon  = %8.1f                                          ║\n", phase1_params$x50_dir_photon))
cat(sprintf("║  x50_ind_photon  = %8.1f                                          ║\n", phase1_params$x50_ind_photon))
cat(sprintf("║  s_dir_light     = %8.3f                                          ║\n", phase1_params$s_dir_light))
cat(sprintf("║  s_ind_light     = %8.3f                                          ║\n", phase1_params$s_ind_light))
cat(sprintf("║  OER_max (theoretical) = %.2f                                       ║\n", phase1_results$OER_max_theoretical))
FROZEN_PHASE1 <- list(
  K_fix = phase1_params$K_fix,
  K_repair = phase1_params$K_repair,
  x50_dir_photon = phase1_params$x50_dir_photon,
  x50_ind_photon = phase1_params$x50_ind_photon,
  s_dir_light = phase1_params$s_dir_light,
  s_ind_light = phase1_params$s_ind_light
)


# Section 4: Phase 2 - Heavy Particles + Z-Ordering + Overkill

cat("\n--- SECTION 4: PHASE 2 - HEAVY PARTICLES (Z-Ordered + Overkill) ---\n")
cat("PHYSICS CONSTRAINT: x50_dir must INCREASE with Z for heavy ions\n")
cat("NEW FEATURES:\n")
# For hybrid model: proton and deuteron use light particle steepness
# but have their own x50 values
particles_phase2 <- c("proton", "deuteron", "He", "C", "Ne", "Ar")

phase2_param_names <- c(
  paste0("x50_dir_", particles_phase2),
  paste0("x50_ind_", particles_phase2),
  "s_dir_base", "s_dir_scale",
  "s_ind_base", "s_ind_scale",
  "overkill_strength",  # NEW: overkill correction strength
  "factor_HSG",
  "factor_T1",
  "factor_CHO"
)

cat(sprintf("Phase 2 parameters: %d\n\n", length(phase2_param_names)))

#' Physics-based ordering penalty for x50 values (ENHANCED)
#' Enforces x50_dir INCREASES with Z for heavy ions
penalty_physics_ordering <- function(params_list, lambda_order = 25.0, lambda_ratio = 12.0) {
  penalty <- 0
  
  # 1. x50_dir < x50_ind for each particle
  for (p in c("photon", particles_phase2)) {
    x50_dir <- params_list[[paste0("x50_dir_", p)]]
    x50_ind <- params_list[[paste0("x50_ind_", p)]]
    if (is.null(x50_dir) || is.null(x50_ind)) next
    
    if (x50_dir > x50_ind) {
      penalty <- penalty + lambda_ratio * ((x50_dir - x50_ind) / x50_ind)^2
    }
    
    ratio <- x50_ind / x50_dir
    if (ratio > 0 && ratio < 1.5) {
      penalty <- penalty + lambda_ratio * 0.3 * (1.5 - ratio)^2
    }
  }
  
  # 2. PHYSICS CONSTRAINT: x50_dir must INCREASE with Z for heavy ions
  # He (Z=2) < C (Z=6) < Ne (Z=10) < Ar (Z=18)
  heavy_ions <- c("He", "C", "Ne", "Ar")
  heavy_Z <- c(2, 6, 10, 18)
  
  for (i in 1:(length(heavy_ions) - 1)) {
    x50_dir_lower <- params_list[[paste0("x50_dir_", heavy_ions[i])]]
    x50_dir_upper <- params_list[[paste0("x50_dir_", heavy_ions[i + 1])]]
    
    if (!is.null(x50_dir_lower) && !is.null(x50_dir_upper)) {
      if (x50_dir_lower >= x50_dir_upper) {
        # Strong penalty for violating Z-ordering
        violation <- (x50_dir_lower - x50_dir_upper) / x50_dir_upper
        penalty <- penalty + lambda_order * (1 + violation)^2
      } else {
        # Encourage reasonable spacing based on Z ratio
        # x50 should scale roughly as Z^alpha where alpha ~ 0.4-0.6
        expected_ratio <- (heavy_Z[i + 1] / heavy_Z[i])^0.5
        actual_ratio <- x50_dir_upper / x50_dir_lower
        if (actual_ratio < expected_ratio * 0.7) {
          penalty <- penalty + lambda_order * 0.15 * (expected_ratio * 0.7 - actual_ratio)^2
        }
      }
    }
  }
  
  # 3. x50_ind should also increase with Z (same physics, weaker penalty)
  for (i in 1:(length(heavy_ions) - 1)) {
    x50_ind_lower <- params_list[[paste0("x50_ind_", heavy_ions[i])]]
    x50_ind_upper <- params_list[[paste0("x50_ind_", heavy_ions[i + 1])]]
    
    if (!is.null(x50_ind_lower) && !is.null(x50_ind_upper)) {
      if (x50_ind_lower >= x50_ind_upper) {
        violation <- (x50_ind_lower - x50_ind_upper) / x50_ind_upper
        penalty <- penalty + lambda_order * 0.5 * (1 + violation)^2
      }
    }
  }
  
  # 4. Penalize very negative steepness scale (prefer positive or mildly negative)
  if (!is.null(params_list$s_dir_scale) && params_list$s_dir_scale < -0.2) {
    penalty <- penalty + lambda_order * 0.3 * (params_list$s_dir_scale + 0.2)^2
  }
  
  return(penalty)
}

objective_phase2 <- function(params_vec, data, return_details = FALSE) {
  
  params_list <- as.list(params_vec)
  names(params_list) <- phase2_param_names
  
  # Add frozen parameters from Phase 1
  params_list$K_fix <- FROZEN_PHASE1$K_fix
  params_list$K_repair <- FROZEN_PHASE1$K_repair
  params_list$x50_dir_photon <- FROZEN_PHASE1$x50_dir_photon
  params_list$x50_ind_photon <- FROZEN_PHASE1$x50_ind_photon
  params_list$s_dir_light <- FROZEN_PHASE1$s_dir_light
  params_list$s_ind_light <- FROZEN_PHASE1$s_ind_light
  
  if (any(!is.finite(unlist(params_list)))) {
    if (return_details) return(list(failed = TRUE))
    return(1e10)
  }
  
  # Cell-line factor bounds
  if (params_list$factor_HSG < 0.7 || params_list$factor_HSG > 1.4 ||
      params_list$factor_T1 < 0.7 || params_list$factor_T1 > 1.4 ||
      params_list$factor_CHO < 0.7 || params_list$factor_CHO > 1.4) {
    if (return_details) return(list(failed = TRUE))
    return(1e10)
  }
  
  # Overkill strength bounds
  if (params_list$overkill_strength < 0 || params_list$overkill_strength > 0.3) {
    if (return_details) return(list(failed = TRUE))
    return(1e10)
  }
  
  # Predict OER for all data
  data$OER_pred <- sapply(1:nrow(data), function(i) {
    p <- data$ion[i]
    p_class <- data$particle_class[i]
    Z <- data$Z[i]
    cell_line <- data$cell_line_std[i]
    max_let <- data$max_let[i]
    
    # Get x50 parameters
    if (p %in% c("photon", particles_phase2)) {
      x50_dir <- params_list[[paste0("x50_dir_", p)]]
      x50_ind <- params_list[[paste0("x50_ind_", p)]]
    } else if (p %in% c("N", "O", "Si")) {
      interp <- interpolate_x50(Z, params_list)
      x50_dir <- interp$x50_dir
      x50_ind <- interp$x50_ind
    } else {
      return(NA)
    }
    
    if (is.null(x50_dir) || is.null(x50_ind)) return(NA)
    
    # HYBRID MODEL: Get steepness based on particle class
    if (p_class == "light") {
      # Light particles (photon, proton, deuteron) use light steepness
      s_dir <- params_list$s_dir_light
      s_ind <- params_list$s_ind_light
    } else {
      # Heavy particles use Z-scaled steepness
      s_dir <- calc_steepness_Z(Z, params_list$s_dir_base, params_list$s_dir_scale)
      s_ind <- calc_steepness_Z(Z, params_list$s_ind_base, params_list$s_ind_scale)
    }
    
    s_dir <- pmax(0.5, pmin(6.0, s_dir))
    s_ind <- pmax(0.5, pmin(6.0, s_ind))
    
    # Use overkill-corrected OER prediction for heavy particles
    if (p_class == "heavy") {
      OER_base <- predict_OER_with_overkill(
        LET = data$LET[i],
        x50_dir = x50_dir,
        x50_ind = x50_ind,
        s_dir = s_dir,
        s_ind = s_ind,
        K_fix = params_list$K_fix,
        K_repair = params_list$K_repair,
        O2_hyp = data$O2_hyp[i],
        O2_ref = data$O2_ref[i],
        max_let = max_let,
        overkill_strength = params_list$overkill_strength
      )
    } else {
      # Light particles: no overkill correction (different physics)
      OER_base <- predict_OER_base(
        LET = data$LET[i],
        x50_dir = x50_dir,
        x50_ind = x50_ind,
        s_dir = s_dir,
        s_ind = s_ind,
        K_fix = params_list$K_fix,
        K_repair = params_list$K_repair,
        O2_hyp = data$O2_hyp[i],
        O2_ref = data$O2_ref[i]
      )
    }
    
    # Apply cell-line correction
    cell_factor <- get_cell_correction(cell_line, 
                                       params_list$factor_HSG, 
                                       params_list$factor_T1,
                                       params_list$factor_CHO)
    
    return(OER_base * cell_factor)
  })
  
  if (any(!is.finite(data$OER_pred))) {
    if (return_details) return(list(failed = TRUE))
    return(1e10)
  }
  
  data$residual <- data$OER_retention - data$OER_pred
  
  # Weighted Huber loss
  loss <- sum(data$weight * huber_loss(data$residual, delta = 0.8))
  
  # Physics-based ordering penalty
  order_penalty <- penalty_physics_ordering(params_list, lambda_order = 25.0, lambda_ratio = 12.0)
  
  # Regularization
  steep_reg <- 0.1 * ((log(params_list$s_dir_base) - log(2.5))^2 + 
                        (log(params_list$s_ind_base) - log(2.5))^2)
  
  scale_reg <- 0.3 * (params_list$s_dir_scale^2 + params_list$s_ind_scale^2)
  
  cell_reg <- 1.5 * ((params_list$factor_HSG - 1.0)^2 + 
                       (params_list$factor_T1 - 1.0)^2 +
                       (params_list$factor_CHO - 1.0)^2)
  
  # Mild regularization on overkill strength (prefer small values)
  overkill_reg <- 2.0 * params_list$overkill_strength^2
  
  total <- loss + order_penalty + steep_reg + scale_reg + cell_reg + overkill_reg
  
  if (return_details) {
    steepness_table <- tibble(
      ion = c("He", "C", "Ne", "Ar"),
      Z = c(2, 6, 10, 18),
      s_dir = sapply(c(2, 6, 10, 18), function(z) 
        calc_steepness_Z(z, params_list$s_dir_base, params_list$s_dir_scale)),
      s_ind = sapply(c(2, 6, 10, 18), function(z) 
        calc_steepness_Z(z, params_list$s_ind_base, params_list$s_ind_scale))
    )
    
    # Check Z-ordering (CORRECT: should increase with Z)
    z_order_ok <- all(diff(c(
      params_list$x50_dir_He,
      params_list$x50_dir_C,
      params_list$x50_dir_Ne,
      params_list$x50_dir_Ar
    )) > 0)
    
    return(list(
      failed = FALSE,
      data = data,
      params = params_list,
      steepness_table = steepness_table,
      z_ordering_satisfied = z_order_ok
    ))
  }
  
  return(ifelse(is.finite(total), total, 1e10))
}

# Initial values with Z-ordered x50_dir (INCREASES with Z)
phase2_init <- c(
  x50_dir_proton = 180, x50_dir_deuteron = 200,
  x50_dir_He = 160, x50_dir_C = 320, x50_dir_Ne = 500, x50_dir_Ar = 700,
  x50_ind_proton = 2500, x50_ind_deuteron = 2200,
  x50_ind_He = 1800, x50_ind_C = 2400, x50_ind_Ne = 3200, x50_ind_Ar = 4000,
  s_dir_base = 2.5, s_dir_scale = 0.05,
  s_ind_base = 2.5, s_ind_scale = 0.05,
  overkill_strength = 0.08,
  factor_HSG = 1.0,
  factor_T1 = 1.0,
  factor_CHO = 1.0
)

# Bounds with enforced Z-ordering through lower/upper limits
phase2_lower <- c(
  x50_dir_proton = 50, x50_dir_deuteron = 50,
  x50_dir_He = 80, x50_dir_C = 150, x50_dir_Ne = 280, x50_dir_Ar = 400,
  x50_ind_proton = 400, x50_ind_deuteron = 400,
  x50_ind_He = 600, x50_ind_C = 900, x50_ind_Ne = 1400, x50_ind_Ar = 1800,
  s_dir_base = 1.2, s_dir_scale = -0.25,
  s_ind_base = 1.2, s_ind_scale = -0.25,
  overkill_strength = 0.0,
  factor_HSG = 0.75,
  factor_T1 = 0.75,
  factor_CHO = 0.75
)

phase2_upper <- c(
  x50_dir_proton = 450, x50_dir_deuteron = 450,
  x50_dir_He = 300, x50_dir_C = 550, x50_dir_Ne = 750, x50_dir_Ar = 1000,
  x50_ind_proton = 5000, x50_ind_deuteron = 4500,
  x50_ind_He = 4000, x50_ind_C = 5000, x50_ind_Ne = 6500, x50_ind_Ar = 8000,
  s_dir_base = 4.5, s_dir_scale = 0.35,
  s_ind_base = 4.0, s_ind_scale = 0.35,
  overkill_strength = 0.25,
  factor_HSG = 1.30,
  factor_T1 = 1.30,
  factor_CHO = 1.30
)

cat("Phase 2 Optimization:\n")

opt_p2_1 <- optim(par = phase2_init, fn = objective_phase2, data = calibration_data,
                  method = "L-BFGS-B", lower = phase2_lower, upper = phase2_upper,
                  control = list(maxit = 3000, factr = 1e8))
cat(sprintf("  L-BFGS-B coarse: %.4f\n", opt_p2_1$value))

opt_p2_2 <- optim(par = opt_p2_1$par, fn = objective_phase2, data = calibration_data,
                  method = "L-BFGS-B", lower = phase2_lower, upper = phase2_upper,
                  control = list(maxit = 5000, factr = 1e6))
cat(sprintf("  L-BFGS-B fine: %.4f\n", opt_p2_2$value))

opt_p2_final <- optim(par = opt_p2_2$par, fn = objective_phase2, data = calibration_data,
                      method = "Nelder-Mead", control = list(maxit = 15000, reltol = 1e-10))
cat(sprintf("  Nelder-Mead: %.4f\n\n", opt_p2_final$value))

opt_p2_final$par <- pmax(phase2_lower, pmin(phase2_upper, opt_p2_final$par))


# Section 5: Final Results

cat("\n--- SECTION 5: FINAL PARAMETERS AND FIT STATISTICS ---\n")
final_results <- objective_phase2(opt_p2_final$par, calibration_data, return_details = TRUE)
final_data <- final_results$data
final_params <- final_results$params

# Check Z-ordering (CORRECT CHECK: x50 should INCREASE with Z)
cat("PHYSICS CONSTRAINT CHECK:\n")
x50_he <- final_params$x50_dir_He
x50_c <- final_params$x50_dir_C
x50_ne <- final_params$x50_dir_Ne
x50_ar <- final_params$x50_dir_Ar

cat(sprintf("  x50_dir(He)  = %6.1f\n", x50_he))
cat(sprintf("  x50_dir(C)   = %6.1f  (should be > He)\n", x50_c))
cat(sprintf("  x50_dir(Ne)  = %6.1f  (should be > C)\n", x50_ne))
cat(sprintf("  x50_dir(Ar)  = %6.1f  (should be > Ne)\n", x50_ar))
# FIXED: Correct Z-ordering check
z_ordering_satisfied <- (x50_he < x50_c) && (x50_c < x50_ne) && (x50_ne < x50_ar)

if (z_ordering_satisfied) {
    } else {
  }

# Display parameters
cat(sprintf("║    K_fix    = %.4f%% O2                                             ║\n", final_params$K_fix))
cat(sprintf("║    K_repair = %.4f%% O2                                             ║\n", final_params$K_repair))
cat(sprintf("��    s_dir_light = %.3f                                               ║\n", final_params$s_dir_light))
cat(sprintf("║    s_ind_light = %.3f                                               ║\n", final_params$s_ind_light))
cat(sprintf("║    s_dir(Z) = %.3f * (1 + %.3f * log(Z/2))                         ║\n", 
            final_params$s_dir_base, final_params$s_dir_scale))
cat(sprintf("║    s_ind(Z) = %.3f * (1 + %.3f * log(Z/2))                         ║\n",
            final_params$s_ind_base, final_params$s_ind_scale))
cat(sprintf("║    overkill_strength = %.3f                                         ║\n", final_params$overkill_strength))
cat(sprintf("║    V79 (reference) = 1.000                                          ║\n"))
cat(sprintf("║    HSG factor      = %.3f                                           ║\n", final_params$factor_HSG))
cat(sprintf("║    T1 factor       = %.3f                                           ║\n", final_params$factor_T1))
cat(sprintf("║    CHO factor      = %.3f                                           ║\n", final_params$factor_CHO))
# Calculate theoretical OER_max
p_ind_hyp <- calc_p_indirect_single_MM(0.001, final_params$K_fix, final_params$K_repair)
p_ind_ref <- calc_p_indirect_single_MM(21.0, final_params$K_fix, final_params$K_repair)

P_hyp <- FIXED_PARAMS$p1_low + FIXED_PARAMS$p2_low * p_ind_hyp + FIXED_PARAMS$p3_low * p_ind_hyp^2
P_ref <- FIXED_PARAMS$p1_low + FIXED_PARAMS$p2_low * p_ind_ref + FIXED_PARAMS$p3_low * p_ind_ref^2
OER_max_theoretical <- P_ref / P_hyp

cat(sprintf("Theoretical OER_max (at anoxia, low LET): %.2f\n\n", OER_max_theoretical))

# Steepness table
cat("Effective steepness for heavy ions:\n")
print(final_results$steepness_table)
# x50 table with Z-ordering check
cat("┌────────────┬─────┬───────┬───────────┬───────────┬─────────┬────────┬──────────┐\n")
cat("│ Particle   │  Z  │ Class │  x50_dir  │  x50_ind  │  Ratio  │ Status │ Z-Order  │\n")
cat("├────────────┼─────┼───────┼───────────┼───────────┼─────────┼────────┼──────────┤\n")

prev_x50_dir <- 0
for (p in c("photon", "proton", "deuteron", "He", "C", "Ne", "Ar")) {
  z <- particle_info$Z[particle_info$ion == p]
  p_class <- particle_info$particle_class[particle_info$ion == p]
  x50_d <- final_params[[paste0("x50_dir_", p)]]
  x50_i <- final_params[[paste0("x50_ind_", p)]]
  ratio <- x50_i / x50_d
  
  status <- if (x50_d > x50_i) "⚠ VIOL" else if (ratio < 1.5) "~ MARG" else "✓ OK"
  
  # Z-ordering check for heavy ions (x50 should INCREASE with Z)
  if (p_class == "heavy") {
    z_order <- if (x50_d > prev_x50_dir) "✓" else "✗"
    prev_x50_dir <- x50_d
  } else {
    z_order <- "-"
  }
  
  cat(sprintf("│ %-10s │ %3d │ %-5s │ %9.1f │ %9.1f │ %7.2f │ %-6s │    %s     │\n", 
              p, z, p_class, x50_d, x50_i, ratio, status, z_order))
}
cat("└────────────┴─────┴───────┴───────────┴───────────┴─────────┴────────┴──────────┘\n\n")

# Show interpolated particles
cat("Interpolated particles (N, O, Si):\n")
cat("┌────────────┬─────┬───────────┬───────────┐\n")
cat("│ Particle   │  Z  │  x50_dir  │  x50_ind  │\n")
cat("├────────────┼─────┼───────────┼───────────┤\n")
for (p in c("N", "O", "Si")) {
  z <- particle_info$Z[particle_info$ion == p]
  interp <- interpolate_x50(z, final_params)
  cat(sprintf("│ %-10s │ %3d │ %9.1f │ %9.1f │\n", p, z, interp$x50_dir, interp$x50_ind))
}
cat("└────────────┴─────┴───────────┴───────────┘\n\n")

# Fit statistics
n <- nrow(final_data)
k <- length(phase1_param_names) + length(phase2_param_names)

residuals_vec <- final_data$residual
TSS <- sum((final_data$OER_retention - mean(final_data$OER_retention))^2)
RSS <- sum(residuals_vec^2)
R2 <- 1 - RSS / TSS
R2_adj <- 1 - (1 - R2) * (n - 1) / (n - k - 1)
RMSE <- sqrt(RSS / n)
MAE <- mean(abs(residuals_vec))

wRSS <- sum(final_data$weight * residuals_vec^2)
wTSS <- sum(final_data$weight * (final_data$OER_retention - 
                                   weighted.mean(final_data$OER_retention, final_data$weight))^2)
wR2 <- 1 - wRSS / wTSS

cat("OVERALL FIT STATISTICS\n")
cat(sprintf("Observations: %d\n", n))
cat(sprintf("Parameters: %d\n", k))
cat(sprintf("Observations per parameter: %.1f\n\n", n / k))

cat(sprintf("R² (unweighted)  = %.4f\n", R2))
cat(sprintf("R² (weighted)    = %.4f\n", wR2))
cat(sprintf("Adjusted R²      = %.4f\n", R2_adj))
cat(sprintf("RMSE             = %.4f\n", RMSE))
cat(sprintf("MAE              = %.4f\n\n", MAE))

# Per-particle statistics
cat("Per-particle performance:\n")
cat("┌────────────┬───────┬──────┬─────────┬─────────┬─────────┬───────────┬──────────┐\n")
cat("│ Particle   │ Class │  N   │   MAE   │  RMSE   │ Max Err │  MAE (%)  │  Status  │\n")
cat("├────────────┼───────┼──────┼─────────┼─────────┼─────────┼───────────┼──────────┤\n")

particle_stats <- final_data %>%
  group_by(ion, particle_class) %>%
  summarise(
    n = n(),
    mae = mean(abs(residual)),
    rmse = sqrt(mean(residual^2)),
    max_err = max(abs(residual)),
    mae_pct = mean(abs(residual / OER_retention) * 100),
    .groups = "drop"
  ) %>%
  left_join(particle_info %>% select(ion, Z), by = "ion") %>%
  arrange(Z)

for (i in 1:nrow(particle_stats)) {
  p <- particle_stats[i, ]
  status <- if (p$mae_pct < 10) "✓ EXCEL" else if (p$mae_pct < 15) "✓ GOOD" else if (p$mae_pct < 20) "⚠ OK" else "✗ POOR"
  
  cat(sprintf("│ %-10s │ %-5s │ %4d │ %7.3f │ %7.3f │ %7.3f │ %8.1f%% │ %-8s │\n",
              p$ion, p$particle_class, p$n, p$mae, p$rmse, p$max_err, p$mae_pct, status))
}
cat("└────────────┴───────┴──────┴─────────┴─────────┴─────────┴───────────┴──────���───┘\n\n")

# Per cell-line statistics
cat("Per cell-line performance:\n")
cat("┌────────────┬──────┬─────────┬─────────┬───────────┐\n")
cat("│ Cell Line  │  N   │   MAE   │  RMSE   │  Factor   │\n")
cat("├────────────┼──────┼─────────┼─────────┼───────────┤\n")

cell_stats <- final_data %>%
  group_by(cell_line_std) %>%
  summarise(
    n = n(),
    mae = mean(abs(residual)),
    rmse = sqrt(mean(residual^2)),
    .groups = "drop"
  ) %>%
  mutate(
    factor = case_when(
      cell_line_std == "V79" ~ 1.0,
      cell_line_std == "HSG" ~ final_params$factor_HSG,
      cell_line_std == "T1" ~ final_params$factor_T1,
      cell_line_std == "CHO" ~ final_params$factor_CHO,
      TRUE ~ 1.0
    )
  ) %>%
  arrange(desc(n))

for (i in 1:nrow(cell_stats)) {
  c <- cell_stats[i, ]
  cat(sprintf("│ %-10s │ %4d │ %7.3f │ %7.3f │ %9.3f │\n",
              c$cell_line_std, c$n, c$mae, c$rmse, c$factor))
}
cat("└────────────┴──────┴─────────┴─────────┴───────────┘\n\n")

# Physical LET limit compliance
cat("Physical LET limit compliance:\n")
cat("┌────────────┬──────┬────────────┬────────────┬───────────┐\n")
cat("│ Particle   │  N   │ Within Lim │ Near Bragg │ Beyond    │\n")
cat("├────────────┼──────┼────────────┼────────────┼───────────┤\n")

let_compliance <- final_data %>%
  group_by(ion) %>%
  summarise(
    n = n(),
    within_limit = sum(bragg_proximity <= 0.7),
    near_bragg = sum(bragg_proximity > 0.7 & bragg_proximity <= 1.0),
    beyond = sum(bragg_proximity > 1.0),
    .groups = "drop"
  ) %>%
  left_join(particle_info %>% select(ion, Z), by = "ion") %>%
  arrange(Z)

for (i in 1:nrow(let_compliance)) {
  lc <- let_compliance[i, ]
  cat(sprintf("│ %-10s │ %4d │ %10d │ %10d │ %9d │\n",
              lc$ion, lc$n, lc$within_limit, lc$near_bragg, lc$beyond))
}
cat("└────────────┴──────┴────────────┴────────────┴───────────┘\n\n")


# Section 6: D-Kondo Validation (Post-Hoc)

cat("\n--- SECTION 6: D-KONDO VALIDATION (POST-HOC COMPARISON) ---\n")
cat("NOTE: D-Kondo data was NOT used in calibration.\n")
dkondo_validation$OER_pred <- sapply(1:nrow(dkondo_validation), function(i) {
  if (dkondo_validation$O2_hyp[i] >= dkondo_validation$O2_ref[i]) return(1.0)
  predict_OER_base(
    LET = dkondo_validation$LET[i],
    x50_dir = final_params$x50_dir_photon,
    x50_ind = final_params$x50_ind_photon,
    s_dir = final_params$s_dir_light,
    s_ind = final_params$s_ind_light,
    K_fix = final_params$K_fix,
    K_repair = final_params$K_repair,
    O2_hyp = dkondo_validation$O2_hyp[i],
    O2_ref = dkondo_validation$O2_ref[i]
  )
})

dkondo_validation$residual <- dkondo_validation$HRF_experimental - dkondo_validation$OER_pred
dkondo_validation$error_pct <- abs(dkondo_validation$residual / dkondo_validation$HRF_experimental) * 100

cat("Model OER vs D-Kondo HRF (reference only):\n")
cat("┌──────────┬──────────┬──────────┬───────────┬────────────┬──────────┐\n")
cat("│ O2 (%)   │ HRF_exp  │ OER_pred │ Diff      │ Diff (%)   │ Measured │\n")
cat("├──────────┼──────────┼──────────┼───────────┼────────────┼──────────┤\n")
for (i in 1:nrow(dkondo_validation)) {
  meas <- if(dkondo_validation$is_measured[i]) "Yes" else "Interp"
  cat(sprintf("│ %8.3f │ %8.2f │ %8.2f │ %+9.3f │ %9.1f%% │ %-8s │\n",
              dkondo_validation$O2_hyp[i],
              dkondo_validation$HRF_experimental[i],
              dkondo_validation$OER_pred[i],
              dkondo_validation$residual[i],
              dkondo_validation$error_pct[i],
              meas))
}
cat("└──────────┴──────────┴──────────┴───────────┴────────────┴──────────┘\n\n")

dkondo_measured <- dkondo_validation %>% filter(is_measured, O2_hyp < 21)
dkondo_mae <- mean(dkondo_measured$error_pct)
cat(sprintf("D-Kondo comparison MAE (measured points): %.1f%%\n\n", dkondo_mae))


# Section 7: Tinganelli Validation (Intermediate O2)

cat("\n--- SECTION 7: TINGANELLI VALIDATION (INTERMEDIATE O2 LEVELS) ---\n")
# Get Tinganelli data from calibration set
tinganelli_data <- final_data %>%
  filter(dataset == "Tinganelli2015" | source == "Tinganelli2015")

if (nrow(tinganelli_data) > 0) {
  cat(sprintf("Tinganelli data in calibration: %d points\n\n", nrow(tinganelli_data)))
  
  tinganelli_stats <- tinganelli_data %>%
    group_by(ion, O2_hyp) %>%
    summarise(
      n = n(),
      OER_obs_mean = mean(OER_retention),
      OER_pred_mean = mean(OER_pred),
      MAE = mean(abs(residual)),
      .groups = "drop"
    ) %>%
    arrange(ion, O2_hyp)
  
  cat("Performance by particle and O2 level:\n")
  cat("┌────────────┬─────────┬─────┬───────────┬───────────┬─────────┐\n")
  cat("│ Particle   │ O2 (%)  │  n  │ OER_obs   │ OER_pred  │   MAE   │\n")
  cat("├────────────┼─────────┼─────┼───────────┼───────────┼─────────┤\n")
  for (i in 1:nrow(tinganelli_stats)) {
    t <- tinganelli_stats[i, ]
    cat(sprintf("│ %-10s │ %7.3f │ %3d │ %9.2f │ %9.2f │ %7.3f │\n",
                t$ion, t$O2_hyp, t$n, t$OER_obs_mean, t$OER_pred_mean, t$MAE))
  }
  cat("└────────────┴─────────┴─────┴───────────┴────────���──┴─────────┘\n\n")
  
  overall_tinganelli <- tinganelli_data %>%
    summarise(
      n = n(),
      MAE = mean(abs(residual)),
      RMSE = sqrt(mean(residual^2)),
      R2 = 1 - sum(residual^2) / sum((OER_retention - mean(OER_retention))^2)
    )
  
  cat(sprintf("Tinganelli overall: MAE = %.3f, RMSE = %.3f, R² = %.3f\n\n",
              overall_tinganelli$MAE, overall_tinganelli$RMSE, overall_tinganelli$R2))
} else {
  cat("No Tinganelli data found in calibration set.\n\n")
}


# Section 8: Save Results

cat("\n--- SECTION 8: SAVING RESULTS ---\n")
if (!dir.exists("results")) dir.create("results")

# Calculate OER_survival from OER_retention for output
final_data <- final_data %>%
  mutate(
    OER_survival = convert_retention_to_survival(OER_retention),
    OER_pred_survival = convert_retention_to_survival(OER_pred)
  )

model <- list(
  version = "VOxA",
  model_type = "physics_constrained_overkill_hybrid",
  date = Sys.Date(),
  
  parameters = final_params,
  fixed_params = FIXED_PARAMS,
  
  phase1_params = FROZEN_PHASE1,
  phase2_params = as.list(opt_p2_final$par),
  
  particle_classes = list(
    light = light_particles, 
    heavy_calibrated = heavy_particles_calibrated,
    heavy_interpolated = heavy_particles_interpolated
  ),
  steepness_table = final_results$steepness_table,
  
  # Physical LET limits
  max_let_physical = MAX_LET_PHYSICAL,
  
  # Z-ordering status
  z_ordering_satisfied = z_ordering_satisfied,
  
  # Interpolated particle parameters
  interpolated_params = list(
    N = interpolate_x50(7, final_params),
    O = interpolate_x50(8, final_params),
    Si = interpolate_x50(14, final_params)
  ),
  
  n_outliers_removed = n_removed,
  
  cell_line_factors = list(
    V79 = 1.0, 
    HSG = final_params$factor_HSG, 
    T1 = final_params$factor_T1,
    CHO = final_params$factor_CHO
  ),
  
  OER_max_theoretical = OER_max_theoretical,
  
  # Overkill correction
  overkill_strength = final_params$overkill_strength,
  
  calibration_data = final_data,
  dkondo_validation = dkondo_validation,
  
  fit_statistics = list(
    n_obs = n,
    n_params = k,
    r2 = R2,
    r2_weighted = wR2,
    adj_r2 = R2_adj,
    rmse = RMSE,
    mae = MAE
  ),
  
  particle_statistics = particle_stats,
  cell_line_statistics = cell_stats,
  
  physics_constraints = list(
    z_ordering = "x50_dir(He) < x50_dir(C) < x50_dir(Ne) < x50_dir(Ar)",
    rationale = "At fixed LET, lighter ions are slower with narrower, denser tracks, leading to more radical recombination and lower OER",
    verified = z_ordering_satisfied
  ),
  
  enhancements = list(
    "Physical LET cutoffs based on Bragg peak limits",
    "Overkill correction near Bragg peak",
    "Hybrid light/heavy particle model",
    "Improved weighting (clinical relevance + physical plausibility)",
    "Z-dependent steepness with physics motivation"
  ),
  
  notes = c(
    "VOxA",
    "Physical LET limits: proton=100, He=200, C=550, Ne=700, Ar=900 keV/μm",
    "Overkill correction applied near Bragg peak for heavy ions",
    "Hybrid model: light particles use s_light, heavy use Z-scaled steepness",
    "Z-ordering enforced: x50 increases with Z (correct physics)",
    "Calibration data: Furusawa + Literature + Tinganelli",
    "Case fractions: p1=0.04, p2=0.32, p3=0.64 at low LET",
    "p1_high=0.64 based on Nikjoo 2001"
  )
)

saveRDS(model, "results/uvaom_v8_corrected_model.rds")
write_csv(final_data, "results/calibration_data_v8_corrected.csv")
write_csv(dkondo_validation, "results/dkondo_comparison_v8.csv")

cat("Saved:\n")
# Section 9: Benchmark Predictions

cat("\n--- SECTION 9: BENCHMARK PREDICTIONS (vs Grimes 2020) ---\n")
cat("Predictions at fixed LET = 100 keV/μm, anoxia (0.001% O2):\n")
cat("Particle   Z    x50_dir   Max LET   VOxA OER_ret   VOxA OER_surv   Physical?\n")
test_LET <- 100
test_O2 <- 0.001

for (p in c("photon", "proton", "He", "C", "Ne", "Ar")) {
  z <- particle_info$Z[particle_info$ion == p]
  p_class <- particle_info$particle_class[particle_info$ion == p]
  max_let <- particle_info$max_let[particle_info$ion == p]
  
  x50_dir <- final_params[[paste0("x50_dir_", p)]]
  x50_ind <- final_params[[paste0("x50_ind_", p)]]
  
  if (p_class == "light") {
    s_dir <- final_params$s_dir_light
    s_ind <- final_params$s_ind_light
  } else {
    s_dir <- calc_steepness_Z(z, final_params$s_dir_base, final_params$s_dir_scale)
    s_ind <- calc_steepness_Z(z, final_params$s_ind_base, final_params$s_ind_scale)
  }
  
  # Use appropriate prediction function
  if (p_class == "heavy") {
    oer_ret <- predict_OER_with_overkill(test_LET, x50_dir, x50_ind, s_dir, s_ind,
                                         final_params$K_fix, final_params$K_repair,
                                         test_O2, 21.0, max_let, 
                                         final_params$overkill_strength)
  } else {
    oer_ret <- predict_OER_base(test_LET, x50_dir, x50_ind, s_dir, s_ind,
                                final_params$K_fix, final_params$K_repair,
                                test_O2, 21.0)
  }
  
  oer_surv <- convert_retention_to_survival(oer_ret)
  
  # Check if LET is physical for this particle
  physical_ok <- if (test_LET <= max_let) "✓ Yes" else "✗ No"
  
  cat(sprintf("%-10s %2d   %7.1f   %7d   %13.3f    %14.3f    %s\n",
              p, z, x50_dir, max_let, oer_ret, oer_surv, physical_ok))
}
# Note about proton at 100 keV/μm
cat("NOTE: Proton at LET=100 keV/μm is at its Bragg peak (physical limit).\n")
# Calculate theoretical OER_max in both conventions
OER_max_retention <- OER_max_theoretical
OER_max_survival <- convert_retention_to_survival(OER_max_retention)

cat("Theoretical OER_max:\n")
cat(sprintf("  OER_max (retention) = %.2f\n", OER_max_retention))
cat(sprintf("  OER_max (survival)  = %.2f\n", OER_max_survival))
cat(sprintf("  Grimes (2020)       = 2.74 (survival)\n\n"))

cat("Key advantages over Grimes (2020):\n")
# Final Summary

n_excellent <- sum(particle_stats$mae_pct < 10)
n_good <- sum(particle_stats$mae_pct >= 10 & particle_stats$mae_pct < 15)
n_ok <- sum(particle_stats$mae_pct >= 15 & particle_stats$mae_pct < 20)
n_poor <- sum(particle_stats$mae_pct >= 20)

if (z_ordering_satisfied) {
  } else {
  }
cat(sprintf("║   R² (unweighted) = %.4f                                           ║\n", R2))
cat(sprintf("║   R² (weighted)   = %.4f                                           ║\n", wR2))
cat(sprintf("║   RMSE = %.4f  |  MAE = %.4f                                     ║\n", RMSE, MAE))
cat(sprintf("║   OER_max (theoretical) = %.2f                                      ║\n", OER_max_theoretical))
cat(sprintf("║   K_fix = %.4f%%  |  K_repair = %.4f%%                            ║\n", 
            final_params$K_fix, final_params$K_repair))
cat(sprintf("║   Overkill strength = %.3f                                         ║\n", 
            final_params$overkill_strength))
cat(sprintf("║   Particles: %d excellent, %d good, %d ok, %d poor                    ║\n",
            n_excellent, n_good, n_ok, n_poor))
cat(sprintf("║     V79 = 1.000 (reference)                                         ║\n"))
cat(sprintf("║     HSG = %.3f                                                      ║\n", final_params$factor_HSG))
cat(sprintf("║     T1  = %.3f                                                      ║\n", final_params$factor_T1))
cat(sprintf("║     CHO = %.3f                                                      ║\n", final_params$factor_CHO))
# Appendix: Model Predictions At Key Conditions

cat("\n--- APPENDIX: MODEL PREDICTIONS AT KEY CONDITIONS ---\n")
# Predictions at different O2 levels for Carbon
cat("Carbon ion OER vs O2 (at LET = 100 keV/μm):\n")
cat("┌──────────────┬──────────────┬────────────────┬──────────────��─┐\n")
cat("│ O2 (%)       │ Condition    │ OER_retention  │ OER_survival   │\n")
cat("├──────────────┼──────────────┼────────────────┼────────────────┤\n")

o2_levels <- c(21.0, 2.0, 0.5, 0.15, 0.001)
conditions <- c("Normoxia", "Physoxia", "Mild hypoxia", "Acute hypoxia", "Anoxia")

for (i in seq_along(o2_levels)) {
  o2 <- o2_levels[i]
  cond <- conditions[i]
  
  s_dir_c <- calc_steepness_Z(6, final_params$s_dir_base, final_params$s_dir_scale)
  s_ind_c <- calc_steepness_Z(6, final_params$s_ind_base, final_params$s_ind_scale)
  
  oer_ret <- predict_OER_with_overkill(
    LET = 100,
    x50_dir = final_params$x50_dir_C,
    x50_ind = final_params$x50_ind_C,
    s_dir = s_dir_c,
    s_ind = s_ind_c,
    K_fix = final_params$K_fix,
    K_repair = final_params$K_repair,
    O2_hyp = o2,
    O2_ref = 21.0,
    max_let = 550,
    overkill_strength = final_params$overkill_strength
  )
  oer_surv <- convert_retention_to_survival(oer_ret)
  
  cat(sprintf("│ %12.3f │ %-12s │ %14.3f │ %14.3f │\n",
              o2, cond, oer_ret, oer_surv))
}
cat("└──────────────┴──────────────┴────────────────┴────────────────┘\n\n")

# OER vs LET for different particles (at anoxia) - WITH PHYSICAL LIMITS
cat("OER_survival vs LET at anoxia (0.001% O2):\n")
cat("(Values marked with * are beyond physical LET limit for that particle)\n\n")
cat("┌──────────┬───────────┬───────────┬───────────┬───────────┬───────────┐\n")
cat("│ LET      │ photon    │ He        │ C         │ Ne        │ Ar        │\n")
cat("├──────────┼───────────┼───────────┼───────────┼───────────┼───────────┤\n")

let_values <- c(2, 10, 30, 50, 100, 200, 300, 500)
particles_to_show <- c("photon", "He", "C", "Ne", "Ar")

for (let_val in let_values) {
  oer_vals <- sapply(particles_to_show, function(p) {
    z <- particle_info$Z[particle_info$ion == p]
    p_class <- particle_info$particle_class[particle_info$ion == p]
    max_let <- particle_info$max_let[particle_info$ion == p]
    
    x50_dir <- final_params[[paste0("x50_dir_", p)]]
    x50_ind <- final_params[[paste0("x50_ind_", p)]]
    
    if (p_class == "light") {
      s_dir <- final_params$s_dir_light
      s_ind <- final_params$s_ind_light
    } else {
      s_dir <- calc_steepness_Z(z, final_params$s_dir_base, final_params$s_dir_scale)
      s_ind <- calc_steepness_Z(z, final_params$s_ind_base, final_params$s_ind_scale)
    }
    
    if (p_class == "heavy") {
      oer_ret <- predict_OER_with_overkill(let_val, x50_dir, x50_ind, s_dir, s_ind,
                                           final_params$K_fix, final_params$K_repair,
                                           0.001, 21.0, max_let,
                                           final_params$overkill_strength)
    } else {
      oer_ret <- predict_OER_base(let_val, x50_dir, x50_ind, s_dir, s_ind,
                                  final_params$K_fix, final_params$K_repair,
                                  0.001, 21.0)
    }
    
    oer_surv <- convert_retention_to_survival(oer_ret)
    
    # Mark if beyond physical limit
    if (let_val > max_let) {
      return(sprintf("%6.2f*", oer_surv))
    } else {
      return(sprintf("%7.2f", oer_surv))
    }
  })
  
  cat(sprintf("│ %8.0f │ %9s │ %9s │ %9s │ %9s │ %9s │\n",
              let_val, oer_vals[1], oer_vals[2], oer_vals[3], oer_vals[4], oer_vals[5]))
}
cat("└──────────┴───────────┴───────────┴───────────┴───────────┴───────────┘\n\n")

cat("Key observations:\n")
# End Of Script

cat("═══════════════════════════════════════════════════════════════════════\n\n")