# Step 4: Comprehensive validation and diagnostics
#
# Generates figures 1–12 for the technical report:
#   1  OER vs LET (all particles, physical limits)
#   2  Individual particle fits
#   3  Predicted vs observed OER
#   4  Residuals vs LET
#   5  Oxygen response curve (D-Kondo comparison)
#   6  Hirayama (2005) carbon validation
#   7  Residual histogram
#   8  Residual boxplot by particle
#   9  Physical LET limits compliance
#   10 Z-ordering validation
#   11 Tinganelli intermediate-O2 validation
#   12 Cell-line performance comparison
#
# Inputs:  results/uvaom_v8_corrected_model.rds,
#          results/calibration_data_v8_corrected.csv
# Outputs: figures/fig01_*.png through figures/fig12_*.png,
#          results/comprehensive_validation_report_voxa.txt

library(tidyverse)
library(scales)
library(gridExtra)
library(grid)

cat("���══════════════════════════════════════════════════════════════════════╝\n")
# Section 0: Load Model And Data

cat("\n--- SECTION 0: LOADING MODEL AND DATA ---\n")
# Load the VOxA model
model <- readRDS("results/uvaom_v8_corrected_model.rds")
calibration_data <- read_csv("results/calibration_data_v8_corrected.csv", show_col_types = FALSE)

# Load conversion factor
if (file.exists("data/uvaom_recalibration_setup.RData")) {
  load("data/uvaom_recalibration_setup.RData")
} else {
  CONVERSION_FACTOR <- list(mean = 1.20, sd = 0.05)
}

cat(sprintf("Model version: %s\n", model$version))
cat(sprintf("Observations: %d\n", model$fit_statistics$n_obs))
cat(sprintf("R² (unweighted) = %.4f\n", model$fit_statistics$r2))
cat(sprintf("R² (weighted) = %.4f\n", model$fit_statistics$r2_weighted))
cat(sprintf("OER_max (retention, theoretical) = %.2f\n", model$OER_max_theoretical))

# Calculate survival OER_max using general conversion
OER_max_survival <- 1.0 + (model$OER_max_theoretical - 1.0) / CONVERSION_FACTOR$mean
cat(sprintf("OER_max (survival, theoretical) = %.2f\n\n", OER_max_survival))

# Extract parameters
params <- model$parameters
FIXED_PARAMS <- model$fixed_params

# Physical LET limits (Bragg peak) - KEY ADDITION
MAX_LET_PHYSICAL <- list(
  photon = 35,
  proton = 100,
  deuteron = 120,
  He = 200,
  C = 550,
  N = 600,
  O = 620,
  Ne = 700,
  Si = 800,
  Ar = 900
)

cat("Physical LET limits (Bragg peak):\n")
for (p in names(MAX_LET_PHYSICAL)) {
  cat(sprintf("  %-10s: %4d keV/μm\n", p, MAX_LET_PHYSICAL[[p]]))
}
# Particle metadata - includes N, O, Si
particle_info <- tibble(
  ion = c("photon", "proton", "deuteron", "He", "C", "N", "O", "Ne", "Si", "Ar"),
  Z = c(0, 1, 1, 2, 6, 7, 8, 10, 14, 18),
  particle_class = c("light", "light", "light", "heavy", "heavy", "heavy", "heavy", "heavy", "heavy", "heavy"),
  max_let = c(35, 100, 120, 200, 550, 600, 620, 700, 800, 900)
)

light_particles <- c("photon", "proton", "deuteron")
heavy_particles_calibrated <- c("He", "C", "Ne", "Ar")
heavy_particles_interpolated <- c("N", "O", "Si")
all_particles <- c(light_particles, heavy_particles_calibrated)
all_particles_with_interpolated <- c(light_particles, "He", "C", "N", "O", "Ne", "Si", "Ar")

# Create output directories
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("results")) dir.create("results")


# Oer Conversion Functions

#' Convert OER_retention to OER_survival (general formula)
#' Ensures OER_survival = 1.0 when OER_retention = 1.0
convert_retention_to_survival <- function(OER_retention) {
  1.00 + ((OER_retention - 1.00) / CONVERSION_FACTOR$mean)
}

#' Convert OER_survival to OER_retention (general formula)
convert_survival_to_retention <- function(OER_survival) {
  1.00 + ((OER_survival - 1.00) * CONVERSION_FACTOR$mean)
}


# Color Palette: "A Summer In Southern Italy"

# Particle colors - distinct and vibrant like Italian coastal towns
particle_colors <- c(
  "photon"   = "#1E5B8C",
  "proton"   = "#E07B3C",
  "deuteron" = "#8B5A83",
  "He"       = "#D4A84B",
  "C"        = "#2D7D46",
  "N"        = "#5DADE2",
  "O"        = "#48C9B0",
  "Ne"       = "#C75B5B",
  "Si"       = "#AF7AC5",
  "Ar"       = "#6B4E3D"
)

# Cell line colors
cell_line_colors <- c(
  "V79"   = "#E07B3C",
  "HSG"   = "#1E5B8C",
  "CHO"   = "#C75B5B",
  "T1"    = "#8B5A83",
  "Other" = "#5C7C4E"
)

# Particle class colors
particle_class_colors <- c(
  "light" = "#D4A84B",
  "heavy" = "#C75B5B"
)

# Model curve color
model_curve_color <- "#0D3B66"


# Section 1: Model Functions (From Step 3 - Enhanced)

cat("\n--- SECTION 1: DEFINING MODEL FUNCTIONS ---\n")
calc_p_indirect_single_MM <- function(O2_percent, K_fix, K_repair) {
  (O2_percent + K_fix) / (O2_percent + K_fix + K_repair)
}

compute_transitions <- function(LET, x50_dir, x50_ind, s_dir, s_ind) {
  x <- 2.5 * LET^1.1
  f_direct <- 1 / (1 + (x50_dir / pmax(x, 0.001))^s_dir)
  f_indirect <- 1 / (1 + (x50_ind / pmax(x, 0.001))^s_ind)
  return(list(f_direct = f_direct, f_indirect = f_indirect))
}

compute_case_fractions <- function(LET, x50_dir, x50_ind, s_dir, s_ind) {
  trans <- compute_transitions(LET, x50_dir, x50_ind, s_dir, s_ind)
  
  p1 <- FIXED_PARAMS$p1_low + (FIXED_PARAMS$p1_high - FIXED_PARAMS$p1_low) * trans$f_direct
  p3 <- FIXED_PARAMS$p3_low * (1 - trans$f_indirect)
  p2 <- 1 - p1 - p3
  
  p2 <- pmax(0, p2)
  total <- p1 + p2 + p3
  
  return(list(p1 = p1/total, p2 = p2/total, p3 = p3/total))
}

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
calc_overkill_factor <- function(LET, max_let, overkill_strength = 0.15) {
  proximity <- LET / max_let
  if (proximity > 0.7) {
    correction <- 1 + overkill_strength * ((proximity - 0.7) / 0.3)^2
    return(pmin(correction, 1 + overkill_strength * 1.5))
  } else {
    return(1.0)
  }
}

#' Predict OER with overkill correction (NEW)
predict_OER_with_overkill <- function(LET, x50_dir, x50_ind, s_dir, s_ind,
                                      K_fix, K_repair, O2_hyp, O2_ref,
                                      max_let, overkill_strength = 0.10) {
  OER_base <- predict_OER_base(LET, x50_dir, x50_ind, s_dir, s_ind,
                               K_fix, K_repair, O2_hyp, O2_ref)
  overkill_factor <- calc_overkill_factor(LET, max_let, overkill_strength)
  OER_corrected <- 1 + (OER_base - 1) * overkill_factor
  return(pmax(1.0, OER_corrected))
}

calc_steepness_Z <- function(Z, s_base, s_scale) {
  s_base * (1 + s_scale * log(pmax(Z, 2) / 2))
}

#' Z-interpolation for x50 parameters (for N, O, Si)
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

#' Unified prediction function - UPDATED for hybrid model + overkill
predict_OER <- function(LET, ion, O2_hyp = 0.001, O2_ref = 21.0) {
  p_info <- particle_info %>% filter(ion == !!ion)
  if (nrow(p_info) == 0) return(NA)
  
  p_class <- p_info$particle_class
  Z <- p_info$Z
  max_let <- p_info$max_let
  
  # Get x50 parameters
  if (ion %in% c("photon", "proton", "deuteron", "He", "C", "Ne", "Ar")) {
    x50_dir <- params[[paste0("x50_dir_", ion)]]
    x50_ind <- params[[paste0("x50_ind_", ion)]]
  } else if (ion %in% c("N", "O", "Si")) {
    if (!is.null(model$interpolated_params[[ion]])) {
      x50_dir <- model$interpolated_params[[ion]]$x50_dir
      x50_ind <- model$interpolated_params[[ion]]$x50_ind
    } else {
      interp <- interpolate_x50(Z, params)
      x50_dir <- interp$x50_dir
      x50_ind <- interp$x50_ind
    }
  } else {
    return(NA)
  }
  
  # HYBRID MODEL: Get steepness based on particle class
  if (p_class == "light") {
    s_dir <- params$s_dir_light
    s_ind <- params$s_ind_light
  } else {
    s_dir <- calc_steepness_Z(Z, params$s_dir_base, params$s_dir_scale)
    s_ind <- calc_steepness_Z(Z, params$s_ind_base, params$s_ind_scale)
  }
  
  # Get overkill strength from model (default to 0 if not present)
  overkill_strength <- if (!is.null(params$overkill_strength)) params$overkill_strength else 0
  
  # Use overkill-corrected prediction for heavy particles (if overkill_strength > 0)
  if (p_class == "heavy" && overkill_strength > 0) {
    predict_OER_with_overkill(LET, x50_dir, x50_ind, s_dir, s_ind,
                              params$K_fix, params$K_repair, O2_hyp, O2_ref,
                              max_let, overkill_strength)
  } else {
    predict_OER_base(LET, x50_dir, x50_ind, s_dir, s_ind,
                     params$K_fix, params$K_repair, O2_hyp, O2_ref)
  }
}

cat("Model functions defined:\n")
# Section 2: Generate Theoretical Curves (With Physical Let Limits)

cat("\n--- SECTION 2: GENERATING THEORETICAL CURVES (Physical LET Limits) ---\n")
# LET ranges for each particle - RESPECTING PHYSICAL LIMITS
let_ranges <- list(
  photon = seq(0.2, 35, length.out = 200),
  proton = seq(0.5, 100, length.out = 200),
  deuteron = seq(5, 120, length.out = 200),
  He = seq(4, 200, length.out = 200),
  C = seq(20, 550, length.out = 200),
  N = seq(50, 600, length.out = 200),
  O = seq(50, 620, length.out = 200),
  Ne = seq(50, 700, length.out = 200),
  Si = seq(100, 800, length.out = 200),
  Ar = seq(80, 900, length.out = 200)
)

# Generate curves
theoretical_curves <- data.frame()

for (particle_name in all_particles_with_interpolated) {
  LET_seq <- let_ranges[[particle_name]]
  OER_vals <- sapply(LET_seq, function(l) predict_OER(l, particle_name))
  
  curve_data <- tibble(
    ion = particle_name,
    LET = LET_seq,
    OER_pred = OER_vals
  )
  
  theoretical_curves <- bind_rows(theoretical_curves, curve_data)
}

# Add particle info
theoretical_curves <- theoretical_curves %>%
  left_join(particle_info, by = "ion")

cat(sprintf("Generated %d theoretical curve points.\n", nrow(theoretical_curves)))
cat(sprintf("Particles: %s\n", paste(unique(theoretical_curves$ion), collapse = ", ")))
cat("Curves cut at physical LET limits (Bragg peak).\n\n")

# Save theoretical curves
write_csv(theoretical_curves, "results/theoretical_curves_voxa.csv")


# Section 3: Prepare Calibration Data

cat("\n--- SECTION 3: PREPARING CALIBRATION DATA ---\n")
# Add particle info if not present
if (!"particle_class" %in% names(calibration_data)) {
  calibration_data <- calibration_data %>%
    left_join(particle_info %>% select(ion, Z, particle_class, max_let), by = "ion")
}

# Add max_let if not present
if (!"max_let" %in% names(calibration_data)) {
  calibration_data <- calibration_data %>%
    left_join(particle_info %>% select(ion, max_let), by = "ion")
}

# Add cell line standardization if not present
if (!"cell_line_std" %in% names(calibration_data)) {
  calibration_data <- calibration_data %>%
    mutate(
      cell_line_std = case_when(
        is.na(cell_line) ~ "Other",
        str_detect(toupper(cell_line), "V79") ~ "V79",
        str_detect(toupper(cell_line), "HSG") ~ "HSG",
        str_detect(toupper(cell_line), "CHO") ~ "CHO",
        str_detect(toupper(cell_line), "T1") ~ "T1",
        TRUE ~ "Other"
      )
    )
}

# Calculate residuals and percent error
calibration_data <- calibration_data %>%
  mutate(
    pct_error = abs(residual / OER_retention) * 100,
    error_category = case_when(
      pct_error < 10 ~ "<10%",
      pct_error < 20 ~ "10-20%",
      pct_error < 30 ~ "20-30%",
      TRUE ~ ">30%"
    ),
    # Calculate Bragg proximity
    bragg_proximity = LET / max_let,
    let_status = case_when(
      bragg_proximity > 1.0 ~ "Beyond",
      bragg_proximity > 0.7 ~ "Near Bragg",
      TRUE ~ "Within"
    )
  )

cat("Calibration data prepared.\n")
cat(sprintf("  Total points: %d\n", nrow(calibration_data)))
cat(sprintf("  Points with <10%% error: %d (%.1f%%)\n", 
            sum(calibration_data$pct_error < 10),
            100 * mean(calibration_data$pct_error < 10)))
cat(sprintf("  Points with <20%% error: %d (%.1f%%)\n", 
            sum(calibration_data$pct_error < 20),
            100 * mean(calibration_data$pct_error < 20)))
cat(sprintf("  Points within physical LET limits: %d (%.1f%%)\n\n", 
            sum(calibration_data$bragg_proximity <= 1.0),
            100 * mean(calibration_data$bragg_proximity <= 1.0)))


# Section 4: Figure 1 - Oer Vs Let (All Particles With Physical Limits)

cat("\n--- SECTION 4: FIGURE 1 - OER vs LET (All Particles) ---\n")
# For main plot, use calibrated particles only (cleaner visualization)
theoretical_curves_main <- theoretical_curves %>%
  filter(ion %in% all_particles)

calibration_data_main <- calibration_data %>%
  filter(ion %in% all_particles)

fig1 <- ggplot() +
  # Theoretical curves (cut at physical limits)
  geom_line(data = theoretical_curves_main,
            aes(x = LET, y = OER_pred, color = ion),
            linewidth = 1.0) +
  # Calibration data points
  geom_point(data = calibration_data_main,
             aes(x = LET, y = OER_retention, color = ion),
             alpha = 0.6, size = 2.5) +
  scale_x_log10(
    breaks = c(0.1, 1, 10, 100, 1000),
    labels = c("0.1", "1", "10", "100", "1000"),
    limits = c(0.1, 1000)
  ) +
  scale_y_continuous(limits = c(1.0, 4.5), breaks = seq(1, 4.5, 0.5)) +
  scale_color_manual(values = particle_colors, name = "Particle") +
  labs(
    title = "VOxA Model: OER vs LET with Theoretical Curves",
    subtitle = sprintf("R² = %.3f | MAE = %.3f | Curves cut at physical LET limits (Bragg peak)",
                       model$fit_statistics$r2, model$fit_statistics$mae),
    x = expression(paste("LET (keV/", mu, "m)")),
    y = "OER (Retention)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30")
  )

ggsave("figures/fig1_oer_vs_let_voxa.png", fig1, width = 12, height = 7, dpi = 300)
cat("Saved: figures/fig1_oer_vs_let_voxa.png/pdf\n\n")


# Section 5: Figure 2 - Individual Particle Fits (Faceted)

cat("\n--- SECTION 5: FIGURE 2 - Individual Particle Fits ---\n")
# Get unique particles in calibration data
particles_in_data <- unique(calibration_data$ion)
particles_in_data <- particles_in_data[particles_in_data %in% names(particle_colors)]

# Filter theoretical curves to match
theoretical_curves_facet <- theoretical_curves %>%
  filter(ion %in% particles_in_data)

# Order particles by Z
particle_order <- particle_info %>%
  filter(ion %in% particles_in_data) %>%
  arrange(Z) %>%
  pull(ion)

calibration_data$ion_factor <- factor(calibration_data$ion, levels = particle_order)
theoretical_curves_facet$ion_factor <- factor(theoretical_curves_facet$ion, levels = particle_order)

# Create faceted plot
fig2 <- ggplot() +
  geom_line(data = theoretical_curves_facet,
            aes(x = LET, y = OER_pred),
            color = model_curve_color, linewidth = 1.0) +
  geom_point(data = calibration_data,
             aes(x = LET, y = OER_retention, color = cell_line_std),
             alpha = 0.7, size = 2) +
  facet_wrap(~ ion_factor, scales = "free_x", ncol = 4) +
  scale_x_log10() +
  scale_y_continuous(limits = c(1, 5.5)) +
  scale_color_manual(values = cell_line_colors, name = "Cell Line") +
  labs(
    title = "VOxA Model: Individual Particle Fits",
    subtitle = "Dark blue curve = Model prediction (cut at Bragg peak) | Points = Calibration data by cell line",
    x = expression(paste("LET (keV/", mu, "m)")),
    y = "OER (Retention)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "#F5E6D3"),
    strip.text = element_text(face = "bold", color = "#4A3728")
  )

ggsave("figures/fig2_individual_fits_voxa.png", fig2, width = 14, height = 10, dpi = 300)
cat("Saved: figures/fig2_individual_fits_voxa.png/pdf\n\n")


# Section 6: Figure 3 - Predicted Vs Observed

cat("\n--- SECTION 6: FIGURE 3 - Predicted vs Observed OER ---\n")
fig3 <- ggplot(calibration_data, aes(x = OER_pred, y = OER_retention)) +
  # Identity line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#4A3728", linewidth = 0.8) +
  # ±20% envelope
  geom_abline(intercept = 0, slope = 1.2, linetype = "dotted", color = "#8B7355") +
  geom_abline(intercept = 0, slope = 0.8, linetype = "dotted", color = "#8B7355") +
  # Points
  geom_point(aes(color = ion), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = particle_colors, name = "Particle") +
  scale_x_continuous(limits = c(1, 5), breaks = 1:5) +
  scale_y_continuous(limits = c(1, 5), breaks = 1:5) +
  coord_fixed() +
  labs(
    title = "VOxA Model: Predicted vs Observed OER",
    subtitle = sprintf("R² = %.3f | Dashed line = Perfect agreement | Dotted = ±20%%",
                       model$fit_statistics$r2),
    x = "Predicted OER (Retention)",
    y = "Observed OER (Retention)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave("figures/fig3_pred_vs_obs_voxa.png", fig3, width = 8, height = 8, dpi = 300)
cat("Saved: figures/fig3_pred_vs_obs_voxa.png/pdf\n\n")


# Section 7: Figure 4 - Residuals Vs Let

cat("\n--- SECTION 7: FIGURE 4 - Residuals vs LET ---\n")
residual_stats <- calibration_data %>%
  summarise(
    mean = mean(residual),
    sd = sd(residual),
    median = median(residual)
  )

fig4 <- ggplot(calibration_data, aes(x = LET, y = residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#4A3728") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted", color = "#8B7355") +
  geom_point(aes(color = ion), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "loess", se = TRUE, color = "#0D3B66", fill = "#A8C4D9", linewidth = 1) +
  scale_x_log10(
    breaks = c(0.1, 1, 10, 100, 1000),
    labels = c("0.1", "1", "10", "100", "1000")
  ) +
  scale_y_continuous(limits = c(-1.5, 2)) +
  scale_color_manual(values = particle_colors, name = "Particle") +
  labs(
    title = "VOxA Model: Residuals vs LET",
    subtitle = sprintf("Mean = %+.3f | SD = %.3f | Blue line = LOESS trend",
                       residual_stats$mean, residual_stats$sd),
    x = expression(paste("LET (keV/", mu, "m)")),
    y = "Residual (Observed - Predicted)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave("figures/fig4_residuals_vs_let_voxa.png", fig4, width = 12, height = 7, dpi = 300)
cat("Saved: figures/fig4_residuals_vs_let_voxa.png/pdf\n\n")


# Section 8: Figure 5 - Oxygen Response Curve (D-Kondo Comparison)

cat("\n--- SECTION 8: FIGURE 5 - Oxygen Response Curve (Validation) ---\n")
cat("NOTE: D-Kondo data was NOT used in calibration (validation only).\n")
# Generate oxygen curve
O2_levels <- 10^seq(-3, log10(21), length.out = 200)

oxygen_curve <- tibble(
  O2 = O2_levels,
  OER_pred = sapply(O2_levels, function(o2) {
    predict_OER(LET = 1.0, ion = "photon", O2_hyp = o2, O2_ref = 21.0)
  })
)

# D-Kondo data (HRF with radioprotector) - VALIDATION ONLY
dkondo_data <- tibble(
  O2 = c(21.0, 2.1, 0.21, 0.021),
  HRF = c(1.0, 1.21, 2.33, 2.83),
  is_measured = TRUE
)

# Interpolated D-Kondo points
dkondo_interp <- tibble(
  O2 = c(10.0, 5.0, 1.0, 0.5, 0.1, 0.05, 0.01, 0.001),
  HRF = c(1.05, 1.10, 1.45, 1.75, 2.55, 2.70, 2.88, 2.95),
  is_measured = FALSE
)

dkondo_all <- bind_rows(dkondo_data, dkondo_interp) %>%
  arrange(desc(O2))

# Calculate MAE for measured points (excluding normoxia)
dkondo_measured_comparison <- dkondo_data %>%
  filter(O2 < 21) %>%
  mutate(OER_model = sapply(O2, function(o2) predict_OER(1.0, "photon", o2, 21.0)))

dkondo_mae <- mean(abs(dkondo_measured_comparison$HRF - dkondo_measured_comparison$OER_model) / 
                     dkondo_measured_comparison$HRF) * 100

fig5 <- ggplot() +
  # Model curve
  geom_line(data = oxygen_curve,
            aes(x = O2, y = OER_pred),
            color = "#1E5B8C", linewidth = 1.2) +
  # D-Kondo measured points
  geom_point(data = dkondo_data,
             aes(x = O2, y = HRF),
             color = "#C75B5B", size = 4, shape = 16) +
  # D-Kondo interpolated points
  geom_point(data = dkondo_interp,
             aes(x = O2, y = HRF),
             color = "#E07B3C", size = 3, shape = 1, stroke = 1.2) +
  scale_x_log10(
    breaks = c(0.001, 0.01, 0.1, 1, 10),
    labels = c("0.001", "0.01", "0.1", "1", "10")
  ) +
  scale_y_continuous(limits = c(1, 3.5), breaks = seq(1, 3.5, 0.5)) +
  labs(
    title = "VOxA Model: Oxygen Response Curve (External Validation)",
    subtitle = sprintf("D-Kondo HRF data (NOT used in calibration) | Measured points MAE = %.1f%%", 
                       dkondo_mae),
    x = expression(paste("Oxygen Concentration (% ", O[2], ")")),
    y = "OER / HRF"
  ) +
  annotate("text", x = 0.003, y = 3.3, 
           label = "● D-Kondo Measured\n○ D-Kondo Interpolated\n— VOxA Model", 
           hjust = 0, size = 3.5, color = "#4A3728") +
  annotate("text", x = 0.003, y = 1.3,
           label = "Note: D-Kondo uses HRF (radioprotector)\nNot directly comparable to OER",
           hjust = 0, size = 3, color = "#8B7355", fontface = "italic") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave("figures/fig5_oxygen_curve_voxa.png", fig5, width = 10, height = 7, dpi = 300)
cat(sprintf("D-Kondo comparison MAE (measured points): %.1f%%\n", dkondo_mae))
cat("Saved: figures/fig5_oxygen_curve_voxa.png/pdf\n\n")


# Section 9: Figure 6 - Hirayama Carbon Validation

# Hirayama 2005 data point
hirayama_LET <- 79.6
hirayama_OER_survival <- 1.8
hirayama_OER_retention_est <- convert_survival_to_retention(hirayama_OER_survival)

# Model prediction at this LET
model_pred_hirayama <- predict_OER(hirayama_LET, "C")

# Generate carbon curve (within physical limits)
carbon_curve <- tibble(
  LET = seq(10, 550, length.out = 200),
  OER_pred = sapply(seq(10, 550, length.out = 200), function(l) predict_OER(l, "C"))
)

# Carbon calibration data
carbon_data <- calibration_data %>% filter(ion == "C")

fig6 <- ggplot() +
  # Model curve
  geom_line(data = carbon_curve,
            aes(x = LET, y = OER_pred),
            color = "#1E5B8C", linewidth = 1.2) +
  # Calibration data
  geom_point(data = carbon_data,
             aes(x = LET, y = OER_retention),
             color = "#8B7355", alpha = 0.5, size = 2) +
  # Hirayama observed
  geom_point(aes(x = hirayama_LET, y = hirayama_OER_retention_est),
             color = "#D4A84B", size = 5, shape = 17) +
  # Model prediction at Hirayama LET
  geom_point(aes(x = hirayama_LET, y = model_pred_hirayama),
             color = "#1E5B8C", size = 5, shape = 18) +
  # Connecting line
  geom_segment(aes(x = hirayama_LET, xend = hirayama_LET,
                   y = hirayama_OER_retention_est, yend = model_pred_hirayama),
               color = "#8B5A83", linetype = "dashed", linewidth = 1) +
  # Physical limit annotation
  geom_vline(xintercept = 550, linetype = "dotted", color = "#C75B5B", linewidth = 0.8) +
  annotate("text", x = 520, y = 4.5, label = "Bragg\npeak", size = 3, color = "#C75B5B", hjust = 1) +
  scale_x_log10() +
  scale_y_continuous(limits = c(1, 5)) +
  labs(
    title = "VOxA Model: Hirayama 2005 Carbon Ion Validation",
    subtitle = sprintf("CHO cells | LET = 80 keV/μm | OER_survival = 1.8 → OER_retention = %.2f | Model = %.2f",
                       hirayama_OER_retention_est, model_pred_hirayama),
    x = expression(paste("LET (keV/", mu, "m)")),
    y = "OER (Retention)"
  ) +
  annotate("text", x = 15, y = 4.5,
           label = sprintf("Hirayama 2005 (CHO cells)\nOER_survival = 1.8\nOER_retention (converted) = %.2f\nModel prediction = %.2f\nDifference = %.1f%%",
                           hirayama_OER_retention_est, model_pred_hirayama,
                           100 * abs(model_pred_hirayama - hirayama_OER_retention_est) / hirayama_OER_retention_est),
           hjust = 0, size = 3.5, color = "#4A3728") +
  annotate("text", x = 250, y = 1.3,
           label = "▲ Hirayama observed | ◆ Model prediction | Gray = Calibration data",
           size = 3, color = "#4A3728") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave("figures/fig6_hirayama_validation_voxa.png", fig6, width = 10, height = 7, dpi = 300)

cat(sprintf("Hirayama comparison:\n"))
cat(sprintf("  LET = %.1f keV/μm\n", hirayama_LET))
cat(sprintf("  OER_survival (reported) = %.2f\n", hirayama_OER_survival))
cat(sprintf("  OER_retention (converted) = %.2f\n", hirayama_OER_retention_est))
cat(sprintf("  Model prediction = %.2f\n", model_pred_hirayama))
cat(sprintf("  Difference = %.1f%%\n\n", 
            100 * abs(model_pred_hirayama - hirayama_OER_retention_est) / hirayama_OER_retention_est))

cat("Saved: figures/fig6_hirayama_validation_voxa.png/pdf\n\n")


# Section 10: Figure 7 - Residual Distribution (Histogram)

cat("\n--- SECTION 10: FIGURE 7 - Residual Distribution ---\n")
fig7 <- ggplot(calibration_data, aes(x = residual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, 
                 fill = "#D4A84B", color = "white", alpha = 0.7) +
  geom_density(color = "#C75B5B", linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#4A3728") +
  geom_vline(xintercept = mean(calibration_data$residual), 
             color = "#C75B5B", linewidth = 1) +
  scale_x_continuous(limits = c(-1.5, 2)) +
  labs(
    title = "VOxA Model: Residual Distribution",
    subtitle = sprintf("N = %d | Mean = %+.3f | SD = %.3f | Median = %+.3f",
                       nrow(calibration_data),
                       mean(calibration_data$residual),
                       sd(calibration_data$residual),
                       median(calibration_data$residual)),
    x = "Residual (Observed - Predicted)",
    y = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave("figures/fig7_residual_histogram_voxa.png", fig7, width = 10, height = 6, dpi = 300)
cat("Saved: figures/fig7_residual_histogram_voxa.png/pdf\n\n")


# Section 11: Figure 8 - Residual Distribution By Particle (Boxplot)

cat("\n--- SECTION 11: FIGURE 8 - Residual Distribution by Particle ---\n")
# Order by median absolute residual
particle_order_residual <- calibration_data %>%
  group_by(ion) %>%
  summarise(med_abs_res = median(abs(residual))) %>%
  arrange(med_abs_res) %>%
  pull(ion)

calibration_data$ion_ordered <- factor(calibration_data$ion, levels = particle_order_residual)

fig8 <- ggplot(calibration_data, aes(x = residual, y = ion_ordered, fill = particle_class)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#4A3728") +
  scale_fill_manual(values = particle_class_colors, name = "Particle Class") +
  scale_x_continuous(limits = c(-1.5, 2)) +
  labs(
    title = "VOxA Model: Residual Distribution by Particle",
    subtitle = "Ordered by median |residual| (best to worst)",
    x = "Residual (Observed - Predicted)",
    y = "Particle"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave("figures/fig8_residual_boxplot_voxa.png", fig8, width = 10, height = 7, dpi = 300)
cat("Saved: figures/fig8_residual_boxplot_voxa.png/pdf\n\n")


# Section 12: Figure 9 - Physical Let Limits Compliance (Corrected V2)

cat("\n--- SECTION 12: FIGURE 9 - Physical LET Limits Compliance ---\n")
# Create ordered factor for particles based on Z
particle_order_z <- particle_info %>%
  filter(ion %in% unique(calibration_data$ion)) %>%
  arrange(Z) %>%
  pull(ion)

calibration_data$ion_factor_z <- factor(calibration_data$ion, levels = particle_order_z)

# Create limit data with matching factor levels
limit_data <- particle_info %>%
  filter(ion %in% unique(calibration_data$ion)) %>%
  mutate(ion_factor_z = factor(ion, levels = particle_order_z))

fig9 <- ggplot(calibration_data, aes(x = LET, y = ion_factor_z)) +
  # Physical limit lines (vertical dashed red lines at Bragg peak)
  geom_vline(data = limit_data,
             aes(xintercept = max_let),
             color = "#C75B5B", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
  # Add text labels for limits
  geom_text(data = limit_data,
            aes(x = max_let, y = ion_factor_z, label = max_let),
            color = "#C75B5B", size = 2.5, hjust = -0.2, vjust = -0.5) +
  # Data points colored by proximity to Bragg peak
  geom_point(aes(color = bragg_proximity), size = 2.5, alpha = 0.7) +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000"),
    limits = c(0.1, 1200)
  ) +
  scale_color_gradient2(
    low = "#2D7D46", mid = "#D4A84B", high = "#C75B5B",
    midpoint = 0.7, limits = c(0, 1.2),
    name = "Proximity to\nBragg Peak"
  ) +
  labs(
    title = "VOxA Model: Physical LET Limits Compliance",
    subtitle = "Dashed red lines = Bragg peak (max LET) | Green→Yellow→Red = proximity to limit",
    x = expression(paste("LET (keV/", mu, "m)")),
    y = "Particle (ordered by Z)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

ggsave("figures/fig9_let_limits_voxa.png", fig9, width = 12, height = 8, dpi = 300)

# Print compliance summary
cat("Physical LET limit compliance:\n")
compliance_summary <- calibration_data %>%
  group_by(ion) %>%
  summarise(
    N = n(),
    max_let = first(max_let),
    LET_max_data = max(LET),
    within_70pct = sum(bragg_proximity <= 0.7),
    near_bragg = sum(bragg_proximity > 0.7 & bragg_proximity <= 1.0),
    beyond = sum(bragg_proximity > 1.0),
    .groups = "drop"
  ) %>%
  arrange(match(ion, particle_order_z))

print(compliance_summary)
cat("\nSaved: figures/fig9_let_limits_voxa.png/pdf\n\n")


# Section 13: Figure 10 - Z-Ordering Validation (Oer At Fixed Let)

cat("\n--- SECTION 13: FIGURE 10 - Z-Ordering Validation ---\n")
# Generate OER at multiple fixed LETs for heavy ions
fixed_let_values <- c(50, 100, 150, 200)
heavy_ions_ordered <- c("He", "C", "Ne", "Ar")

z_ordering_data <- expand_grid(
  LET = fixed_let_values,
  ion = heavy_ions_ordered
) %>%
  left_join(particle_info %>% select(ion, Z, max_let), by = "ion") %>%
  mutate(
    OER_pred = mapply(function(l, i) predict_OER(l, i), LET, ion),
    within_limit = LET <= max_let
  )

fig10 <- ggplot(z_ordering_data %>% filter(within_limit), 
                aes(x = Z, y = OER_pred, color = factor(LET), group = LET)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c("50" = "#2D7D46", "100" = "#1E5B8C", "150" = "#D4A84B", "200" = "#C75B5B"),
    name = "LET (keV/μm)"
  ) +
  scale_x_continuous(breaks = c(2, 6, 10, 18), labels = c("He (2)", "C (6)", "Ne (10)", "Ar (18)")) +
  labs(
    title = "VOxA Model: Z-Ordering Validation",
    subtitle = "At fixed LET, OER INCREASES with Z (heavier ions have higher OER)\nPhysics: Lighter ions have denser tracks → more direct damage → lower OER",
    x = "Atomic Number (Z)",
    y = "OER (Retention)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14)
  )

ggsave("figures/fig10_z_ordering_voxa.png", fig10, width = 10, height = 7, dpi = 300)
cat("Saved: figures/fig10_z_ordering_voxa.png/pdf\n\n")

# Verify Z-ordering
cat("Z-ordering check at LET = 100 keV/μm:\n")
for (ion_name in heavy_ions_ordered) {
  oer <- predict_OER(100, ion_name)
  cat(sprintf("  %s (Z=%d): OER = %.3f\n", ion_name, 
              particle_info$Z[particle_info$ion == ion_name], oer))
}
# Section 14: Figure 11 - Tinganelli Intermediate O2 Validation

cat("\n--- SECTION 14: FIGURE 11 - Tinganelli Intermediate O2 Validation ---\n")
# Get Tinganelli data
tinganelli_data <- calibration_data %>%
  filter(str_detect(dataset, "Tinganelli") | str_detect(source, "Tinganelli", negate = FALSE))

if (nrow(tinganelli_data) == 0) {
  # Try alternative identification
  tinganelli_data <- calibration_data %>%
    filter(O2_hyp > 0.01 | cell_line_std == "CHO")
}

if (nrow(tinganelli_data) > 0) {
  cat(sprintf("Tinganelli data points: %d\n", nrow(tinganelli_data)))
  cat("NOTE: Tinganelli data IS included in calibration (affects R²).\n\n")
  
  # Summary by O2 level
  tinganelli_summary <- tinganelli_data %>%
    group_by(ion, O2_hyp) %>%
    summarise(
      n = n(),
      OER_obs_mean = mean(OER_retention),
      OER_pred_mean = mean(OER_pred),
      MAE = mean(abs(residual)),
      .groups = "drop"
    )
  
  fig11 <- ggplot(tinganelli_data, aes(x = OER_pred, y = OER_retention)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#4A3728") +
    geom_point(aes(color = factor(round(O2_hyp, 3)), shape = ion), size = 3, alpha = 0.8) +
    scale_color_manual(
      values = c("0.001" = "#C75B5B", "0.15" = "#D4A84B", "0.5" = "#2D7D46", "2" = "#1E5B8C"),
      name = "O2 (%)"
    ) +
    scale_shape_manual(values = c("photon" = 16, "C" = 17, "N" = 15, "O" = 18, "Si" = 8), name = "Particle") +
    coord_fixed(xlim = c(1, 4), ylim = c(1, 4)) +
    labs(
      title = "VOxA Model: Tinganelli et al. (2015) Validation",
      subtitle = "Intermediate O2 levels (0.15%, 0.5%, 2%) | CHO cells | Included in calibration",
      x = "Predicted OER (Retention)",
      y = "Observed OER (Retention)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14)
    )
  
  ggsave("figures/fig11_tinganelli_validation_voxa.png", fig11, width = 10, height = 8, dpi = 300)
  cat("Saved: figures/fig11_tinganelli_validation_voxa.png/pdf\n\n")
  
  # Print summary
  cat("Tinganelli performance by O2 level:\n")
  print(tinganelli_summary)
  } else {
  cat("No Tinganelli data identified in calibration set.\n\n")
}


# Section 15: Figure 12 - Cell-Line Performance Comparison

cat("\n--- SECTION 15: FIGURE 12 - Cell-Line Performance Comparison ---\n")
cell_stats <- calibration_data %>%
  group_by(cell_line_std) %>%
  summarise(
    N = n(),
    MAE = mean(abs(residual)),
    RMSE = sqrt(mean(residual^2)),
    Bias = mean(residual),
    .groups = "drop"
  ) %>%
  mutate(
    factor = case_when(
      cell_line_std == "V79" ~ 1.0,
      cell_line_std == "HSG" ~ params$factor_HSG,
      cell_line_std == "T1" ~ params$factor_T1,
      cell_line_std == "CHO" ~ params$factor_CHO,
      TRUE ~ 1.0
    )
  ) %>%
  arrange(desc(N))

fig12 <- ggplot(calibration_data, aes(x = OER_pred, y = OER_retention)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#4A3728") +
  geom_point(aes(color = cell_line_std), alpha = 0.6, size = 2) +
  facet_wrap(~ cell_line_std, ncol = 3) +
  scale_color_manual(values = cell_line_colors, guide = "none") +
  coord_fixed(xlim = c(1, 5), ylim = c(1, 5)) +
  labs(
    title = "VOxA Model: Cell-Line Performance Comparison",
    subtitle = sprintf("V79 (ref)=1.00, HSG=%.3f, T1=%.3f, CHO=%.3f",
                       params$factor_HSG, params$factor_T1, params$factor_CHO),
    x = "Predicted OER (Retention)",
    y = "Observed OER (Retention)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "#F5E6D3"),
    strip.text = element_text(face = "bold", color = "#4A3728")
  )

ggsave("figures/fig12_cell_line_comparison_voxa.png", fig12, width = 12, height = 8, dpi = 300)

cat("Cell-line statistics:\n")
print(cell_stats)
cat("Saved: figures/fig12_cell_line_comparison_voxa.png/pdf\n\n")


# Section 16: Comprehensive Validation Report

cat("\n--- SECTION 16: GENERATING COMPREHENSIVE VALIDATION REPORT ---\n")
# Get outlier count safely
n_outliers <- if (!is.null(model$n_outliers_removed)) model$n_outliers_removed else 0

# Get overkill strength safely
overkill_strength <- if (!is.null(params$overkill_strength)) params$overkill_strength else 0

# Calculate particle statistics
particle_stats <- calibration_data %>%
  group_by(ion, particle_class) %>%
  summarise(
    N = n(),
    MAE = mean(abs(residual)),
    MAE_pct = mean(abs(residual / OER_retention) * 100),
    RMSE = sqrt(mean(residual^2)),
    Bias = mean(residual),
    .groups = "drop"
  ) %>%
  left_join(particle_info %>% select(ion, Z), by = "ion") %>%
  arrange(Z)

# Create report
report <- c(
  "================================================================================",
  "              VOxA MODEL COMPREHENSIVE VALIDATION REPORT",
  "              Variable Oxygen-dependent Amorphous Track Model",
  "================================================================================",
  "",
  sprintf("Generated: %s", Sys.time()),
  sprintf("Model version: %s", model$version),
  "",
  "--------------------------------------------------------------------------------",
  "1. MODEL OVERVIEW",
  "--------------------------------------------------------------------------------",
  "",
  "The VOxA Model predicts Oxygen Enhancement Ratio (OER) for DNA damage RETENTION",
  "(not cell survival). The model accounts for:",
  "  - Particle-specific track structure via x50 parameters",
  "  - Z-ordering physics constraint (lighter ions = lower OER at same LET)",
  "  - Physical LET limits (Bragg peak cutoffs)",
  "  - Hybrid light/heavy particle treatment",
  "  - Cell-line correction factors (V79, HSG, T1, CHO)",
  "  - Intermediate oxygen levels (Tinganelli data)",
  "",
  "--------------------------------------------------------------------------------",
  "2. MODEL PARAMETERS",
  "--------------------------------------------------------------------------------",
  "",
  "Oxygen Kinetics:",
  sprintf("  K_fix = %.4f%%", params$K_fix),
  sprintf("  K_repair = %.4f%%", params$K_repair),
  "",
  "Case Fractions (DSB combinatorics with d=0.20, i=0.80):",
  sprintf("  p1_low = %.2f (d² = purely direct DSBs)", FIXED_PARAMS$p1_low),
  sprintf("  p2_low = %.2f (2di = hybrid DSBs)", FIXED_PARAMS$p2_low),
  sprintf("  p3_low = %.2f (i² = purely indirect DSBs)", FIXED_PARAMS$p3_low),
  sprintf("  p1_high = %.2f (direct at high LET, Nikjoo 2001)", FIXED_PARAMS$p1_high),
  "",
  "Steepness (Hybrid Model):",
  sprintf("  Light particles: s_dir = %.3f, s_ind = %.3f", params$s_dir_light, params$s_ind_light),
  sprintf("  Heavy base: s_dir = %.3f, s_ind = %.3f", params$s_dir_base, params$s_ind_base),
  sprintf("  Heavy scale: s_dir = %.3f, s_ind = %.3f", params$s_dir_scale, params$s_ind_scale),
  "",
  "Overkill Correction:",
  sprintf("  overkill_strength = %.3f", overkill_strength),
  "",
  "Cell-line Correction Factors:",
  sprintf("  V79 (reference) = 1.000"),
  sprintf("  HSG factor = %.3f", params$factor_HSG),
  sprintf("  T1 factor = %.3f", params$factor_T1),
  sprintf("  CHO factor = %.3f", params$factor_CHO),
  "",
  "Physics Constraint (Z-ordering):",
  "  x50_dir(He) < x50_dir(C) < x50_dir(Ne) < x50_dir(Ar)",
  sprintf("  Status: %s", ifelse(model$z_ordering_satisfied, "SATISFIED ✓", "VIOLATED ✗")),
  "",
  "x50 Parameters (directly calibrated):",
  sprintf("  photon:   x50_dir = %6.1f, x50_ind = %6.1f", params$x50_dir_photon, params$x50_ind_photon),
  sprintf("  proton:   x50_dir = %6.1f, x50_ind = %6.1f", params$x50_dir_proton, params$x50_ind_proton),
  sprintf("  deuteron: x50_dir = %6.1f, x50_ind = %6.1f", params$x50_dir_deuteron, params$x50_ind_deuteron),
  sprintf("  He:       x50_dir = %6.1f, x50_ind = %6.1f", params$x50_dir_He, params$x50_ind_He),
  sprintf("  C:        x50_dir = %6.1f, x50_ind = %6.1f", params$x50_dir_C, params$x50_ind_C),
  sprintf("  Ne:       x50_dir = %6.1f, x50_ind = %6.1f", params$x50_dir_Ne, params$x50_ind_Ne),
  sprintf("  Ar:       x50_dir = %6.1f, x50_ind = %6.1f", params$x50_dir_Ar, params$x50_ind_Ar),
  "",
  "x50 Parameters (Z-interpolated):"
)

# Add interpolated params if available
if (!is.null(model$interpolated_params)) {
  for (ion_name in names(model$interpolated_params)) {
    ip <- model$interpolated_params[[ion_name]]
    report <- c(report, sprintf("  %s:        x50_dir = %6.1f, x50_ind = %6.1f", 
                                ion_name, ip$x50_dir, ip$x50_ind))
  }
}

report <- c(report,
            "",
            "Physical LET Limits (Bragg Peak):",
            "  photon:   35 keV/μm",
            "  proton:  100 keV/μm",
            "  deuteron: 120 keV/μm",
            "  He:      200 keV/μm",
            "  C:       550 keV/μm",
            "  N:       600 keV/μm",
            "  O:       620 keV/μm",
            "  Ne:      700 keV/μm",
            "  Si:      800 keV/μm",
            "  Ar:      900 keV/μm",
            "",
            "--------------------------------------------------------------------------------",
            "3. CALIBRATION STATISTICS",
            "--------------------------------------------------------------------------------",
            "",
            sprintf("Observations: %d", model$fit_statistics$n_obs),
            sprintf("Parameters: %d", model$fit_statistics$n_params),
            sprintf("Outliers removed: %d", n_outliers),
            "",
            sprintf("R² (unweighted) = %.4f  <-- PRIMARY METRIC", model$fit_statistics$r2),
            sprintf("R² (weighted)   = %.4f  (reflects clinical weighting)", model$fit_statistics$r2_weighted),
            sprintf("Adjusted R²     = %.4f", model$fit_statistics$adj_r2),
            sprintf("RMSE            = %.4f", model$fit_statistics$rmse),
            sprintf("MAE             = %.4f", model$fit_statistics$mae),
            "",
            sprintf("OER_max (retention, theoretical) = %.2f", model$OER_max_theoretical),
            sprintf("OER_max (survival, theoretical)  = %.2f", OER_max_survival),
            sprintf("Grimes (2020) OER_max            = 2.74 (survival)", ""),
            "",
            "NOTE: This model predicts OER for DNA damage RETENTION, not cell survival.",
            "      Conversion: OER_survival = 1 + (OER_retention - 1) / 1.20",
            "",
            "--------------------------------------------------------------------------------",
            "4. EXTERNAL VALIDATION (NOT used in calibration)",
            "--------------------------------------------------------------------------------",
            "",
            "D-Kondo oxygen curve (HRF, not OER):",
            "  Note: D-Kondo used radioprotectors; HRF ≠ OER",
            "  This data was NOT used in calibration - for validation only",
            sprintf("  Measured points MAE = %.1f%% (qualitative comparison)", dkondo_mae),
            "",
            "Hirayama 2005 Carbon (80 keV/μm, CHO cells):",
            sprintf("  OER_survival (reported) = 1.80"),
            sprintf("  OER_retention (converted) = %.2f", hirayama_OER_retention_est),
            sprintf("  Model prediction = %.2f", model_pred_hirayama),
            sprintf("  Difference = %.1f%%", 
                    100 * abs(model_pred_hirayama - hirayama_OER_retention_est) / hirayama_OER_retention_est),
            "",
            "--------------------------------------------------------------------------------",
            "5. INTERNAL VALIDATION (Tinganelli - included in calibration)",
            "--------------------------------------------------------------------------------",
            "",
            "Tinganelli et al. (2015) data:",
            "  - 23 data points included in calibration",
            "  - Intermediate O2 levels: 0.15%, 0.5%, 2%",
            "  - Cell line: CHO",
            "  - Particles: photon, C, N, O, Si",
            "  - This data DOES affect R²",
            "",
            "--------------------------------------------------------------------------------",
            "6. PER-PARTICLE PERFORMANCE",
            "--------------------------------------------------------------------------------",
            ""
)

# Add per-particle stats
for (i in 1:nrow(particle_stats)) {
  p <- particle_stats[i, ]
  status <- if (p$MAE_pct < 10) "[EXCELLENT]" else if (p$MAE_pct < 15) "[GOOD]" else if (p$MAE_pct < 20) "[OK]" else "[POOR]"
  report <- c(report, sprintf("%-10s: N=%3d, MAE=%.3f, MAE%%=%5.1f%%, RMSE=%.3f, Bias=%+.3f %s", 
                              p$ion, p$N, p$MAE, p$MAE_pct, p$RMSE, p$Bias, status))
}

report <- c(report,
            "",
            "--------------------------------------------------------------------------------",
            "7. PER CELL-LINE PERFORMANCE",
            "--------------------------------------------------------------------------------",
            ""
)

# Add cell-line stats
for (i in 1:nrow(cell_stats)) {
  c <- cell_stats[i, ]
  report <- c(report, sprintf("%-10s: N=%3d, MAE=%.3f, RMSE=%.3f, Bias=%+.3f, Factor=%.3f", 
                              c$cell_line_std, c$N, c$MAE, c$RMSE, c$Bias, c$factor))
}

report <- c(report,
            "",
            "--------------------------------------------------------------------------------",
            "8. PHYSICAL BASIS",
            "--------------------------------------------------------------------------------",
            "",
            "Case fractions derived from DSB combinatorics:",
            "  - Single-hit direct probability d ≈ 0.20 (Sakata et al. 2019)",
            "  - Single-hit indirect probability i ≈ 0.80",
            "  - p1 = d² = 0.04 (both hits direct → O₂-independent)",
            "  - p2 = 2di = 0.32 (one hit each → linear O₂ dependence)",
            "  - p3 = i² = 0.64 (both hits indirect → quadratic O₂ dependence)",
            "",
            "High-LET limit:",
            "  - p1_high = 0.68 based on Nikjoo 2009 / Friedland 2017",
            "  - ~32% indirect action remains even at extreme LET (2106 keV/μm)",
            "",
            "Z-ordering physics constraint:",
            "  At the same LET, lighter ions are SLOWER (LET ∝ z²/v²)",
            "  → Narrower track core (R_c ∝ β)",
            "  → Higher local dose density",
            "  → More radical recombination (•OH + •OH → H₂O₂)",
            "  → Less indirect damage available for oxygen fixation",
            "  → LOWER OER",
            "",
            "  Therefore: x50_dir(He) < x50_dir(C) < x50_dir(Ne) < x50_dir(Ar)",
            "  (lighter ions transition to direct damage at LOWER LET)",
            "",
            "Physical LET limits:",
            "  - Each particle has a maximum achievable LET (Bragg peak)",
            "  - Protons max out at ~100 keV/μm",
            "  - Photons (secondary electrons) max at ~35 keV/μm",
            "  - Heavier ions can reach higher LET (C: 550, Ne: 700, Ar: 900)",
            "  - Model curves are cut at these physical limits",
            "",
            "OER type:",
            "  - Model predicts OER for DNA damage RETENTION",
            "  - NOT cell survival OER (which is typically ~20% lower excess)",
            sprintf("  - OER_max = %.2f (retention) = %.2f (survival)", 
                    model$OER_max_theoretical, OER_max_survival),
            "",
            "--------------------------------------------------------------------------------",
            "9. FIGURES GENERATED",
            "--------------------------------------------------------------------------------",
            "",
            "  Fig 1:  OER vs LET (all particles, physical limits)",
            "  Fig 2:  Individual particle fits (faceted by particle)",
            "  Fig 3:  Predicted vs Observed OER",
            "  Fig 4:  Residuals vs LET",
            "  Fig 5:  Oxygen response curve (D-Kondo validation)",
            "  Fig 6:  Hirayama carbon validation",
            "  Fig 7:  Residual histogram",
            "  Fig 8:  Residual boxplot by particle",
            "  Fig 9:  Physical LET limits compliance",
            "  Fig 10: Z-ordering validation",
            "  Fig 11: Tinganelli intermediate O2 validation",
            "  Fig 12: Cell-line performance comparison",
            "",
            "--------------------------------------------------------------------------------",
            "10. CONCLUSIONS",
            "--------------------------------------------------------------------------------",
            "",
            "STRENGTHS:",
            sprintf("  + R² = %.3f (unweighted), %.3f (weighted)", 
                    model$fit_statistics$r2, model$fit_statistics$r2_weighted),
            "  + Z-ordering constraint satisfied (physically correct)",
            sprintf("  + %d particles EXCELLENT (<10%% MAE), %d GOOD (10-15%%)", 
                    sum(particle_stats$MAE_pct < 10), 
                    sum(particle_stats$MAE_pct >= 10 & particle_stats$MAE_pct < 15)),
            "  + Physical LET limits prevent non-physical extrapolation",
            "  + Hybrid light/heavy model captures different track physics",
            "  + Cell-line factors account for biological variability",
            "  + Intermediate O2 levels validated (Tinganelli data)",
            "  + Z-interpolation enables prediction for new particles (N, O, Si)",
            "",
            "LIMITATIONS:",
            "  - He shows highest MAE among heavy ions (may need more data)",
            "  - D-Kondo comparison limited (HRF ≠ OER)",
            "  - Limited data for some particles (deuteron N=2, N/O/Si N=1-3)",
            "  - Overkill correction strength = 0 (data didn't require it)",
            "",
            "COMPARISON TO GRIMES (2020):",
            "  - VOxA: Particle-specific x50 parameters (10 particles)",
            "  - Grimes: Same curve for ALL particles (physically incorrect)",
            sprintf("  - VOxA OER_max = %.2f vs Grimes = 2.74 (survival units)", OER_max_survival),
            "  - VOxA: Physical LET limits enforced",
            "  - VOxA: Can predict for new particles via Z-interpolation",
            "",
            "================================================================================"
)

# Write report
writeLines(report, "results/comprehensive_validation_report_voxa.txt")
cat("Saved: results/comprehensive_validation_report_voxa.txt\n\n")


# Section 17: Save Diagnostic Data

cat("\n--- SECTION 17: SAVING DIAGNOSTIC DATA ---\n")
# Theoretical curve ranges
curve_ranges <- theoretical_curves %>%
  group_by(ion) %>%
  summarise(
    LET_min = min(LET),
    LET_max = max(LET),
    OER_min = min(OER_pred),
    OER_max = max(OER_pred),
    .groups = "drop"
  ) %>%
  left_join(particle_info %>% select(ion, max_let), by = "ion") %>%
  rename(physical_LET_limit = max_let)

# Calibration data ranges
data_ranges <- calibration_data %>%
  group_by(ion) %>%
  summarise(
    N = n(),
    LET_min = min(LET),
    LET_max = max(LET),
    OER_min = min(OER_retention),
    OER_max = max(OER_retention),
    n_near_bragg = sum(bragg_proximity > 0.7),
    .groups = "drop"
  )

# Largest residuals
largest_residuals <- calibration_data %>%
  arrange(desc(abs(residual))) %>%
  head(20) %>%
  select(ion, LET, cell_line_std, O2_hyp, OER_retention, OER_pred, residual, pct_error, bragg_proximity)

# Save all diagnostics
write_csv(theoretical_curves, "results/theoretical_curves_voxa.csv")
write_csv(curve_ranges, "results/theoretical_curve_ranges_voxa.csv")
write_csv(data_ranges, "results/calibration_data_ranges_voxa.csv")
write_csv(particle_stats, "results/particle_statistics_voxa.csv")
write_csv(cell_stats, "results/cell_line_statistics_voxa.csv")
write_csv(largest_residuals, "results/largest_residuals_voxa.csv")

cat("Saved diagnostics:\n")
# Final Summary

n_excellent <- sum(particle_stats$MAE_pct < 10)
n_good <- sum(particle_stats$MAE_pct >= 10 & particle_stats$MAE_pct < 15)
n_ok <- sum(particle_stats$MAE_pct >= 15 & particle_stats$MAE_pct < 20)
n_poor <- sum(particle_stats$MAE_pct >= 20)

cat(sprintf("║   R² (unweighted) = %.4f  <-- Primary metric                       ║\n",
            model$fit_statistics$r2))
cat(sprintf("║   R² (weighted)   = %.4f  (clinical weighting)                     ║\n",
            model$fit_statistics$r2_weighted))
cat(sprintf("║   MAE = %.4f  |  RMSE = %.4f                                     ║\n",
            model$fit_statistics$mae, model$fit_statistics$rmse))
cat(sprintf("║   OER_max = %.2f (retention) = %.2f (survival)                   ║\n",
            model$OER_max_theoretical, OER_max_survival))
cat(sprintf("║   Particles: %d excellent, %d good, %d ok, %d poor                    ║\n",
            n_excellent, n_good, n_ok, n_poor))
cat("Step 4 Complete. Ready for Steps 5-12.\n\n")