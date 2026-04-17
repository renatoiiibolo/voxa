# Step 5: Extended validation figures (13–20)
#
#   13 Case fraction evolution with LET
#   14 K-value sensitivity
#   15 P_indirect vs O2
#   16 OER heatmap (O2 × LET) for carbon
#   17 Particle-specific OER curves (overlaid)
#   18 Model component breakdown
#   19 x50 comparison with Z-ordering
#   20 x50 vs Z scaling
#
# Inputs:  results/uvaom_v8_corrected_model.rds,
#          results/calibration_data_v8_corrected.csv

library(tidyverse)
library(scales)
library(gridExtra)
library(grid)
library(viridis)

# Create output directories
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("results")) dir.create("results")


# Color Palette: "A Summer In Southern Italy"

# Case fraction colors (for p1, p2, p3)
case_fraction_colors <- c(
  "p1 (Direct)"            = "#C75B5B",
  "p2 (Hybrid)"            = "#1E5B8C",
  "p3 (Indirect)"          = "#2D7D46"
)

# Particle colors - consistent with Step 4
particle_colors <- c(
  "Photon"   = "#1E5B8C",
  "Proton"   = "#E07B3C",
  "Deuteron" = "#8B5A83",
  "Helium"   = "#D4A84B",
  "Carbon"   = "#2D7D46",
  "Nitrogen" = "#5DADE2",
  "Oxygen"   = "#48C9B0",
  "Neon"     = "#C75B5B",
  "Silicon"  = "#AF7AC5",
  "Argon"    = "#6B4E3D"
)

# x50 comparison colors
x50_colors <- c(
  "x50_dir (Direct)"   = "#C75B5B",
  "x50_ind (Indirect)" = "#D4A84B"
)

# Component breakdown colors (stacked bar)
component_colors <- c(
  "p1 (Direct)"              = "#C75B5B",
  "p2 × p_ind (Hybrid)"      = "#1E5B8C",
  "p3 × p_ind² (Indirect)"   = "#D4A84B"
)

# Accent colors
mediterranean_blue <- "#1E5B8C"
pompeii_red <- "#C75B5B"
terracotta <- "#E07B3C"
lemon_yellow <- "#D4A84B"
olive_green <- "#2D7D46"
ancient_stone <- "#4A3728"
warm_stone <- "#8B7355"
sandy_beige <- "#F5E6D3"
italian_cream <- "#FDF8F0"


# Physical Let Limits (Bragg Peak)

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


# Load Model

cat("\n--- LOADING MODEL ---\n")
model <- readRDS("results/uvaom_v8_corrected_model.rds")
params <- model$parameters
FIXED_PARAMS <- model$fixed_params

# Load conversion factor
if (file.exists("data/uvaom_recalibration_setup.RData")) {
  load("data/uvaom_recalibration_setup.RData")
} else {
  CONVERSION_FACTOR <- list(mean = 1.20, sd = 0.05)
}

cat(sprintf("Model: %s\n", model$version))
cat(sprintf("R² (unweighted) = %.4f\n", model$fit_statistics$r2))
cat(sprintf("R² (weighted) = %.4f\n", model$fit_statistics$r2_weighted))

# Calculate survival OER_max
OER_max_survival <- 1.0 + (model$OER_max_theoretical - 1.0) / CONVERSION_FACTOR$mean
cat(sprintf("OER_max (retention) = %.2f\n", model$OER_max_theoretical))
cat(sprintf("OER_max (survival) = %.2f\n\n", OER_max_survival))

# Extract key parameters
K_fix <- params$K_fix
K_repair <- params$K_repair
factor_HSG <- params$factor_HSG
factor_T1 <- params$factor_T1
factor_CHO <- params$factor_CHO
overkill_strength <- if (!is.null(params$overkill_strength)) params$overkill_strength else 0

cat(sprintf("K_fix = %.4f%%\n", K_fix))
cat(sprintf("K_repair = %.4f%%\n", K_repair))
cat(sprintf("Overkill strength = %.3f\n\n", overkill_strength))


# Model Functions

calc_p_indirect <- function(O2, K_fix, K_repair) {
  (O2 + K_fix) / (O2 + K_fix + K_repair)
}

calc_steepness_Z <- function(Z, s_base, s_scale) {
  s_base * (1 + s_scale * log(max(Z, 2) / 2))
}

compute_transitions <- function(LET, x50_dir, x50_ind, s_dir, s_ind) {
  x <- 2.5 * LET^1.1
  f_direct <- 1 / (1 + (x50_dir / max(x, 0.001))^s_dir)
  f_indirect <- 1 / (1 + (x50_ind / max(x, 0.001))^s_ind)
  return(list(f_direct = f_direct, f_indirect = f_indirect))
}

compute_case_fractions <- function(LET, x50_dir, x50_ind, s_dir, s_ind) {
  trans <- compute_transitions(LET, x50_dir, x50_ind, s_dir, s_ind)
  
  p1 <- FIXED_PARAMS$p1_low + (FIXED_PARAMS$p1_high - FIXED_PARAMS$p1_low) * trans$f_direct
  p3 <- FIXED_PARAMS$p3_low * (1 - trans$f_indirect)
  p2 <- 1 - p1 - p3
  p2 <- max(0, p2)
  
  total <- p1 + p2 + p3
  return(list(p1 = p1/total, p2 = p2/total, p3 = p3/total))
}

calc_P_DSB <- function(O2, LET, x50_dir, x50_ind, s_dir, s_ind, K_fix, K_repair) {
  fracs <- compute_case_fractions(LET, x50_dir, x50_ind, s_dir, s_ind)
  p_ind <- calc_p_indirect(O2, K_fix, K_repair)
  P_DSB <- fracs$p1 + fracs$p2 * p_ind + fracs$p3 * p_ind^2
  return(P_DSB)
}

calc_OER <- function(O2_hyp, O2_ref, LET, x50_dir, x50_ind, s_dir, s_ind, K_fix, K_repair) {
  P_hyp <- calc_P_DSB(O2_hyp, LET, x50_dir, x50_ind, s_dir, s_ind, K_fix, K_repair)
  P_ref <- calc_P_DSB(O2_ref, LET, x50_dir, x50_ind, s_dir, s_ind, K_fix, K_repair)
  return(max(1.0, P_ref / P_hyp))
}

#' Overkill correction factor
calc_overkill_factor <- function(LET, max_let, overkill_strength = 0.15) {
  proximity <- LET / max_let
  if (proximity > 0.7) {
    correction <- 1 + overkill_strength * ((proximity - 0.7) / 0.3)^2
    return(min(correction, 1 + overkill_strength * 1.5))
  } else {
    return(1.0)
  }
}

#' OER with overkill correction
calc_OER_with_overkill <- function(O2_hyp, O2_ref, LET, x50_dir, x50_ind, s_dir, s_ind, 
                                   K_fix, K_repair, max_let, overkill_strength) {
  OER_base <- calc_OER(O2_hyp, O2_ref, LET, x50_dir, x50_ind, s_dir, s_ind, K_fix, K_repair)
  overkill_factor <- calc_overkill_factor(LET, max_let, overkill_strength)
  OER_corrected <- 1 + (OER_base - 1) * overkill_factor
  return(max(1.0, OER_corrected))
}


# Figure 13: Case Fraction Evolution With Let

cat("\n--- FIGURE 13: Case Fraction Evolution with LET ---\n")
# Generate case fractions for each particle across LET range (within physical limits)
particles_for_fracs <- list(
  list(name = "Photon", ion = "photon", Z = 0, class = "light", max_let = 35),
  list(name = "Proton", ion = "proton", Z = 1, class = "light", max_let = 100),
  list(name = "Carbon", ion = "C", Z = 6, class = "heavy", max_let = 550),
  list(name = "Neon", ion = "Ne", Z = 10, class = "heavy", max_let = 700)
)

case_fraction_data <- data.frame()

for (p in particles_for_fracs) {
  # Use physical LET limits
  LET_min <- ifelse(p$class == "light", 0.1, 10)
  LET_seq <- 10^seq(log10(LET_min), log10(p$max_let), length.out = 100)
  
  x50_dir <- params[[paste0("x50_dir_", p$ion)]]
  x50_ind <- params[[paste0("x50_ind_", p$ion)]]
  
  if (p$class == "light") {
    s_dir <- params$s_dir_light
    s_ind <- params$s_ind_light
  } else {
    s_dir <- calc_steepness_Z(p$Z, params$s_dir_base, params$s_dir_scale)
    s_ind <- calc_steepness_Z(p$Z, params$s_ind_base, params$s_ind_scale)
  }
  
  for (LET in LET_seq) {
    fracs <- compute_case_fractions(LET, x50_dir, x50_ind, s_dir, s_ind)
    case_fraction_data <- rbind(case_fraction_data, data.frame(
      particle = p$name,
      LET = LET,
      max_let = p$max_let,
      p1 = fracs$p1,
      p2 = fracs$p2,
      p3 = fracs$p3
    ))
  }
}

# Pivot for plotting
case_fraction_long <- case_fraction_data %>%
  pivot_longer(cols = c(p1, p2, p3), names_to = "case", values_to = "fraction") %>%
  mutate(
    case_label = case_when(
      case == "p1" ~ "p1 (Direct)",
      case == "p2" ~ "p2 (Hybrid)",
      case == "p3" ~ "p3 (Indirect)"
    )
  )

fig13 <- ggplot(case_fraction_long, aes(x = LET, y = fraction, color = case_label)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ particle, scales = "free_x", ncol = 2) +
  scale_x_log10() +
  scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
  scale_color_manual(
    values = case_fraction_colors,
    name = "Case Fraction"
  ) +
  labs(
    title = "VOxA Model: Case Fraction Evolution with LET",
    subtitle = "p1 = d² (direct), p2 = 2di (hybrid), p3 = i² (indirect) | Curves cut at physical LET limits",
    x = expression(paste("LET (keV/", mu, "m)")),
    y = "Fraction"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = sandy_beige),
    strip.text = element_text(face = "bold", color = ancient_stone)
  ) +
  # Add horizontal reference lines for low-LET values
  geom_hline(yintercept = FIXED_PARAMS$p1_low, linetype = "dotted", color = pompeii_red, alpha = 0.5) +
  geom_hline(yintercept = FIXED_PARAMS$p2_low, linetype = "dotted", color = mediterranean_blue, alpha = 0.5) +
  geom_hline(yintercept = FIXED_PARAMS$p3_low, linetype = "dotted", color = olive_green, alpha = 0.5)

ggsave("figures/fig13_case_fraction_evolution_voxa.png", fig13, width = 12, height = 10, dpi = 300)
cat("Saved: figures/fig13_case_fraction_evolution_voxa.png\n\n")


# Figure 14: K-Value Sensitivity Analysis

cat("\n--- FIGURE 14: K-Value Sensitivity Analysis ---\n")
# Vary K_fix and K_repair and see effect on OER_max
K_fix_range <- seq(0.05, 0.50, length.out = 50)
K_repair_range <- seq(0.10, 0.60, length.out = 50)

sensitivity_grid <- expand.grid(K_fix = K_fix_range, K_repair = K_repair_range)

sensitivity_grid <- sensitivity_grid %>%
  filter(K_fix < K_repair) %>%  # Physical constraint
  rowwise() %>%
  mutate(
    p_ind_anoxia = calc_p_indirect(0.001, K_fix, K_repair),
    p_ind_normoxia = calc_p_indirect(21.0, K_fix, K_repair),
    P_anoxia = FIXED_PARAMS$p1_low + FIXED_PARAMS$p2_low * p_ind_anoxia + FIXED_PARAMS$p3_low * p_ind_anoxia^2,
    P_normoxia = FIXED_PARAMS$p1_low + FIXED_PARAMS$p2_low * p_ind_normoxia + FIXED_PARAMS$p3_low * p_ind_normoxia^2,
    OER_max = P_normoxia / P_anoxia
  ) %>%
  ungroup()

# Mark current model values
current_K_fix <- params$K_fix
current_K_repair <- params$K_repair

fig14 <- ggplot(sensitivity_grid, aes(x = K_fix, y = K_repair, fill = OER_max)) +
  geom_tile() +
  geom_contour(aes(z = OER_max), color = "white", alpha = 0.5, bins = 10) +
  # Mark current values
  geom_point(aes(x = current_K_fix, y = current_K_repair), 
             color = terracotta, size = 5, shape = 4, stroke = 2) +
  annotate("text", x = current_K_fix + 0.03, y = current_K_repair + 0.03,
           label = sprintf("Current\n(%.3f, %.3f)\nOER=%.2f", 
                           current_K_fix, current_K_repair, model$OER_max_theoretical),
           color = terracotta, size = 3, hjust = 0, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#0D3B66", "#1E5B8C", "#8B5A83", "#C75B5B", "#E07B3C", "#D4A84B"),
    name = expression(OER[max])
  ) +
  labs(
    title = "VOxA Model: K-Value Sensitivity Analysis",
    subtitle = expression(paste("Effect of ", K[fix], " and ", K[repair], " on theoretical ", OER[max], " (retention)")),
    x = expression(paste(K[fix], " (% ", O[2], ")")),
    y = expression(paste(K[repair], " (% ", O[2], ")"))
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right") +
  coord_fixed()

ggsave("figures/fig14_k_sensitivity_voxa.png", fig14, width = 10, height = 8, dpi = 300)
cat("Saved: figures/fig14_k_sensitivity_voxa.png\n\n")


# Figure 15: P_Ind Vs O2 Curve

cat("\n--- FIGURE 15: Oxygen Fixation Probability (p_ind vs O2) ---\n")
O2_range <- 10^seq(-3, log10(25), length.out = 200)

p_ind_curve <- tibble(
  O2 = O2_range,
  p_ind = sapply(O2_range, function(o2) calc_p_indirect(o2, K_fix, K_repair))
)

# Key points
key_points <- tibble(
  O2 = c(0.001, 0.021, 0.21, 2.1, 21),
  condition = c("Anoxia", "Severe", "Moderate", "Mild", "Normoxia")
) %>%
  mutate(p_ind = sapply(O2, function(o2) calc_p_indirect(o2, K_fix, K_repair)))

fig15 <- ggplot(p_ind_curve, aes(x = O2, y = p_ind)) +
  geom_line(color = mediterranean_blue, linewidth = 1.5) +
  geom_point(data = key_points, aes(x = O2, y = p_ind), 
             color = pompeii_red, size = 4) +
  geom_text(data = key_points, aes(x = O2, y = p_ind, label = condition),
            vjust = -1.5, size = 3.5, color = ancient_stone) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = warm_stone) +
  scale_x_log10(
    breaks = c(0.001, 0.01, 0.1, 1, 10),
    labels = c("0.001", "0.01", "0.1", "1", "10")
  ) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
  annotate("text", x = 0.003, y = 0.52, 
           label = "50% fixation", color = warm_stone, size = 3) +
  annotate("label", x = 0.1, y = 0.2,
           label = sprintf("K_fix = %.4f%%\nK_repair = %.4f%%\np_ind(anoxia) = %.3f\np_ind(normoxia) = %.3f",
                           K_fix, K_repair,
                           key_points$p_ind[1], key_points$p_ind[5]),
           size = 3.5, hjust = 0, fill = italian_cream, color = ancient_stone) +
  labs(
    title = "VOxA Model: Oxygen Fixation Probability",
    subtitle = expression(paste(p[ind], " = (", O[2], " + ", K[fix], ") / (", O[2], " + ", K[fix], " + ", K[repair], ")")),
    x = expression(paste("Oxygen Concentration (% ", O[2], ")")),
    y = expression(paste("Fixation Probability (", p[ind], ")"))
  ) +
  theme_bw(base_size = 12)

ggsave("figures/fig15_p_ind_curve_voxa.png", fig15, width = 10, height = 7, dpi = 300)
cat("Saved: figures/fig15_p_ind_curve_voxa.png\n\n")


# Figure 16: Oer Heatmap (O2 Vs Let) For Carbon

cat("\n--- FIGURE 16: OER Heatmap (O2 vs LET) for Carbon ---\n")
# Use carbon as representative particle (within physical limit)
x50_dir_C <- params$x50_dir_C
x50_ind_C <- params$x50_ind_C
s_dir_C <- calc_steepness_Z(6, params$s_dir_base, params$s_dir_scale)
s_ind_C <- calc_steepness_Z(6, params$s_ind_base, params$s_ind_scale)
max_let_C <- 550

O2_grid <- 10^seq(-2, log10(21), length.out = 50)
LET_grid <- 10^seq(log10(10), log10(max_let_C), length.out = 50)

heatmap_data <- expand.grid(O2 = O2_grid, LET = LET_grid) %>%
  rowwise() %>%
  mutate(
    OER = calc_OER(O2, 21.0, LET, x50_dir_C, x50_ind_C, s_dir_C, s_ind_C, K_fix, K_repair)
  ) %>%
  ungroup()

fig16 <- ggplot(heatmap_data, aes(x = LET, y = O2, fill = OER)) +
  geom_tile() +
  geom_contour(aes(z = OER), color = "white", alpha = 0.7, bins = 8) +
  scale_x_log10() +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10),
    labels = c("0.01", "0.1", "1", "10")
  ) +
  scale_fill_gradientn(
    colors = c("#0D3B66", "#1E5B8C", "#5B8C8C", "#D4A84B", "#E07B3C", "#C75B5B"),
    name = "OER\n(Retention)", 
    limits = c(1, 3.5)
  ) +
  # Add physical limit annotation
  geom_vline(xintercept = max_let_C, linetype = "dashed", color = "white", linewidth = 0.8) +
  annotate("text", x = max_let_C * 0.85, y = 0.015, label = "Bragg peak", 
           color = "white", size = 3, angle = 90, vjust = 0) +
  labs(
    title = "VOxA Model: OER Heatmap for Carbon Ions",
    subtitle = expression(paste("OER (retention) as function of ", O[2], " and LET | LET limited to 550 keV/μm")),
    x = expression(paste("LET (keV/", mu, "m)")),
    y = expression(paste("Oxygen (% ", O[2], ")"))
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave("figures/fig16_oer_heatmap_carbon_voxa.png", fig16, width = 10, height = 8, dpi = 300)
cat("Saved: figures/fig16_oer_heatmap_carbon_voxa.png\n\n")


# Figure 17: Particle-Specific Oer Curves (With Physical Let Limits)

cat("\n--- FIGURE 17: Particle-Specific OER Curves (with Physical LET Limits) ---\n")
particle_list <- list(
  list(name = "Photon", ion = "photon", Z = 0, class = "light", max_let = 35),
  list(name = "Proton", ion = "proton", Z = 1, class = "light", max_let = 100),
  list(name = "Helium", ion = "He", Z = 2, class = "heavy", max_let = 200),
  list(name = "Carbon", ion = "C", Z = 6, class = "heavy", max_let = 550),
  list(name = "Neon", ion = "Ne", Z = 10, class = "heavy", max_let = 700),
  list(name = "Argon", ion = "Ar", Z = 18, class = "heavy", max_let = 900)
)

oer_overlay_data <- data.frame()

for (p in particle_list) {
  # Generate LET range up to physical limit
  LET_min <- ifelse(p$class == "light", 0.1, 5)
  LET_seq <- 10^seq(log10(LET_min), log10(p$max_let), length.out = 150)
  
  x50_dir <- params[[paste0("x50_dir_", p$ion)]]
  x50_ind <- params[[paste0("x50_ind_", p$ion)]]
  
  if (is.null(x50_dir) || is.null(x50_ind)) next
  
  if (p$class == "light") {
    s_dir <- params$s_dir_light
    s_ind <- params$s_ind_light
  } else {
    s_dir <- calc_steepness_Z(p$Z, params$s_dir_base, params$s_dir_scale)
    s_ind <- calc_steepness_Z(p$Z, params$s_ind_base, params$s_ind_scale)
  }
  
  OER_vals <- sapply(LET_seq, function(l) {
    if (p$class == "heavy" && overkill_strength > 0) {
      calc_OER_with_overkill(0.001, 21.0, l, x50_dir, x50_ind, s_dir, s_ind, 
                             K_fix, K_repair, p$max_let, overkill_strength)
    } else {
      calc_OER(0.001, 21.0, l, x50_dir, x50_ind, s_dir, s_ind, K_fix, K_repair)
    }
  })
  
  oer_overlay_data <- rbind(oer_overlay_data, data.frame(
    particle = p$name,
    Z = p$Z,
    LET = LET_seq,
    OER = OER_vals,
    max_let = p$max_let
  ))
}

# Add endpoint markers
endpoint_data <- oer_overlay_data %>%
  group_by(particle) %>%
  filter(LET == max(LET)) %>%
  ungroup()

fig17 <- ggplot(oer_overlay_data, aes(x = LET, y = OER, color = particle)) +
  geom_line(linewidth = 1.2) +
  # Mark endpoints (Bragg peak)
  geom_point(data = endpoint_data, aes(x = LET, y = OER, color = particle), 
             size = 3, shape = 16) +
  scale_x_log10(
    breaks = c(0.1, 1, 10, 100, 1000),
    labels = c("0.1", "1", "10", "100", "1000")
  ) +
  scale_y_continuous(limits = c(1, 3.6), breaks = seq(1, 3.5, 0.5)) +
  scale_color_manual(values = particle_colors, name = "Particle") +
  geom_hline(yintercept = model$OER_max_theoretical, linetype = "dashed", color = warm_stone) +
  annotate("text", x = 0.15, y = model$OER_max_theoretical + 0.08, 
           label = sprintf("OER_max = %.2f (retention)", model$OER_max_theoretical),
           size = 3.5, color = ancient_stone) +
  labs(
    title = "VOxA Model: Particle-Specific OER Curves",
    subtitle = "OER (retention) at anoxia vs normoxia | Curves end at physical LET limits (Bragg peak)",
    x = expression(paste("LET (keV/", mu, "m)")),
    y = "OER (Retention)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave("figures/fig17_particle_oer_overlay_voxa.png", fig17, width = 12, height = 7, dpi = 300)
cat("Saved: figures/fig17_particle_oer_overlay_voxa.png\n\n")


# Figure 18: Model Component Breakdown

cat("\n--- FIGURE 18: Model Component Breakdown (P_DSB contributions) ---\n")
O2_levels <- c(0.001, 0.021, 0.21, 2.1, 21)
O2_labels <- c("Anoxia\n(0.001%)", "Severe\n(0.021%)", "Moderate\n(0.21%)", 
               "Mild\n(2.1%)", "Normoxia\n(21%)")

component_data <- data.frame()

for (i in seq_along(O2_levels)) {
  O2 <- O2_levels[i]
  p_ind <- calc_p_indirect(O2, K_fix, K_repair)
  
  p1 <- FIXED_PARAMS$p1_low
  p2 <- FIXED_PARAMS$p2_low
  p3 <- FIXED_PARAMS$p3_low
  
  contrib_p1 <- p1
  contrib_p2 <- p2 * p_ind
  contrib_p3 <- p3 * p_ind^2
  
  P_DSB <- contrib_p1 + contrib_p2 + contrib_p3
  
  component_data <- rbind(component_data, data.frame(
    O2 = O2,
    O2_label = O2_labels[i],
    component = c("p1 (Direct)", "p2 × p_ind (Hybrid)", "p3 × p_ind² (Indirect)"),
    contribution = c(contrib_p1, contrib_p2, contrib_p3),
    P_DSB = P_DSB
  ))
}

component_data$O2_label <- factor(component_data$O2_label, levels = O2_labels)
component_data$component <- factor(component_data$component, 
                                   levels = c("p3 × p_ind² (Indirect)", 
                                              "p2 × p_ind (Hybrid)", 
                                              "p1 (Direct)"))

fig18 <- ggplot(component_data, aes(x = O2_label, y = contribution, fill = component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(
    values = component_colors,
    name = "Component"
  ) +
  scale_y_continuous(limits = c(0, 1.05), labels = percent_format()) +
  geom_text(data = component_data %>% group_by(O2_label) %>% summarise(P_DSB = first(P_DSB)),
            aes(x = O2_label, y = P_DSB + 0.03, label = sprintf("%.3f", P_DSB)),
            inherit.aes = FALSE, size = 3.5, fontface = "bold", color = ancient_stone) +
  labs(
    title = "VOxA Model: P_DSB Component Breakdown at Low LET",
    subtitle = "P_DSB = p1 + p2·p_ind + p3·p_ind² | Values shown at top of bars",
    x = "Oxygen Condition",
    y = expression(paste("Contribution to ", P[DSB]))
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 10)
  )

ggsave("figures/fig18_component_breakdown_voxa.png", fig18, width = 10, height = 7, dpi = 300)
cat("Saved: figures/fig18_component_breakdown_voxa.png\n\n")


# Figure 19: X50 Values Comparison Across Particles

cat("\n--- FIGURE 19: x50 Values Comparison Across Particles ---\n")
x50_data <- tibble(
  particle = c("Photon", "Proton", "Deuteron", "Helium", "Carbon", "Nitrogen", "Oxygen", "Neon", "Silicon", "Argon"),
  Z = c(0, 1, 1, 2, 6, 7, 8, 10, 14, 18),
  x50_dir = c(
    params$x50_dir_photon, params$x50_dir_proton, params$x50_dir_deuteron,
    params$x50_dir_He, params$x50_dir_C, 
    model$interpolated_params$N$x50_dir,
    model$interpolated_params$O$x50_dir,
    params$x50_dir_Ne, 
    model$interpolated_params$Si$x50_dir,
    params$x50_dir_Ar
  ),
  x50_ind = c(
    params$x50_ind_photon, params$x50_ind_proton, params$x50_ind_deuteron,
    params$x50_ind_He, params$x50_ind_C,
    model$interpolated_params$N$x50_ind,
    model$interpolated_params$O$x50_ind,
    params$x50_ind_Ne,
    model$interpolated_params$Si$x50_ind,
    params$x50_ind_Ar
  ),
  calib_type = c("calibrated", "calibrated", "calibrated", "calibrated", "calibrated",
                 "interpolated", "interpolated", "calibrated", "interpolated", "calibrated")
) %>%
  mutate(ratio = x50_ind / x50_dir)

x50_long <- x50_data %>%
  pivot_longer(cols = c(x50_dir, x50_ind), names_to = "type", values_to = "x50") %>%
  mutate(
    type_label = ifelse(type == "x50_dir", "x50_dir (Direct)", "x50_ind (Indirect)")
  )

x50_data$particle <- factor(x50_data$particle, 
                            levels = c("Photon", "Proton", "Deuteron", "Helium", "Carbon", 
                                       "Nitrogen", "Oxygen", "Neon", "Silicon", "Argon"))
x50_long$particle <- factor(x50_long$particle, 
                            levels = c("Photon", "Proton", "Deuteron", "Helium", "Carbon",
                                       "Nitrogen", "Oxygen", "Neon", "Silicon", "Argon"))

fig19 <- ggplot(x50_long, aes(x = particle, y = x50, fill = type_label)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = x50_colors,
    name = "Transition"
  ) +
  scale_y_log10() +
  geom_text(data = x50_data, 
            aes(x = particle, y = x50_ind * 1.3, label = sprintf("%.1f×", ratio)),
            inherit.aes = FALSE, size = 2.5, vjust = 0, color = ancient_stone) +
  # Mark interpolated particles
  geom_point(data = x50_data %>% filter(calib_type == "interpolated"),
             aes(x = particle, y = 50), 
             inherit.aes = FALSE, shape = 8, size = 2, color = warm_stone) +
  labs(
    title = "VOxA Model: x50 Transition Parameters by Particle",
    subtitle = "x50_dir < x50_ind ensures direct transition occurs before indirect | ★ = Z-interpolated | Numbers = ratio",
    x = "Particle",
    y = expression(paste(x[50], " (radiation quality units)"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("figures/fig19_x50_comparison_voxa.png", fig19, width = 12, height = 7, dpi = 300)
cat("Saved: figures/fig19_x50_comparison_voxa.png\n\n")


# Figure 20: X50_Dir Vs Z (Z-Ordering Validation)

cat("\n--- FIGURE 20: x50_dir vs Z (Z-Ordering Validation) ---\n")
# Focus on heavy ions to show Z-ordering
heavy_x50 <- x50_data %>%
  filter(Z >= 2) %>%
  arrange(Z)

fig20 <- ggplot(heavy_x50, aes(x = Z, y = x50_dir)) +
  geom_line(color = pompeii_red, linewidth = 1.2) +
  geom_point(aes(shape = calib_type), color = pompeii_red, size = 4) +
  scale_shape_manual(values = c("calibrated" = 16, "interpolated" = 8), 
                     name = "Type") +
  geom_text(aes(label = particle), vjust = -1, size = 3.5, color = ancient_stone) +
  scale_x_continuous(breaks = c(2, 6, 7, 8, 10, 14, 18)) +
  labs(
    title = "VOxA Model: Z-Ordering of x50_dir (Physics Validation)",
    subtitle = "x50_dir INCREASES with Z | Physics: At same LET, lighter ions have lower OER (denser tracks)",
    x = "Atomic Number (Z)",
    y = expression(paste(x[50]^{dir}, " (radiation quality units)"))
  ) +
  annotate("text", x = 5, y = max(heavy_x50$x50_dir) * 0.9,
           label = "x50_dir(He) < x50_dir(C) < x50_dir(Ne) < x50_dir(Ar)\n\nLighter ions transition to direct damage\nat LOWER LET (smaller x50)",
           hjust = 0, size = 3.5, color = ancient_stone) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("figures/fig20_x50_z_ordering_voxa.png", fig20, width = 10, height = 7, dpi = 300)
cat("Saved: figures/fig20_x50_z_ordering_voxa.png\n\n")


# Save Extended Data

cat("\n--- SAVING EXTENDED DATA ---\n")
write_csv(case_fraction_data, "results/case_fraction_evolution_voxa.csv")
write_csv(sensitivity_grid, "results/k_sensitivity_grid_voxa.csv")
write_csv(p_ind_curve, "results/p_ind_curve_voxa.csv")
write_csv(heatmap_data, "results/oer_heatmap_carbon_voxa.csv")
write_csv(oer_overlay_data, "results/particle_oer_curves_voxa.csv")
write_csv(x50_data, "results/x50_parameters_voxa.csv")

cat("Saved:\n")
# Extended Validation Report

cat("\n--- GENERATING EXTENDED VALIDATION REPORT ---\n")
# Calculate additional statistics for extended report
p_ind_anoxia <- calc_p_indirect(0.001, K_fix, K_repair)
p_ind_normoxia <- calc_p_indirect(21.0, K_fix, K_repair)

# Case fractions at low and high LET for carbon
fracs_low_let <- compute_case_fractions(20, params$x50_dir_C, params$x50_ind_C, 
                                        calc_steepness_Z(6, params$s_dir_base, params$s_dir_scale),
                                        calc_steepness_Z(6, params$s_ind_base, params$s_ind_scale))
fracs_high_let <- compute_case_fractions(500, params$x50_dir_C, params$x50_ind_C,
                                         calc_steepness_Z(6, params$s_dir_base, params$s_dir_scale),
                                         calc_steepness_Z(6, params$s_ind_base, params$s_ind_scale))

# Steepness values for each heavy ion
steepness_data <- tibble(
  ion = c("He", "C", "Ne", "Ar"),
  Z = c(2, 6, 10, 18),
  s_dir = sapply(c(2, 6, 10, 18), function(z) calc_steepness_Z(z, params$s_dir_base, params$s_dir_scale)),
  s_ind = sapply(c(2, 6, 10, 18), function(z) calc_steepness_Z(z, params$s_ind_base, params$s_ind_scale))
)

# Create extended report
extended_report <- c(
  "================================================================================",
  "         VOxA MODEL EXTENDED VALIDATION REPORT (Step 5 Supplement)",
  "         Variable Oxygen-dependent Amorphous Track Model",
  "================================================================================",
  "",
  sprintf("Generated: %s", Sys.time()),
  sprintf("Model version: %s", model$version),
  sprintf("This report extends: comprehensive_validation_report_voxa.txt"),
  "",
  "================================================================================",
  "PART A: OXYGEN FIXATION KINETICS",
  "================================================================================",
  "",
  "1. MICHAELIS-MENTEN OXYGEN KINETICS",
  "--------------------------------------------------------------------------------",
  "",
  "The probability of oxygen fixation (indirect damage retention) follows:",
  "",
  "  p_ind = (O2 + K_fix) / (O2 + K_fix + K_repair)",
  "",
  "Calibrated parameters:",
  sprintf("  K_fix    = %.4f%% O2", K_fix),
  sprintf("  K_repair = %.4f%% O2", K_repair),
  sprintf("  Ratio K_repair/K_fix = %.2f", K_repair / K_fix),
  "",
  "Oxygen fixation probability at key O2 levels:",
  "┌────────────────┬────────────┬────────────┐",
  "│ Condition      │ O2 (%)     │ p_ind      │",
  "├────────────────┼────────────┼────────────┤",
  sprintf("│ Anoxia         │ %10.4f │ %10.4f │", 0.001, p_ind_anoxia),
  sprintf("│ Severe hypoxia │ %10.4f │ %10.4f │", 0.021, calc_p_indirect(0.021, K_fix, K_repair)),
  sprintf("│ Moderate       │ %10.4f │ %10.4f │", 0.21, calc_p_indirect(0.21, K_fix, K_repair)),
  sprintf("│ Mild hypoxia   │ %10.4f │ %10.4f │", 2.1, calc_p_indirect(2.1, K_fix, K_repair)),
  sprintf("│ Normoxia       │ %10.4f │ %10.4f │", 21.0, p_ind_normoxia),
  "└────────────────┴────────────┴────────────┘",
  "",
  "Physical interpretation:",
  sprintf("  At anoxia:   %.1f%% of indirect damage is fixed (retained)", p_ind_anoxia * 100),
  sprintf("  At normoxia: %.1f%% of indirect damage is fixed (retained)", p_ind_normoxia * 100),
  sprintf("  Difference:  %.1f percentage points", (p_ind_normoxia - p_ind_anoxia) * 100),
  "",
  "================================================================================",
  "PART B: CASE FRACTION ANALYSIS",
  "================================================================================",
  "",
  "2. DSB COMBINATORICS (Sakata et al. 2019)",
  "--------------------------------------------------------------------------------",
  "",
  "DSBs form when two SSBs occur within ~10 bp. The probability of each SSB",
  "being direct (d) or indirect (i) determines the case fractions:",
  "",
  "  d ≈ 0.20 (single-hit direct probability)",
  "  i ≈ 0.80 (single-hit indirect probability)",
  "",
  "Low-LET case fractions (baseline):",
  sprintf("  p1_low = d²   = %.2f (both hits direct → O2-independent)", FIXED_PARAMS$p1_low),
  sprintf("  p2_low = 2di  = %.2f (one hit each → linear O2 dependence)", FIXED_PARAMS$p2_low),
  sprintf("  p3_low = i²   = %.2f (both hits indirect → quadratic O2 dependence)", FIXED_PARAMS$p3_low),
  "",
  "High-LET limit (Nikjoo 2001):",
  sprintf("  p1_high = %.2f (direct fraction at extreme LET)", FIXED_PARAMS$p1_high),
  "  Note: ~32% indirect action remains even at 2106 keV/μm",
  "",
  "3. CASE FRACTION EVOLUTION WITH LET (Carbon example)",
  "--------------------------------------------------------------------------------",
  "",
  "At LET = 20 keV/μm (entrance region):",
  sprintf("  p1 = %.4f (%.1f%%)", fracs_low_let$p1, fracs_low_let$p1 * 100),
  sprintf("  p2 = %.4f (%.1f%%)", fracs_low_let$p2, fracs_low_let$p2 * 100),
  sprintf("  p3 = %.4f (%.1f%%)", fracs_low_let$p3, fracs_low_let$p3 * 100),
  "",
  "At LET = 500 keV/μm (near Bragg peak):",
  sprintf("  p1 = %.4f (%.1f%%)", fracs_high_let$p1, fracs_high_let$p1 * 100),
  sprintf("  p2 = %.4f (%.1f%%)", fracs_high_let$p2, fracs_high_let$p2 * 100),
  sprintf("  p3 = %.4f (%.1f%%)", fracs_high_let$p3, fracs_high_let$p3 * 100),
  "",
  "Change from low to high LET:",
  sprintf("  Δp1 = %+.4f (direct damage INCREASES)", fracs_high_let$p1 - fracs_low_let$p1),
  sprintf("  Δp2 = %+.4f (hybrid damage)", fracs_high_let$p2 - fracs_low_let$p2),
  sprintf("  Δp3 = %+.4f (indirect damage DECREASES)", fracs_high_let$p3 - fracs_low_let$p3),
  "",
  "================================================================================",
  "PART C: PARTICLE-SPECIFIC PARAMETERS",
  "================================================================================",
  "",
  "4. x50 TRANSITION PARAMETERS",
  "--------------------------------------------------------------------------------",
  "",
  "The x50 parameters control when the direct/indirect transitions occur.",
  "x50_dir < x50_ind ensures direct damage transition happens first.",
  "",
  "Calibrated x50 values:",
  "┌────────────┬─────┬───────────┬───────────┬─────────┬──────────────┐",
  "│ Particle   │  Z  │  x50_dir  │  x50_ind  │  Ratio  │ Type         │",
  "├────────────┼─────┼───────────┼───────────┼─────────┼──────────────┤"
)

# Add x50 data rows
for (i in 1:nrow(x50_data)) {
  row <- x50_data[i, ]
  extended_report <- c(extended_report,
                       sprintf("│ %-10s │ %3d │ %9.1f │ %9.1f │ %7.1f │ %-12s │",
                               row$particle, row$Z, row$x50_dir, row$x50_ind, row$ratio, row$calib_type))
}

extended_report <- c(extended_report,
                     "└────────────┴─────┴───────────┴───────────┴─────────┴──────────────┘",
                     "",
                     "5. Z-ORDERING PHYSICS CONSTRAINT",
                     "--------------------------------------------------------------------------------",
                     "",
                     "CONSTRAINT: x50_dir must INCREASE with Z for heavy ions",
                     "",
                     "Physical rationale:",
                     "  At the same LET, lighter ions are SLOWER (since LET ∝ z²/v²)",
                     "  → Slower ions have narrower track cores (R_core ∝ β)",
                     "  → Higher local ionization density in smaller volume",
                     "  → More radical-radical recombination (•OH + •OH → H₂O₂)",
                     "  → Less indirect damage available for oxygen fixation",
                     "  → LOWER OER at the same LET",
                     "",
                     "  Therefore, lighter ions transition to direct damage at LOWER LET,",
                     "  which means SMALLER x50_dir values.",
                     "",
                     "Verification (heavy ions only):",
                     sprintf("  x50_dir(He)  = %6.1f (Z=2)", params$x50_dir_He),
                     sprintf("  x50_dir(C)   = %6.1f (Z=6)   [%.1f× He]", params$x50_dir_C, params$x50_dir_C / params$x50_dir_He),
                     sprintf("  x50_dir(Ne)  = %6.1f (Z=10)  [%.1f× He]", params$x50_dir_Ne, params$x50_dir_Ne / params$x50_dir_He),
                     sprintf("  x50_dir(Ar)  = %6.1f (Z=18)  [%.1f× He]", params$x50_dir_Ar, params$x50_dir_Ar / params$x50_dir_He),
                     "",
                     sprintf("  Z-ORDERING STATUS: %s", ifelse(model$z_ordering_satisfied, "✓ SATISFIED", "✗ VIOLATED")),
                     "",
                     "6. STEEPNESS PARAMETERS (HYBRID MODEL)",
                     "--------------------------------------------------------------------------------",
                     "",
                     "Light particles (photon, proton, deuteron) use fixed steepness:",
                     sprintf("  s_dir_light = %.3f", params$s_dir_light),
                     sprintf("  s_ind_light = %.3f", params$s_ind_light),
                     "",
                     "Heavy particles use Z-scaled steepness:",
                     sprintf("  s_dir(Z) = %.3f × (1 + %.3f × log(Z/2))", params$s_dir_base, params$s_dir_scale),
                     sprintf("  s_ind(Z) = %.3f × (1 + %.3f × log(Z/2))", params$s_ind_base, params$s_ind_scale),
                     "",
                     "Effective steepness for heavy ions:",
                     "┌──────────┬─────┬─────────┬─────────┐",
                     "│ Ion      │  Z  │  s_dir  │  s_ind  │",
                     "├──────────┼─────┼─────────┼─────────┤"
)

for (i in 1:nrow(steepness_data)) {
  row <- steepness_data[i, ]
  extended_report <- c(extended_report,
                       sprintf("│ %-8s │ %3d │ %7.3f │ %7.3f │", row$ion, row$Z, row$s_dir, row$s_ind))
}

extended_report <- c(extended_report,
                     "└──────────┴─────┴─────────┴─────────┘",
                     "",
                     "================================================================================",
                     "PART D: PHYSICAL LET LIMITS",
                     "================================================================================",
                     "",
                     "7. BRAGG PEAK CUTOFFS",
                     "--------------------------------------------------------------------------------",
                     "",
                     "Each particle has a maximum achievable LET determined by its Bragg peak.",
                     "The VOxA model enforces these physical limits:",
                     "",
                     "┌────────────┬────────────────┬──────────────────────────────────────┐",
                     "│ Particle   │ Max LET        │ Notes                                │",
                     "│            │ (keV/μm)       │                                      │",
                     "├────────────┼────────────────┼──────────────────────────────────────┤",
                     "│ Photon     │             35 │ Secondary electrons (δ-rays)         │",
                     "│ Proton     │            100 │ Clinical SOBP endpoint               │",
                     "│ Deuteron   │            120 │ Slightly higher than proton          │",
                     "│ Helium     │            200 │ He-4 ion                             │",
                     "│ Carbon     │            550 │ Clinical carbon therapy              │",
                     "│ Nitrogen   │            600 │ Interpolated                         │",
                     "│ Oxygen     │            620 │ Interpolated                         │",
                     "│ Neon       │            700 │ Research beams                       │",
                     "│ Silicon    │            800 │ Interpolated                         │",
                     "│ Argon      │            900 │ Research beams                       │",
                     "└────────────┴────────────────┴──────────────────────────────────────┘",
                     "",
                     "Implications:",
                     "  • Model predictions are only valid up to these LET values",
                     "  • Theoretical curves are cut at Bragg peak limits",
                     "  • Extrapolation beyond limits is physically meaningless",
                     "",
                     "8. OVERKILL CORRECTION",
                     "--------------------------------------------------------------------------------",
                     "",
                     sprintf("Overkill strength parameter: %.3f", overkill_strength),
                     "",
                     ifelse(overkill_strength > 0,
                            paste0("The overkill correction accounts for reduced biological effectiveness\n",
                                   "near the Bragg peak due to 'wasted' dose in already-killed cells.\n",
                                   "Applied when LET > 70% of physical limit for heavy ions."),
                            paste0("Overkill correction is disabled (strength = 0).\n",
                                   "The calibration data did not require this correction.")),
                     "",
                     "================================================================================",
                     "PART E: K-VALUE SENSITIVITY ANALYSIS",
                     "================================================================================",
                     "",
                     "9. EFFECT OF K VALUES ON OER_max",
                     "--------------------------------------------------------------------------------",
                     "",
                     "The theoretical OER_max depends on K_fix and K_repair through p_ind:",
                     "",
                     "  OER_max = P_DSB(normoxia) / P_DSB(anoxia)",
                     "",
                     "where P_DSB = p1 + p2·p_ind + p3·p_ind²",
                     "",
                     "Current model:",
                     sprintf("  K_fix = %.4f%%, K_repair = %.4f%%", K_fix, K_repair),
                     sprintf("  OER_max (retention) = %.2f", model$OER_max_theoretical),
                     sprintf("  OER_max (survival)  = %.2f", OER_max_survival),
                     "",
                     "Sensitivity (at low LET):",
                     "  • Increasing K_fix → DECREASES OER_max (more fixation at anoxia)",
                     "  • Increasing K_repair → INCREASES OER_max (less fixation at anoxia)",
                     "  • Constraint: K_fix < K_repair (physical requirement)",
                     "",
                     "Comparison to literature:",
                     sprintf("  VOxA OER_max (survival)  = %.2f", OER_max_survival),
                     "  Grimes (2020)            = 2.74",
                     "  Wenzl & Wilkens (2011)   = 2.5-3.0 (depending on model)",
                     "  Experimental consensus   = 2.5-3.0 (photons, anoxia)",
                     "",
                     "================================================================================",
                     "PART F: P_DSB COMPONENT BREAKDOWN",
                     "================================================================================",
                     "",
                     "10. CONTRIBUTIONS TO DNA DAMAGE RETENTION",
                     "--------------------------------------------------------------------------------",
                     "",
                     "P_DSB represents the probability of DNA damage being retained (not repaired).",
                     "It has three components with different O2 dependencies:",
                     "",
                     "  P_DSB = p1 + p2·p_ind + p3·p_ind²",
                     "",
                     "At low LET (baseline case fractions):",
                     "",
                     "┌────────────────┬────────────┬────────────┬────────────┬────────────┐",
                     "│ Condition      │ p1         │ p2·p_ind   │ p3·p_ind²  │ P_DSB      │",
                     "│                │ (O2-indep) │ (linear)   │ (quadratic)│ (total)    │",
                     "├────────────────┼────────────┼────────────┼────────────┼────────────┤"
)

# Add component breakdown rows
for (cond in c("Anoxia", "Severe", "Moderate", "Mild", "Normoxia")) {
  O2_val <- switch(cond,
                   "Anoxia" = 0.001,
                   "Severe" = 0.021,
                   "Moderate" = 0.21,
                   "Mild" = 2.1,
                   "Normoxia" = 21.0
  )
  p_ind_val <- calc_p_indirect(O2_val, K_fix, K_repair)
  p1_contrib <- FIXED_PARAMS$p1_low
  p2_contrib <- FIXED_PARAMS$p2_low * p_ind_val
  p3_contrib <- FIXED_PARAMS$p3_low * p_ind_val^2
  P_DSB_val <- p1_contrib + p2_contrib + p3_contrib
  
  extended_report <- c(extended_report,
                       sprintf("│ %-14s │ %10.4f │ %10.4f │ %10.4f │ %10.4f │",
                               cond, p1_contrib, p2_contrib, p3_contrib, P_DSB_val))
}

extended_report <- c(extended_report,
                     "└────────────────┴────────────┴────────────┴────────────┴────────────┘",
                     "",
                     "Key observations:",
                     "  • p1 (direct) is constant across all O2 levels",
                     "  • p2·p_ind (hybrid) increases linearly with O2",
                     "  �� p3·p_ind² (indirect) increases quadratically with O2",
                     "  • At anoxia, indirect damage is minimally fixed (low p_ind)",
                     "  • At normoxia, indirect damage is maximally fixed (high p_ind)",
                     "",
                     "================================================================================",
                     "PART G: FIGURES GENERATED IN STEP 5",
                     "================================================================================",
                     "",
                     "11. FIGURE LIST",
                     "--------------------------------------------------------------------------------",
                     "",
                     "Figure 13: Case Fraction Evolution with LET",
                     "  - Shows how p1, p2, p3 change with LET for different particles",
                     "  - Demonstrates transition from indirect-dominated to direct-dominated",
                     "  - Curves cut at physical LET limits",
                     "",
                     "Figure 14: K-Value Sensitivity Analysis",
                     "  - Heatmap showing OER_max as function of K_fix and K_repair",
                     "  - Current model values marked",
                     "  - Contour lines for equal OER_max",
                     "",
                     "Figure 15: Oxygen Fixation Probability (p_ind vs O2)",
                     "  - Michaelis-Menten curve shape",
                     "  - Key oxygen conditions marked",
                     "  - Shows steep transition in hypoxic region",
                     "",
                     "Figure 16: OER Heatmap (O2 vs LET) for Carbon",
                     "  - 2D visualization of OER across O2-LET space",
                     "  - Demonstrates both O2 and LET dependencies",
                     "  - Physical LET limit annotated",
                     "",
                     "Figure 17: Particle-Specific OER Curves",
                     "  - All particles overlaid on same axes",
                     "  - Each curve ends at its physical LET limit",
                     "  - Demonstrates particle-specific track structure effects",
                     "",
                     "Figure 18: Model Component Breakdown",
                     "  - Stacked bar chart of P_DSB contributions",
                     "  - Shows relative importance of each component vs O2",
                     "  - Demonstrates why OER changes with O2",
                     "",
                     "Figure 19: x50 Values Comparison",
                     "  - Bar chart comparing x50_dir and x50_ind across particles",
                     "  - Interpolated particles marked",
                     "  - Ratio (x50_ind/x50_dir) shown",
                     "",
                     "Figure 20: x50_dir vs Z (Z-Ordering Validation)",
                     "  - Line plot showing x50_dir increases with Z",
                     "  - Validates physics constraint",
                     "  - Key for understanding particle-specific OER",
                     "",
                     "================================================================================",
                     "PART H: SUMMARY",
                     "================================================================================",
                     "",
                     sprintf("Model: VOxA (Variable Oxygen-dependent Amorphous Track Model)"),
                     sprintf("R² (unweighted): %.4f", model$fit_statistics$r2),
                     sprintf("R² (weighted):   %.4f", model$fit_statistics$r2_weighted),
                     "",
                     "Key physics captured:",
                     "  ✓ DSB combinatorics (Sakata 2019)",
                     "  ✓ Michaelis-Menten oxygen kinetics",
                     "  ✓ LET-dependent case fraction evolution",
                     "  ✓ Particle-specific track structure (x50 parameters)",
                     "  ✓ Z-ordering constraint (lighter ions = lower OER at same LET)",
                     "  ✓ Physical LET limits (Bragg peak cutoffs)",
                     "  ✓ Hybrid light/heavy particle model",
                     "",
                     "Advantages over previous models (e.g., Grimes 2020):",
                     "  • Particle-specific parameters (not one curve for all)",
                     "  • Physical LET limits enforced",
                     "  • Z-interpolation for unmeasured particles",
                     "  • Intermediate O2 levels validated (Tinganelli data)",
                     "  • Explicit distinction between retention and survival OER",
                     "",
                     "================================================================================"
)

# Write extended report
writeLines(extended_report, "results/comprehensive_validation_report_voxa_extended.txt")
cat("Saved: results/comprehensive_validation_report_voxa_extended.txt\n\n")


# Final Summary

cat("╚══════════════════════════════════════════════════════════════════════╝\n")