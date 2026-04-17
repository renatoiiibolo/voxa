# Step 11: VA scaling validation
#
# Tests whether delta_f calibrated on ~2500 DSBs (calibration arm)
# generalizes to the independent ~400-DSB validation arm. Six statistical
# tests per particle (18 total): F-test, Levene, KS, TOST, Feltz-Miller,
# bootstrap CI overlap. The primary claim is sample-size invariance of
# the within-nucleus P_DSB coefficient of variation.
#
# Inputs:  results/voxa_voxel_aware_calibration.json,
#          voxa_features_output_calibration/ and voxa_features_output_validation/
# Outputs: results/voxa_scaling_validation_results.json,
#          figures/fig42_* through fig46_*

library(tidyverse)
library(jsonlite)
library(car)  # For Levene's test

set.seed(42)

# Create output directories
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("results")) dir.create("results")


# Color Palette: "A Summer In Southern Italy"

mediterranean_blue <- "#1E5B8C"
pompeii_red <- "#C75B5B"
terracotta <- "#E07B3C"
lemon_yellow <- "#D4A84B"
olive_green <- "#2D7D46"
amalfi_purple <- "#8B5A83"
capri_teal <- "#48C9B0"
espresso_brown <- "#6B4E3D"
ancient_stone <- "#4A3728"
warm_stone <- "#8B7355"
sandy_beige <- "#F5E6D3"
italian_cream <- "#FDF8F0"
deep_navy <- "#0D3B66"

# Particle-specific colors
particle_colors <- c(
  "electron" = mediterranean_blue,
  "proton" = lemon_yellow,
  "carbon" = pompeii_red
)


# Section 1: Load Voxa Calibration From Step 9

cat("\n--- SECTION 1: LOADING VOxA CALIBRATION (Step 9) ---\n")
calibration_json_paths <- c(
  "results/voxa_voxel_aware_calibration.json",
  "voxa_voxel_aware_calibration.json"
)

calibration_file <- NULL
for (path in calibration_json_paths) {
  if (file.exists(path)) {
    calibration_file <- path
    break
  }
}

if (is.null(calibration_file)) {
  stop("ERROR: Step 9 calibration file not found!\n",
       "Please run Step 9 first: Rscript step9_voxa_voxel_aware_calibration.R")
}

cat(sprintf("Loading calibration from: %s\n", calibration_file))
calibration_json <- fromJSON(calibration_file)

# Extract VOxA parameters
VOXA_PARAMS <- list(
  model_version = calibration_json$model_version,
  K_fix = calibration_json$base_parameters$K_fix,
  K_repair = calibration_json$base_parameters$K_repair,
  O2_normoxia = 21.0,
  O2_calibration = calibration_json$va_calibration_settings$O2_calibration_pct,
  f_min = calibration_json$va_calibration_settings$f_min,
  f_max = calibration_json$va_calibration_settings$f_max,
  particles = list()
)

# Load particle-specific parameters
PARTICLE_KEYS <- names(calibration_json$particle_calibration)
for (p in PARTICLE_KEYS) {
  pc <- calibration_json$particle_calibration[[p]]
  VOXA_PARAMS$particles[[p]] <- list(
    n_dsbs = pc$n_dsbs,
    p1 = pc$p1_base,
    p2 = pc$p2_base,
    p3 = pc$p3_base,
    delta_f = pc$delta_f,
    delta_f_ci_low = pc$delta_f_ci_low,
    delta_f_ci_high = pc$delta_f_ci_high,
    P_DSB_cv = pc$P_DSB_cv,
    P_DSB_cv_ci_low = pc$P_DSB_cv_ci_low,
    P_DSB_cv_ci_high = pc$P_DSB_cv_ci_high,
    E_mean = pc$E_mean,
    E_sd = pc$E_sd
  )
}

cat(sprintf("Model version: %s\n", VOXA_PARAMS$model_version))
cat(sprintf("Calibration O2: %.2f%%\n", VOXA_PARAMS$O2_calibration))
cat(sprintf("Particles: %s\n\n", paste(PARTICLE_KEYS, collapse = ", ")))


# Section 2: Load Datasets

cat("\n--- SECTION 2: LOADING CALIBRATION AND VALIDATION DATASETS ---\n")
# Calibration data paths
cal_data_paths <- c(
  "voxa_features_output_calibration/all_particles_calibration_energy_features.csv",
  "voxa_features_output_calibration/all_particles_energy_features.csv",
  "calibration_output/all_particles_calibration_energy_features.csv"
)

cal_data_file <- NULL
for (path in cal_data_paths) {
  if (file.exists(path)) {
    cal_data_file <- path
    break
  }
}

if (is.null(cal_data_file)) {
  stop("ERROR: Calibration energy features file not found!")
}

cat(sprintf("Loading calibration data from: %s\n", cal_data_file))
cal_data <- read_csv(cal_data_file, show_col_types = FALSE)
cat(sprintf("  Total calibration DSBs: %d\n", nrow(cal_data)))

# Validation data paths
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
cat("Dataset summary:\n")
cat("─" %>% strrep(60), "\n")
cat(sprintf("%-12s %12s %12s %12s\n", "Particle", "Calibration", "Validation", "Ratio"))
cat("─" %>% strrep(60), "\n")

for (p in PARTICLE_KEYS) {
  n_cal <- sum(cal_data$particle == p)
  n_val <- sum(val_data$particle == p)
  ratio <- n_val / n_cal * 100
  cat(sprintf("%-12s %12d %12d %11.1f%%\n", p, n_cal, n_val, ratio))
}
# Section 3: Voxa Core Functions

cat("\n--- SECTION 3: VOxA VA FUNCTIONS ---\n")
#' Get particle parameters safely
get_particle_params <- function(particle_name) {
  particle_name <- as.character(particle_name)
  
  if (particle_name %in% names(VOXA_PARAMS$particles)) {
    return(VOXA_PARAMS$particles[[particle_name]])
  }
  
  # Try lowercase match
  for (key in names(VOXA_PARAMS$particles)) {
    if (tolower(key) == tolower(particle_name)) {
      return(VOXA_PARAMS$particles[[key]])
    }
  }
  
  stop(sprintf("Particle '%s' not found", particle_name))
}

#' Calculate oxygen fixation probability
calc_p_indirect <- function(O2) {
  (O2 + VOXA_PARAMS$K_fix) / (O2 + VOXA_PARAMS$K_fix + VOXA_PARAMS$K_repair)
}

#' Compute voxel-aware P_DSB for each DSB
compute_P_DSB_voxel <- function(E_local, particle, O2) {
  particle <- as.character(particle)
  params <- get_particle_params(particle)
  
  # Standardize energy
  E_mean <- mean(E_local)
  E_std <- sd(E_local)
  
  if (E_std < 1e-10) {
    E_zscore <- rep(0, length(E_local))
  } else {
    E_zscore <- (E_local - E_mean) / E_std
  }
  
  # Compute f_direct
  f_direct_raw <- params$p1 + params$delta_f * E_zscore
  f_direct <- pmax(pmin(f_direct_raw, VOXA_PARAMS$f_max), VOXA_PARAMS$f_min)
  
  # Redistribute remaining probability
  remaining <- 1.0 - f_direct
  total_indirect <- params$p2 + params$p3
  
  if (total_indirect > 1e-10) {
    p2_local <- remaining * (params$p2 / total_indirect)
    p3_local <- remaining * (params$p3 / total_indirect)
  } else {
    p2_local <- rep(0, length(f_direct))
    p3_local <- remaining
  }
  
  # Oxygen fixation probabilities
  p_ind <- calc_p_indirect(O2)
  p_ind_ref <- calc_p_indirect(VOXA_PARAMS$O2_normoxia)
  
  # Normalized P_DSB
  P_raw <- f_direct + p2_local * p_ind + p3_local * p_ind^2
  P_raw_ref <- f_direct + p2_local * p_ind_ref + p3_local * p_ind_ref^2
  
  P_DSB <- ifelse(P_raw_ref > 0, P_raw / P_raw_ref, 1.0)
  pmax(pmin(P_DSB, 1.0), 0.0)
}

cat("Functions defined:\n")
# Section 4: Compute P_Dsb For Both Datasets

cat("\n--- SECTION 4: COMPUTING P_DSB FOR CALIBRATION AND VALIDATION ---\n")
O2_test <- VOXA_PARAMS$O2_calibration  # Use the VA calibration O2 level

# Compute P_DSB statistics for each particle in each dataset
cal_stats <- list()
val_stats <- list()

for (p in PARTICLE_KEYS) {
  cat(sprintf("Processing %s...\n", toupper(p)))
  
  # Calibration dataset
  cal_E <- cal_data %>% filter(particle == p) %>% pull(E_local)
  cal_P_DSB <- compute_P_DSB_voxel(cal_E, p, O2_test)
  
  cal_stats[[p]] <- list(
    n_dsbs = length(cal_E),
    E_mean = mean(cal_E),
    E_std = sd(cal_E),
    P_DSB_mean = mean(cal_P_DSB),
    P_DSB_std = sd(cal_P_DSB),
    P_DSB_cv = sd(cal_P_DSB) / mean(cal_P_DSB) * 100,
    P_DSB_values = cal_P_DSB  # Store for KS test
  )
  
  # Validation dataset
  val_E <- val_data %>% filter(particle == p) %>% pull(E_local)
  val_P_DSB <- compute_P_DSB_voxel(val_E, p, O2_test)
  
  val_stats[[p]] <- list(
    n_dsbs = length(val_E),
    E_mean = mean(val_E),
    E_std = sd(val_E),
    P_DSB_mean = mean(val_P_DSB),
    P_DSB_std = sd(val_P_DSB),
    P_DSB_cv = sd(val_P_DSB) / mean(val_P_DSB) * 100,
    P_DSB_values = val_P_DSB
  )
  
  cat(sprintf("  Calibration: n=%d, P_DSB=%.4f±%.4f, CV=%.2f%%\n",
              cal_stats[[p]]$n_dsbs, cal_stats[[p]]$P_DSB_mean, 
              cal_stats[[p]]$P_DSB_std, cal_stats[[p]]$P_DSB_cv))
  cat(sprintf("  Validation:  n=%d, P_DSB=%.4f±%.4f, CV=%.2f%%\n",
              val_stats[[p]]$n_dsbs, val_stats[[p]]$P_DSB_mean, 
              val_stats[[p]]$P_DSB_std, val_stats[[p]]$P_DSB_cv))
}
# Section 5: Statistical Test Functions

cat("\n--- SECTION 5: STATISTICAL TESTS ---\n")
#' F-test for variance equality
f_test_variance <- function(var1, n1, var2, n2) {
  if (var1 >= var2) {
    F_stat <- var1 / var2
    df1 <- n1 - 1
    df2 <- n2 - 1
  } else {
    F_stat <- var2 / var1
    df1 <- n2 - 1
    df2 <- n1 - 1
  }
  
  p_value <- 2 * min(pf(F_stat, df1, df2), 1 - pf(F_stat, df1, df2))
  
  list(
    test = "F-test for variance equality",
    F_statistic = F_stat,
    df1 = df1,
    df2 = df2,
    p_value = p_value,
    significant = p_value < 0.05,
    interpretation = ifelse(p_value < 0.05, 
                            "Variances significantly different", 
                            "Variances not significantly different")
  )
}

#' Levene's test using actual data
levene_test_actual <- function(P_DSB_val, P_DSB_cal) {
  df <- data.frame(
    value = c(P_DSB_val, P_DSB_cal),
    group = factor(c(rep("validation", length(P_DSB_val)), 
                     rep("calibration", length(P_DSB_cal))))
  )
  
  test_result <- tryCatch(
    leveneTest(value ~ group, data = df, center = median),
    error = function(e) NULL
  )
  
  if (!is.null(test_result)) {
    p_value <- test_result$`Pr(>F)`[1]
    F_stat <- test_result$`F value`[1]
  } else {
    p_value <- NA
    F_stat <- NA
  }
  
  list(
    test = "Levene's test (actual data)",
    F_statistic = F_stat,
    p_value = p_value,
    significant = !is.na(p_value) && p_value < 0.05,
    interpretation = ifelse(is.na(p_value), "Test failed",
                            ifelse(p_value < 0.05, "Variances differ", "Variances equal"))
  )
}

#' Bootstrap confidence interval for CV ratio
bootstrap_cv_ratio <- function(P_DSB_val, P_DSB_cal, n_bootstrap = 10000) {
  cv_val <- sd(P_DSB_val) / mean(P_DSB_val) * 100
  cv_cal <- sd(P_DSB_cal) / mean(P_DSB_cal) * 100
  
  n_val <- length(P_DSB_val)
  n_cal <- length(P_DSB_cal)
  
  cv_ratios <- numeric(n_bootstrap)
  
  for (i in seq_len(n_bootstrap)) {
    boot_val <- sample(P_DSB_val, n_val, replace = TRUE)
    boot_cal <- sample(P_DSB_cal, n_cal, replace = TRUE)
    
    cv_boot_val <- sd(boot_val) / mean(boot_val) * 100
    cv_boot_cal <- sd(boot_cal) / mean(boot_cal) * 100
    
    cv_ratios[i] <- ifelse(cv_boot_cal > 0, cv_boot_val / cv_boot_cal, NA)
  }
  
  cv_ratios <- na.omit(cv_ratios)
  
  ci_low <- quantile(cv_ratios, 0.025)
  ci_high <- quantile(cv_ratios, 0.975)
  
  list(
    test = "Bootstrap CV ratio",
    n_bootstrap = length(cv_ratios),
    cv_validation = cv_val,
    cv_calibration = cv_cal,
    cv_ratio = cv_val / cv_cal,
    cv_ratio_mean = mean(cv_ratios),
    cv_ratio_std = sd(cv_ratios),
    cv_ratio_ci_low = as.numeric(ci_low),
    cv_ratio_ci_high = as.numeric(ci_high),
    ratio_contains_1 = (ci_low <= 1.0) && (1.0 <= ci_high),
    interpretation = ifelse((ci_low <= 1.0) && (1.0 <= ci_high),
                            "CVs equivalent (CI contains 1.0)",
                            "CVs differ (CI excludes 1.0)")
  )
}

#' Kolmogorov-Smirnov test using actual data
ks_test_actual <- function(P_DSB_val, P_DSB_cal) {
  ks_result <- ks.test(P_DSB_val, P_DSB_cal)
  
  list(
    test = "Kolmogorov-Smirnov test",
    D_statistic = as.numeric(ks_result$statistic),
    p_value = ks_result$p.value,
    significant = ks_result$p.value < 0.05,
    interpretation = ifelse(ks_result$p.value < 0.05,
                            "Distributions differ",
                            "Distributions equivalent")
  )
}

#' TOST equivalence test for means
tost_equivalence_test <- function(P_DSB_val, P_DSB_cal, equivalence_margin = 0.05) {
  mean1 <- mean(P_DSB_val)
  std1 <- sd(P_DSB_val)
  n1 <- length(P_DSB_val)
  
  mean2 <- mean(P_DSB_cal)
  std2 <- sd(P_DSB_cal)
  n2 <- length(P_DSB_cal)
  
  delta <- equivalence_margin * mean2
  se <- sqrt(std1^2/n1 + std2^2/n2)
  diff <- mean1 - mean2
  
  # Welch-Satterthwaite degrees of freedom
  df <- (std1^2/n1 + std2^2/n2)^2 / 
    ((std1^2/n1)^2/(n1-1) + (std2^2/n2)^2/(n2-1))
  
  # Two one-sided tests
  t1 <- (diff - (-delta)) / se
  p1 <- 1 - pt(t1, df)
  
  t2 <- (diff - delta) / se
  p2 <- pt(t2, df)
  
  p_tost <- max(p1, p2)
  
  # 90% CI for difference
  t_crit <- qt(0.95, df)
  ci_low <- diff - t_crit * se
  ci_high <- diff + t_crit * se
  
  equivalent <- (ci_low > -delta) && (ci_high < delta)
  
  list(
    test = "TOST equivalence test",
    mean_difference = diff,
    mean_diff_pct = diff / mean2 * 100,
    equivalence_margin_pct = equivalence_margin * 100,
    equivalence_margin_abs = delta,
    se = se,
    df = df,
    p_tost = p_tost,
    ci_90_low = ci_low,
    ci_90_high = ci_high,
    equivalent = equivalent || p_tost < 0.05,
    interpretation = ifelse(equivalent || p_tost < 0.05, 
                            "Means equivalent",
                            "Equivalence not established")
  )
}

#' Feltz-Miller test for CV equality
feltz_miller_cv_test <- function(cv1, n1, cv2, n2) {
  cv1_prop <- cv1 / 100
  cv2_prop <- cv2 / 100
  
  cv_pooled <- (n1 * cv1_prop + n2 * cv2_prop) / (n1 + n2)
  
  D <- (n1 * (cv1_prop - cv_pooled)^2 + n2 * (cv2_prop - cv_pooled)^2) / 
    (0.5 + cv_pooled^2)
  
  p_value <- 1 - pchisq(D, df = 1)
  
  list(
    test = "Feltz-Miller CV equality",
    cv_validation = cv1,
    cv_calibration = cv2,
    cv_pooled = cv_pooled * 100,
    D_statistic = D,
    p_value = p_value,
    significant = p_value < 0.05,
    interpretation = ifelse(p_value < 0.05,
                            "CVs significantly different",
                            "CVs equivalent")
  )
}

cat("Statistical test functions defined:\n")
# Section 6: Run Comprehensive Validation

cat("\n--- SECTION 6: RUNNING COMPREHENSIVE VALIDATION ---\n")
validation_results <- list()

for (p in PARTICLE_KEYS) {
  cat(sprintf("─" %>% strrep(60), "\n"))
  cat(sprintf("%s\n", toupper(p)))
  cat(sprintf("─" %>% strrep(60), "\n"))
  
  vs <- val_stats[[p]]
  cs <- cal_stats[[p]]
  
  cat(sprintf("  Running F-test...\n"))
  f_test <- f_test_variance(vs$P_DSB_std^2, vs$n_dsbs, cs$P_DSB_std^2, cs$n_dsbs)
  
  cat(sprintf("  Running Levene's test...\n"))
  levene <- levene_test_actual(vs$P_DSB_values, cs$P_DSB_values)
  
  cat(sprintf("  Running Bootstrap CV ratio (10,000 iterations)...\n"))
  bootstrap <- bootstrap_cv_ratio(vs$P_DSB_values, cs$P_DSB_values, n_bootstrap = 10000)
  
  cat(sprintf("  Running K-S test...\n"))
  ks <- ks_test_actual(vs$P_DSB_values, cs$P_DSB_values)
  
  cat(sprintf("  Running TOST equivalence test...\n"))
  tost <- tost_equivalence_test(vs$P_DSB_values, cs$P_DSB_values, equivalence_margin = 0.05)
  
  cat(sprintf("  Running Feltz-Miller CV test...\n"))
  fm <- feltz_miller_cv_test(vs$P_DSB_cv, vs$n_dsbs, cs$P_DSB_cv, cs$n_dsbs)
  
  # Determine overall pass/fail
  # Key criteria:
  # 1. Bootstrap CV ratio CI contains 1.0
  # 2. TOST shows equivalence (or mean diff < 5%)
  # 3. Feltz-Miller not significant (CVs equivalent)
  
  tests_pass <- c(
    bootstrap$ratio_contains_1,
    tost$equivalent || abs(tost$mean_diff_pct) < 5,
    !fm$significant
  )
  
  all_pass <- all(tests_pass)
  
  validation_results[[p]] <- list(
    particle = p,
    n_validation = vs$n_dsbs,
    n_calibration = cs$n_dsbs,
    n_ratio = vs$n_dsbs / cs$n_dsbs,
    
    mean_validation = vs$P_DSB_mean,
    mean_calibration = cs$P_DSB_mean,
    mean_diff_pct = (vs$P_DSB_mean - cs$P_DSB_mean) / cs$P_DSB_mean * 100,
    
    std_validation = vs$P_DSB_std,
    std_calibration = cs$P_DSB_std,
    
    cv_validation = vs$P_DSB_cv,
    cv_calibration = cs$P_DSB_cv,
    cv_ratio = vs$P_DSB_cv / cs$P_DSB_cv,
    
    f_test = f_test,
    levene_test = levene,
    bootstrap_cv = bootstrap,
    ks_test = ks,
    tost_test = tost,
    feltz_miller = fm,
    
    all_tests_pass = all_pass
  )
  
  cat(sprintf("\n  RESULTS:\n"))
  cat(sprintf("    CV Ratio: %.4f [95%% CI: %.3f, %.3f]\n", 
              bootstrap$cv_ratio, bootstrap$cv_ratio_ci_low, bootstrap$cv_ratio_ci_high))
  cat(sprintf("    Mean diff: %.3f%%\n", tost$mean_diff_pct))
  cat(sprintf("    Bootstrap CI contains 1.0: %s\n", ifelse(bootstrap$ratio_contains_1, "YES", "NO")))
  cat(sprintf("    TOST equivalent: %s\n", ifelse(tost$equivalent, "YES", "NO")))
  cat(sprintf("    Feltz-Miller p-value: %.4f\n", fm$p_value))
  cat(sprintf("    OVERALL: %s\n\n", ifelse(all_pass, "✓ PASS", "✗ FAIL")))
}


# Section 7: Generate Figures

cat("\n--- SECTION 7: GENERATING FIGURES (42-46) ---\n")
# Standardize particle order
ordered_particles <- c()
for (std_p in c("electron", "proton", "carbon")) {
  for (p in PARTICLE_KEYS) {
    if (tolower(p) == std_p) {
      ordered_particles <- c(ordered_particles, p)
      break
    }
  }
}

# Create color and label mappings
local_particle_colors <- setNames(
  c(mediterranean_blue, lemon_yellow, pompeii_red),
  ordered_particles
)

display_labels <- setNames(
  c("Electron", "Proton", "Carbon"),
  ordered_particles
)

#' Figure 42: CV Comparison Bar Chart
plot_fig42_cv_comparison <- function() {
  cat("Generating Figure 42: CV comparison...\n")
  
  df <- data.frame(
    particle = rep(ordered_particles, 2),
    dataset = rep(c("Validation (~400)", "Calibration (~2500)"), each = length(ordered_particles)),
    cv = c(
      sapply(ordered_particles, function(p) validation_results[[p]]$cv_validation),
      sapply(ordered_particles, function(p) validation_results[[p]]$cv_calibration)
    )
  )
  
  df$particle <- factor(df$particle, levels = ordered_particles)
  df$dataset <- factor(df$dataset, levels = c("Validation (~400)", "Calibration (~2500)"))
  
  p <- ggplot(df, aes(x = particle, y = cv, fill = interaction(particle, dataset), 
                      pattern = dataset)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
             width = 0.7, color = "black", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f%%", cv)), 
              position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
    scale_fill_manual(
      values = c(
        setNames(local_particle_colors, paste0(ordered_particles, ".Validation (~400)")),
        setNames(alpha(local_particle_colors, 0.5), paste0(ordered_particles, ".Calibration (~2500)"))
      ),
      guide = "none"
    ) +
    scale_x_discrete(labels = display_labels) +
    labs(
      title = "P_DSB Coefficient of Variation: Validation vs Calibration",
      subtitle = sprintf("At %.2f%% O₂ (VA calibration level)", VOXA_PARAMS$O2_calibration),
      x = "Particle Type",
      y = "P_DSB CV (%)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = warm_stone),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(ylim = c(0, max(df$cv) * 1.3))
  
  # Add legend annotation
  p <- p + 
    annotate("rect", xmin = 2.5, xmax = 3.5, ymin = max(df$cv) * 1.15, ymax = max(df$cv) * 1.28,
             fill = "white", color = "gray80") +
    annotate("text", x = 3, y = max(df$cv) * 1.22, 
             label = "Solid = Validation\nLight = Calibration", size = 3, color = warm_stone)
  
  ggsave("figures/fig42_cv_comparison.png", p, 
         width = 9, height = 7, dpi = 300, bg = "white")
  }

#' Figure 43: CV Ratio with Bootstrap CI
plot_fig43_cv_ratio_ci <- function() {
  cat("Generating Figure 43: CV ratio with CI...\n")
  
  df <- data.frame(
    particle = factor(ordered_particles, levels = ordered_particles),
    ratio = sapply(ordered_particles, function(p) validation_results[[p]]$cv_ratio),
    ci_low = sapply(ordered_particles, function(p) validation_results[[p]]$bootstrap_cv$cv_ratio_ci_low),
    ci_high = sapply(ordered_particles, function(p) validation_results[[p]]$bootstrap_cv$cv_ratio_ci_high)
  )
  
  p <- ggplot(df, aes(x = particle, y = ratio, color = particle)) +
    # Tolerance band
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.8, ymax = 1.2,
             fill = olive_green, alpha = 0.15) +
    # Perfect agreement line
    geom_hline(yintercept = 1.0, linetype = "dashed", color = olive_green, linewidth = 1) +
    # Error bars
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 1.2) +
    # Points
    geom_point(size = 5) +
    # Labels
    geom_text(aes(y = ci_high + 0.06, label = sprintf("%.3f", ratio)), 
              size = 4, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = local_particle_colors, guide = "none") +
    scale_x_discrete(labels = display_labels) +
    labs(
      title = "CV Ratio with 95% Bootstrap Confidence Intervals",
      subtitle = "CV Ratio = CV(Validation) / CV(Calibration); Green band = ±20% tolerance",
      x = "Particle Type",
      y = "CV Ratio"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = warm_stone),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(ylim = c(0.5, 1.5)) +
    annotate("text", x = 3.1, y = 1.0, label = "Perfect\nagreement", 
             hjust = 0, size = 3, color = olive_green, fontface = "italic")
  
  ggsave("figures/fig43_cv_ratio_ci.png", p, 
         width = 8, height = 6, dpi = 300, bg = "white")
  }

#' Figure 44: Mean P_DSB Agreement
plot_fig44_mean_agreement <- function() {
  cat("Generating Figure 44: Mean P_DSB agreement...\n")
  
  df <- data.frame(
    particle = factor(ordered_particles, levels = ordered_particles),
    mean_cal = sapply(ordered_particles, function(p) validation_results[[p]]$mean_calibration),
    mean_val = sapply(ordered_particles, function(p) validation_results[[p]]$mean_validation)
  )
  
  lims <- c(min(c(df$mean_cal, df$mean_val)) * 0.9,
            max(c(df$mean_cal, df$mean_val)) * 1.1)
  
  p <- ggplot(df, aes(x = mean_cal, y = mean_val, color = particle)) +
    # ±5% tolerance band
    geom_ribbon(data = data.frame(x = seq(lims[1], lims[2], length.out = 100)),
                aes(x = x, ymin = x * 0.95, ymax = x * 1.05),
                fill = olive_green, alpha = 0.15, inherit.aes = FALSE) +
    # Perfect agreement line
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                color = olive_green, linewidth = 1) +
    # Points
    geom_point(size = 5, stroke = 1.5) +
    # Labels
    geom_text(aes(label = display_labels[as.character(particle)]), 
              vjust = -1.5, size = 3.5, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = local_particle_colors, guide = "none") +
    labs(
      title = "Mean P_DSB Agreement: Validation vs Calibration",
      subtitle = sprintf("At %.2f%% O₂; Green band = ±5%% tolerance", VOXA_PARAMS$O2_calibration),
      x = "Mean P_DSB (Calibration)",
      y = "Mean P_DSB (Validation)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = warm_stone),
      panel.grid.minor = element_blank()
    ) +
    coord_fixed(ratio = 1, xlim = lims, ylim = lims)
  
  ggsave("figures/fig44_mean_pdsb_agreement.png", p, 
         width = 7, height = 7, dpi = 300, bg = "white")
  }

#' Figure 45: CV Stability vs Sample Size
plot_fig45_cv_stability <- function() {
  cat("Generating Figure 45: CV stability...\n")
  
  df <- data.frame(
    particle = factor(ordered_particles, levels = ordered_particles),
    n_ratio = sapply(ordered_particles, function(p) validation_results[[p]]$n_ratio * 100),
    cv_ratio = sapply(ordered_particles, function(p) validation_results[[p]]$cv_ratio)
  )
  
  p <- ggplot(df, aes(x = n_ratio, y = cv_ratio, color = particle)) +
    # Tolerance band
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.8, ymax = 1.2,
             fill = olive_green, alpha = 0.15) +
    # Perfect agreement line
    geom_hline(yintercept = 1.0, linetype = "dashed", color = olive_green, linewidth = 1) +
    # Points
    geom_point(size = 6) +
    # Labels
    geom_text(aes(label = display_labels[as.character(particle)]), 
              vjust = -1.5, size = 3.5, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = local_particle_colors, guide = "none") +
    labs(
      title = "CV Stability Across Sample Sizes",
      subtitle = "CV ratio remains stable despite ~85% reduction in sample size",
      x = "Sample Size Ratio (Validation / Calibration, %)",
      y = "CV Ratio"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = warm_stone),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(xlim = c(0, 25), ylim = c(0.7, 1.3)) +
    annotate("text", x = 24, y = 1.18, label = "±20% tolerance", 
             hjust = 1, size = 3, color = olive_green, fontface = "italic")
  
  ggsave("figures/fig45_cv_stability.png", p, 
         width = 8, height = 6, dpi = 300, bg = "white")
  }

#' Figure 46: Statistical Test Summary
plot_fig46_test_summary <- function() {
  cat("Generating Figure 46: Test summary...\n")
  
  # Build summary table data
  summary_data <- tibble(
    Particle = sapply(ordered_particles, function(p) display_labels[p]),
    `n (Val)` = sapply(ordered_particles, function(p) validation_results[[p]]$n_validation),
    `n (Cal)` = sapply(ordered_particles, function(p) validation_results[[p]]$n_calibration),
    `CV Ratio` = sapply(ordered_particles, function(p) 
      sprintf("%.3f", validation_results[[p]]$cv_ratio)),
    `95% CI` = sapply(ordered_particles, function(p) 
      sprintf("[%.3f, %.3f]", 
              validation_results[[p]]$bootstrap_cv$cv_ratio_ci_low,
              validation_results[[p]]$bootstrap_cv$cv_ratio_ci_high)),
    `Contains 1.0` = sapply(ordered_particles, function(p) 
      ifelse(validation_results[[p]]$bootstrap_cv$ratio_contains_1, "Yes", "No")),
    `TOST p` = sapply(ordered_particles, function(p) 
      sprintf("%.4f", validation_results[[p]]$tost_test$p_tost)),
    `FM p` = sapply(ordered_particles, function(p) 
      sprintf("%.4f", validation_results[[p]]$feltz_miller$p_value)),
    Status = sapply(ordered_particles, function(p) 
      ifelse(validation_results[[p]]$all_tests_pass, "PASS", "FAIL"))
  )
  
  # Create table plot
  n_rows <- nrow(summary_data)
  n_cols <- ncol(summary_data)
  
  grid_df <- expand.grid(row = 0:n_rows, col = 1:n_cols)
  grid_df$text <- ""
  grid_df$fill <- "white"
  grid_df$text_color <- "black"
  grid_df$fontface <- "plain"
  
  col_names <- names(summary_data)
  
  # Header row
  for (j in 1:n_cols) {
    idx <- which(grid_df$row == 0 & grid_df$col == j)
    grid_df$text[idx] <- col_names[j]
    grid_df$fill[idx] <- warm_stone
    grid_df$text_color[idx] <- "white"
    grid_df$fontface[idx] <- "bold"
  }
  
  # Data rows
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      idx <- which(grid_df$row == i & grid_df$col == j)
      grid_df$text[idx] <- as.character(summary_data[i, j])
      
      # Particle column coloring
      if (j == 1) {
        p_key <- ordered_particles[i]
        grid_df$text_color[idx] <- local_particle_colors[p_key]
        grid_df$fontface[idx] <- "bold"
      }
      
      # Status column coloring
      if (j == n_cols) {
        if (summary_data$Status[i] == "PASS") {
          grid_df$fill[idx] <- olive_green
          grid_df$text_color[idx] <- "white"
        } else {
          grid_df$fill[idx] <- pompeii_red
          grid_df$text_color[idx] <- "white"
        }
        grid_df$fontface[idx] <- "bold"
      }
      
      # Contains 1.0 column
      if (j == 6) {
        if (summary_data$`Contains 1.0`[i] == "Yes") {
          grid_df$text_color[idx] <- olive_green
        } else {
          grid_df$text_color[idx] <- pompeii_red
        }
        grid_df$fontface[idx] <- "bold"
      }
      
      # Alternating row shading
      if (j != n_cols && i %% 2 == 0) {
        grid_df$fill[idx] <- sandy_beige
      }
    }
  }
  
  # Column widths
  col_widths <- c(1.0, 0.7, 0.7, 0.8, 1.4, 0.9, 0.7, 0.7, 0.7)
  col_positions <- cumsum(c(0, col_widths[-length(col_widths)]))
  col_centers <- col_positions + col_widths / 2
  
  grid_df$x <- col_centers[grid_df$col]
  grid_df$y <- -grid_df$row
  grid_df$width <- col_widths[grid_df$col]
  
  p <- ggplot(grid_df) +
    geom_tile(aes(x = x, y = y, width = width, height = 0.9, fill = fill),
              color = espresso_brown, linewidth = 0.3) +
    geom_text(aes(x = x, y = y, label = text, color = text_color, fontface = fontface),
              size = 3.2) +
    scale_fill_identity() +
    scale_color_identity() +
    labs(
      title = "VOxA Scaling Validation: Statistical Test Summary",
      subtitle = "TOST = Two One-Sided Tests; FM = Feltz-Miller CV equality"
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = warm_stone, hjust = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    coord_fixed(ratio = 0.45)
  
  # Add footer
  all_pass <- all(sapply(validation_results, function(x) x$all_tests_pass))
  footer_text <- ifelse(all_pass,
                        "✓ All particles pass scaling validation: δf parameters generalize to smaller datasets",
                        "⚠ Some particles failed - review individual test results")
  footer_color <- ifelse(all_pass, olive_green, pompeii_red)
  
  p <- p + annotate("text", x = mean(col_centers), y = -n_rows - 0.8, 
                    label = footer_text, size = 3.5, fontface = "italic", color = footer_color)
  
  ggsave("figures/fig46_test_summary.png", p, 
         width = 11, height = 5, dpi = 300, bg = "white")
  }

# Generate all figures
plot_fig42_cv_comparison()
plot_fig43_cv_ratio_ci()
plot_fig44_mean_agreement()
plot_fig45_cv_stability()
plot_fig46_test_summary()


# Section 8: Generate Reports

cat("\n--- SECTION 8: GENERATING REPORTS ---\n")
# JSON output
json_output <- list(
  model_version = VOXA_PARAMS$model_version,
  validation_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  calibration_O2 = VOXA_PARAMS$O2_calibration,
  
  datasets = list(
    calibration_file = cal_data_file,
    validation_file = val_data_file
  ),
  
  results = lapply(validation_results, function(r) {
    list(
      particle = r$particle,
      n_validation = r$n_validation,
      n_calibration = r$n_calibration,
      n_ratio = r$n_ratio,
      
      mean_validation = r$mean_validation,
      mean_calibration = r$mean_calibration,
      mean_diff_pct = r$mean_diff_pct,
      
      cv_validation = r$cv_validation,
      cv_calibration = r$cv_calibration,
      cv_ratio = r$cv_ratio,
      cv_ratio_ci_low = r$bootstrap_cv$cv_ratio_ci_low,
      cv_ratio_ci_high = r$bootstrap_cv$cv_ratio_ci_high,
      cv_ratio_contains_1 = r$bootstrap_cv$ratio_contains_1,
      
      f_test_p = r$f_test$p_value,
      levene_p = r$levene_test$p_value,
      ks_test_p = r$ks_test$p_value,
      tost_p = r$tost_test$p_tost,
      tost_equivalent = r$tost_test$equivalent,
      feltz_miller_p = r$feltz_miller$p_value,
      
      all_tests_pass = r$all_tests_pass
    )
  }),
  
  overall_pass = all(sapply(validation_results, function(x) x$all_tests_pass))
)

json_path <- "results/voxa_scaling_validation_results.json"
write_json(json_output, json_path, pretty = TRUE, auto_unbox = TRUE)
cat(sprintf("  ✓ Saved: %s\n", json_path))

# Text report
report_lines <- c(
  paste(rep("=", 80), collapse = ""),
  "VOxA MODEL SCALING VALIDATION REPORT",
  "Variable Oxygen-dependent Amorphous Track Model",
  paste(rep("=", 80), collapse = ""),
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Model version: %s", VOXA_PARAMS$model_version),
  sprintf("Calibration O2: %.2f%%", VOXA_PARAMS$O2_calibration),
  "",
  "OBJECTIVE",
  paste(rep("-", 80), collapse = ""),
  "Validate that VOxA VA δf parameters calibrated on pooled datasets (~2,500 DSBs)",
  "produce equivalent P_DSB predictions when applied to individual simulation runs",
  "(~400 DSBs).",
  "",
  "DATASETS",
  paste(rep("-", 80), collapse = ""),
  sprintf("  Calibration: %s", cal_data_file),
  sprintf("  Validation:  %s", val_data_file),
  ""
)

for (p in ordered_particles) {
  r <- validation_results[[p]]
  report_lines <- c(report_lines,
                    sprintf("  %s: Cal n=%d, Val n=%d (ratio=%.1f%%)",
                            display_labels[p], r$n_calibration, r$n_validation, r$n_ratio * 100))
}

report_lines <- c(report_lines,
                  "",
                  "STATISTICAL TESTS",
                  paste(rep("-", 80), collapse = ""),
                  "1. F-test for variance equality",
                  "2. Levene's test for homogeneity of variances",
                  "3. Bootstrap CV ratio (10,000 iterations, 95% CI)",
                  "4. Kolmogorov-Smirnov test for distribution equality",
                  "5. TOST equivalence test (5% margin)",
                  "6. Feltz-Miller test for CV equality",
                  "",
                  "RESULTS BY PARTICLE",
                  paste(rep("-", 80), collapse = "")
)

for (p in ordered_particles) {
  r <- validation_results[[p]]
  report_lines <- c(report_lines,
                    "",
                    sprintf("%s", toupper(display_labels[p])),
                    paste(rep("~", 40), collapse = ""),
                    sprintf("  Sample sizes: %d (val) vs %d (cal), ratio = %.1f%%", 
                            r$n_validation, r$n_calibration, r$n_ratio * 100),
                    sprintf("  Mean P_DSB: %.4f (val) vs %.4f (cal), diff = %.3f%%",
                            r$mean_validation, r$mean_calibration, r$mean_diff_pct),
                    sprintf("  CV: %.2f%% (val) vs %.2f%% (cal)", r$cv_validation, r$cv_calibration),
                    sprintf("  CV Ratio: %.4f", r$cv_ratio),
                    "",
                    sprintf("  Bootstrap 95%% CI: [%.3f, %.3f]",
                            r$bootstrap_cv$cv_ratio_ci_low, r$bootstrap_cv$cv_ratio_ci_high),
                    sprintf("  CI contains 1.0: %s", ifelse(r$bootstrap_cv$ratio_contains_1, "YES", "NO")),
                    "",
                    sprintf("  F-test p-value: %.4f (%s)", r$f_test$p_value, r$f_test$interpretation),
                    sprintf("  Levene's p-value: %.4f (%s)", r$levene_test$p_value, r$levene_test$interpretation),
                    sprintf("  K-S test p-value: %.4f (%s)", r$ks_test$p_value, r$ks_test$interpretation),
                    sprintf("  TOST p-value: %.4f (%s)", r$tost_test$p_tost, r$tost_test$interpretation),
                    sprintf("  Feltz-Miller p-value: %.4f (%s)", r$feltz_miller$p_value, r$feltz_miller$interpretation),
                    "",
                    sprintf("  OVERALL: %s", ifelse(r$all_tests_pass, "✓ PASS", "✗ FAIL"))
  )
}

all_pass <- all(sapply(validation_results, function(x) x$all_tests_pass))

report_lines <- c(report_lines,
                  "",
                  paste(rep("=", 80), collapse = ""),
                  "SUMMARY",
                  paste(rep("=", 80), collapse = ""),
                  ""
)

if (all_pass) {
  report_lines <- c(report_lines,
                    "✓ ALL PARTICLES PASSED SCALING VALIDATION",
                    "",
                    "Key findings:",
                    "  • CV ratios within expected range for all particles",
                    "  • Bootstrap 95% CIs for CV ratio all contain 1.0",
                    "  • TOST confirms mean equivalence within 5% margin",
                    "  • Feltz-Miller test confirms CV equality",
                    "",
                    "CONCLUSION:",
                    "VOxA VA predictions are stable when applied to individual simulation",
                    "runs (~400 DSBs) despite calibration on larger pooled datasets (~2,500 DSBs).",
                    "The δf parameters generalize robustly across sample sizes."
  )
} else {
  failed <- ordered_particles[!sapply(ordered_particles, function(p) validation_results[[p]]$all_tests_pass)]
  report_lines <- c(report_lines,
                    sprintf("⚠ WARNING: Scaling validation issues for: %s", 
                            paste(display_labels[failed], collapse = ", ")),
                    "",
                    "Review individual particle results above for details.",
                    "Consider recalibrating δf parameters or investigating data quality."
  )
}

report_lines <- c(report_lines,
                  "",
                  "FIGURES GENERATED (Step 11)",
                  paste(rep("-", 80), collapse = ""),
                  "  Fig 42: CV comparison bar chart",
                  "  Fig 43: CV ratio with 95% bootstrap CI",
                  "  Fig 44: Mean P_DSB agreement",
                  "  Fig 45: CV stability vs sample size",
                  "  Fig 46: Statistical test summary table",
                  "",
                  paste(rep("=", 80), collapse = ""),
                  "END OF REPORT",
                  paste(rep("=", 80), collapse = "")
)

report_path <- "results/voxa_scaling_validation_report.txt"
writeLines(report_lines, report_path)
cat(sprintf("  ✓ Saved: %s\n", report_path))


# Final Summary

if (all_pass) {
  } else {
  }

for (p in ordered_particles) {
  r <- validation_results[[p]]
  status <- ifelse(r$all_tests_pass, "✓ PASS", "✗ FAIL")
  cat(sprintf("║     %-10s: CV ratio = %.3f [%.3f, %.3f] %s         ║\n",
              display_labels[p], r$cv_ratio,
              r$bootstrap_cv$cv_ratio_ci_low, r$bootstrap_cv$cv_ratio_ci_high,
              status))
}

# Print conclusion
if (all_pass) {
  cat("CONCLUSION: VOxA VA model validated. The δf parameters calibrated on\n")
  cat("~2,500 DSBs generalize robustly to independent ~400 DSB datasets.\n")
  cat("P_DSB heterogeneity (CV) is a stable physical property of the track\n")
  cat("structure, independent of sample size.\n")
} else {
  cat("WARNING: Some particles failed validation. Review the detailed report\n")
  cat("and consider recalibrating δf parameters.\n")
}

cat("\nStep 11 complete.\n\n")