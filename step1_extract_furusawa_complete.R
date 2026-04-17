# Step 1: Furusawa et al. (2000) data extraction
#
# Extracts paired aerobic/hypoxic D10 values for He, C, and Ne ions from
# Furusawa 2000 Tables 2-4. Only experiments with both conditions are kept.
#
# OER = D10_hypoxic / D10_aerobic (dose-ratio definition at 10% survival).
# For high-LET hypoxic data where beta≈0, D10_hyp = 2.303 / alpha_hyp.
#
# Output: data/furusawa_oer_data.csv

library(tidyverse)

# Try to load setup, create if not exists
if (file.exists("data/uvaom_recalibration_setup.RData")) {
  load("data/uvaom_recalibration_setup.RData")
} else {
  # Define conversion factor if setup not available
  CONVERSION_FACTOR <- list(
    mean = 1.20,
    sd = 0.05,
    range = c(1.15, 1.25),
    source = "Hirayama2005"
  )
}

# Complete Furusawa Data With D10 Values
# 
# Extracted from Tables 2-4 of Furusawa et al. 2000
# Only paired experiments (both aerobic AND hypoxic) are included
# 
# Note: For hypoxic conditions, the original table's "D10" column appears
# to contain α values. We extract:
#   - D10_aer: directly from table (aerobic D10)
#   - alpha_hyp: from hypoxic α column
#   - D10_hyp: calculated as 2.303/alpha_hyp (valid for high LET where β≈0)
#

furusawa_paired <- tribble(
  ~ion, ~LET, ~Z_star2_beta2, ~cell_line, ~exp_aer, ~alpha_aer, ~D10_aer, ~exp_hyp, ~alpha_hyp,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # ³He IONS - V79 CELLS (18 paired experiments)
  # ═══════════════════════════════════════════════════════════════════════════
  "He", 18.6, 196, "V79", "V79-304", 0.225, 5.17, "V79-305", 0.099,
  "He", 18.6, 196, "V79", "V79-324", 0.297, 5.04, "V79-325", 0.087,
  "He", 23.0, 255, "V79", "V79-198", 0.255, 5.04, "V79-199", 0.165,
  "He", 23.8, 266, "V79", "V79-306", 0.261, 4.47, "V79-307", 0.079,
  "He", 24.0, 269, "V79", "V79-336", 0.168, 4.98, "V79-337", 0.108,
  "He", 29.9, 354, "V79", "V79-308", 0.162, 4.76, "V79-309", 0.062,
  "He", 29.9, 354, "V79", "V79-326", 0.376, 4.13, "V79-327", 0.157,
  "He", 38.1, 482, "V79", "V79-200", 0.534, 3.25, "V79-201", 0.197,
  "He", 39.2, 500, "V79", "V79-310", 0.529, 3.36, "V79-311", 0.132,
  "He", 39.4, 504, "V79", "V79-328", 0.486, 3.53, "V79-329", 0.134,
  "He", 50.0, 692, "V79", "V79-338", 0.561, 2.70, "V79-339", 0.153,
  "He", 51.9, 729, "V79", "V79-312", 0.405, 2.91, "V79-313", 0.177,
  "He", 61.9, 922, "V79", "V79-314", 0.560, 2.47, "V79-315", 0.327,
  "He", 73.9, 1179, "V79", "V79-340", 1.040, 1.93, "V79-341", 0.329,
  "He", 74.6, 1195, "V79", "V79-316", 1.105, 1.71, "V79-317", 0.271,
  "He", 74.6, 1195, "V79", "V79-332", 0.822, 1.95, "V79-333", 0.303,
  "He", 90.8, 1578, "V79", "V79-318", 0.956, 1.79, "V79-319", 0.541,
  "He", 90.8, 1578, "V79", "V79-334", 1.117, 1.77, "V79-335", 0.565,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # ³He IONS - HSG CELLS (6 paired experiments)
  # ═══════════════════════════════════════════════════════════════════════════
  "He", 18.5, 195, "HSG", "HSG-64", 0.543, 3.12, "HSG-65", 0.213,
  "He", 23.6, 263, "HSG", "HSG-66", 0.662, 2.71, "HSG-67", 0.234,
  "He", 32.4, 392, "HSG", "HSG-68", 1.221, 1.89, "HSG-69", 0.442,
  "He", 46.0, 617, "HSG", "HSG-70", 1.328, 1.55, "HSG-71", 0.557,
  "He", 53.8, 765, "HSG", "HSG-72", 1.705, 1.30, "HSG-73", 0.818,
  "He", 70.3, 1100, "HSG", "HSG-74", 2.283, 1.01, "HSG-75", 1.248,
  
  # ════════════════════════════════��══════════════════════════════════════════
  # ¹²C IONS - V79 CELLS (19 paired experiments)
  # ═══════════════════════════════════════════════════════════════════════════
  "C", 30.0, 239, "V79", "V79-387", 0.295, 4.80, "V79-388", 0.072,
  "C", 31.0, 248, "V79", "V79-253", 0.427, 4.07, "V79-254", 0.077,
  "C", 40.1, 337, "V79", "V79-389", 0.283, 3.01, "V79-390", 0.044,
  "C", 40.6, 342, "V79", "V79-431", 0.619, 3.28, "V79-432", 0.141,
  "C", 50.3, 446, "V79", "V79-466", 0.841, 2.56, "V79-467", 0.231,
  "C", 57.6, 520, "V79", "V79-231", 0.549, 2.99, "V79-232", 0.187,
  "C", 60.0, 544, "V79", "V79-255", 0.631, 2.65, "V79-256", 0.218,
  "C", 78.5, 736, "V79", "V79-433", 1.027, 2.15, "V79-434", 0.520,
  "C", 80.6, 756, "V79", "V79-393", 0.781, 2.26, "V79-394", 0.434,
  "C", 88.0, 836, "V79", "V79-395", 0.971, 2.11, "V79-396", 0.449,
  "C", 102, 989, "V79", "V79-233", 1.066, 1.95, "V79-234", 0.390,
  "C", 117, 1160, "V79", "V79-257", 1.315, 1.69, "V79-258", 0.543,
  "C", 127, 1280, "V79", "V79-235", 1.080, 2.05, "V79-236", 0.646,
  "C", 142, 1460, "V79", "V79-468", 1.446, 1.59, "V79-469", 0.915,
  "C", 206, 2280, "V79", "V79-435", 1.323, 1.74, "V79-436", 1.198,
  "C", 232, 2640, "V79", "V79-298", 1.064, 2.16, "V79-299", 1.130,
  "C", 255, 2980, "V79", "V79-272", 1.070, 2.15, "V79-273", 0.979,
  "C", 276, 3290, "V79", "V79-223", 1.091, 1.95, "V79-224", 1.044,
  "C", 360, 4670, "V79", "V79-446", 1.070, 2.15, "V79-447", 0.944,
  "C", 432, 6060, "V79", "V79-225", 0.844, 2.73, "V79-226", 0.755,
  "C", 493, 7300, "V79", "V79-444", 0.837, 2.75, "V79-445", 0.793,
  "C", 502, 7500, "V79", "V79-448", 0.779, 2.96, "V79-449", 0.792,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # ¹²C IONS - HSG CELLS (16 paired experiments)
  # ═══════════════════════════════════════════════════════════════════════════
  "C", 22.5, 169, "HSG", "HSG-132", 0.161, 3.71, "HSG-133", 0.115,
  "C", 30.0, 239, "HSG", "HSG-134", 0.429, 3.30, "HSG-135", 0.202,
  "C", 30.3, 242, "HSG", "HSG-212", 0.359, 3.22, "HSG-213", 0.094,
  "C", 30.6, 244, "HSG", "HSG-168", 0.483, 2.88, "HSG-169", 0.137,
  "C", 40.2, 338, "HSG", "HSG-136", 0.720, 2.25, "HSG-137", 0.181,
  "C", 42.5, 361, "HSG", "HSG-170", 0.805, 2.25, "HSG-171", 0.168,
  "C", 50.0, 443, "HSG", "HSG-105", 1.137, 1.86, "HSG-106", 0.418,
  "C", 54.5, 489, "HSG", "HSG-138", 1.010, 1.84, "HSG-139", 0.353,
  "C", 80.6, 758, "HSG", "HSG-140", 1.449, 1.56, "HSG-141", 0.596,
  "C", 88.0, 836, "HSG", "HSG-142", 1.680, 1.36, "HSG-143", 0.761,
  "C", 137, 1400, "HSG", "HSG-144", 1.900, 1.21, "HSG-145", 1.168,
  "C", 144, 1490, "HSG", "HSG-216", 1.507, 1.53, "HSG-217", 1.049,
  "C", 147, 1530, "HSG", "HSG-146", 1.742, 1.32, "HSG-147", 1.147,
  "C", 199, 2190, "HSG", "HSG-107", 1.712, 1.34, "HSG-108", 1.558,
  "C", 247, 2860, "HSG", "HSG-148", 1.423, 1.62, "HSG-149", 1.114,
  "C", 360, 4670, "HSG", "HSG-194", 1.245, 1.85, "HSG-195", 1.199,
  "C", 467, 6760, "HSG", "HSG-174", 1.075, 2.14, "HSG-175", 0.874,
  "C", 493, 7300, "HSG", "HSG-188", 1.035, 2.22, "HSG-189", 0.980,
  "C", 502, 7500, "HSG", "HSG-196", 0.923, 2.50, "HSG-197", 0.826,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # ²⁰Ne IONS - V79 CELLS (13 paired experiments)
  # ═══════════════════════════════════════════════════════════════════════════
  "Ne", 62.1, 466, "V79", "V79-259", 0.511, 3.21, "V79-260", 0.168,
  "Ne", 62.2, 467, "V79", "V79-450", 0.634, 3.17, "V79-451", 0.144,
  "Ne", 80.0, 631, "V79", "V79-282", 0.555, 2.87, "V79-283", 0.199,
  "Ne", 84.6, 675, "V79", "V79-452", 0.672, 2.73, "V79-453", 0.237,
  "Ne", 96.9, 793, "V79", "V79-261", 0.761, 2.36, "V79-262", 0.200,
  "Ne", 110, 922, "V79", "V79-284", 0.774, 2.51, "V79-258", 0.381,
  "Ne", 146, 1310, "V79", "V79-408", 1.087, 2.01, "V79-409", 0.666,
  "Ne", 158, 1430, "V79", "V79-286", 1.225, 1.88, "V79-287", 0.821,
  "Ne", 178, 1620, "V79", "V79-263", 1.259, 1.83, "V79-264", 0.880,
  "Ne", 219, 2050, "V79", "V79-288", 1.150, 2.00, "V79-289", 0.896,
  "Ne", 239, 2260, "V79", "V79-456", 1.299, 1.77, "V79-457", 0.986,
  "Ne", 287, 2790, "V79", "V79-410", 1.212, 1.90, "V79-411", 0.969,
  "Ne", 528, 5820, "V79", "V79-412", 0.675, 3.41, "V79-413", 0.676,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # ²⁰Ne IONS - HSG CELLS (20 paired experiments)
  # ═══════���═══════════════════════════════════════════════════════════════════
  "Ne", 62.0, 465, "HSG", "HSG-198", 0.739, 2.47, "HSG-199", 0.267,
  "Ne", 72.0, 556, "HSG", "HSG-89", 1.201, 1.68, "HSG-90", 0.431,
  "Ne", 74.2, 577, "HSG", "HSG-150", 1.135, 1.71, "HSG-151", 0.352,
  "Ne", 82.0, 650, "HSG", "HSG-93", 1.735, 1.30, "HSG-94", 0.495,
  "Ne", 83.7, 666, "HSG", "HSG-117", 1.516, 1.42, "HSG-118", 0.517,
  "Ne", 103, 853, "HSG", "HSG-85", 2.000, 1.15, "HSG-86", 0.734,
  "Ne", 116, 978, "HSG", "HSG-109", 1.878, 1.23, "HSG-110", 0.703,
  "Ne", 166, 1510, "HSG", "HSG-119", 1.947, 1.18, "HSG-153", 0.999,
  "Ne", 169, 1530, "HSG", "HSG-95", 2.009, 1.15, "HSG-96", 1.435,
  "Ne", 216, 2020, "HSG", "HSG-111", 1.849, 1.25, "HSG-112", 1.418,
  "Ne", 287, 2790, "HSG", "HSG-154", 1.338, 1.72, "HSG-155", 1.281,
  "Ne", 339, 3400, "HSG", "HSG-202", 1.253, 1.84, "HSG-203", 1.065,
  "Ne", 340, 3410, "HSG", "HSG-224", 1.211, 1.90, "HSG-225", 1.071,
  "Ne", 343, 3440, "HSG", "HSG-87", 1.253, 1.84, "HSG-88", 1.087,
  "Ne", 347, 3490, "HSG", "HSG-121", 1.315, 1.75, "HSG-122", 1.150,
  "Ne", 361, 3660, "HSG", "HSG-230", 1.068, 2.16, "HSG-231", 0.903,
  "Ne", 373, 3810, "HSG", "HSG-200", 0.998, 2.31, "HSG-201", 0.834,
  "Ne", 445, 4700, "HSG", "HSG-115a", 0.871, 2.64, "HSG-115b", 0.818,
  "Ne", 528, 5820, "HSG", "HSG-154a", 0.829, 2.78, "HSG-157", 0.797,
  "Ne", 654, 7640, "HSG", "HSG-91", 0.700, 3.29, "HSG-92", 0.620
)

cat(sprintf("Raw paired experiments from Furusawa tables: %d\n\n", nrow(furusawa_paired)))


# Calculate D10 And Oer

cat("\n--- CALCULATING D10 AND OER VALUES ---\n")
# Function to estimate D10 from alpha
estimate_D10_from_alpha <- function(alpha, D10_aer_ref = NULL, alpha_aer_ref = NULL) {
  if (is.na(alpha) || alpha <= 0) return(NA)
  
  if (!is.null(D10_aer_ref) && !is.null(alpha_aer_ref) && alpha_aer_ref > 0) {
    beta_aer <- (log(10) - alpha_aer_ref * D10_aer_ref) / (D10_aer_ref^2)
    beta_hyp <- max(0, beta_aer * 0.8)
    
    if (beta_hyp > 0.001) {
      discriminant <- alpha^2 + 4 * beta_hyp * log(10)
      if (discriminant > 0) {
        D10 <- (-alpha + sqrt(discriminant)) / (2 * beta_hyp)
        return(D10)
      }
    }
  }
  
  return(2.303 / alpha)
}

furusawa_oer <- furusawa_paired %>%
  rowwise() %>%
  mutate(
    D10_hyp = estimate_D10_from_alpha(alpha_hyp, D10_aer, alpha_aer),
    OER_survival = D10_hyp / D10_aer,
    OER_alpha_ratio = alpha_aer / alpha_hyp,
    OER_method_diff = abs(OER_survival - OER_alpha_ratio) / OER_survival * 100
  ) %>%
  ungroup() %>%
  filter(!is.na(OER_survival) & OER_survival > 0)

cat(sprintf("After OER calculation: %d points\n", nrow(furusawa_oer)))


# Compare Methods

cat("\nComparison of OER calculation methods:\n")
cat(sprintf("  Mean difference (D10 vs α ratio): %.1f%%\n", 
            mean(furusawa_oer$OER_method_diff, na.rm = TRUE)))
cat(sprintf("  Max difference: %.1f%%\n", 
            max(furusawa_oer$OER_method_diff, na.rm = TRUE)))
cat(sprintf("  Median difference: %.1f%%\n\n", 
            median(furusawa_oer$OER_method_diff, na.rm = TRUE)))


# Convert To Retention Oer

furusawa_oer <- furusawa_oer %>%
  mutate(
    OER_retention = OER_survival * CONVERSION_FACTOR$mean,
    D10_relative_error = 0.10,
    OER_survival_se = OER_survival * sqrt(2) * D10_relative_error,
    OER_retention_se = sqrt(
      (OER_survival * CONVERSION_FACTOR$sd)^2 + 
        (OER_survival_se * CONVERSION_FACTOR$mean)^2
    ),
    weight = 1 / OER_retention_se^2,
    source = "Furusawa2000"
  )


# Outlier Identification

furusawa_oer <- furusawa_oer %>%
  mutate(
    outlier_flag = case_when(
      OER_survival < 0.95 ~ "OER < 1 (physically implausible)",
      OER_survival > 4.0 ~ "OER > 4 (exceeds typical maximum)",
      LET > 200 & OER_survival > 2.0 ~ "High LET with high OER",
      LET < 30 & OER_survival < 1.5 ~ "Low LET with low OER",
      TRUE ~ NA_character_
    ),
    is_outlier = !is.na(outlier_flag)
  )

# Show outliers
cat("Identified outliers:\n")
outliers_df <- furusawa_oer %>%
  filter(is_outlier) %>%
  select(ion, LET, cell_line, OER_survival, OER_retention, outlier_flag)

if (nrow(outliers_df) > 0) {
  print(outliers_df %>% arrange(desc(OER_survival)), n = 30)
} else {
  }

cat(sprintf("\nTotal outliers: %d / %d (%.1f%%)\n",
            sum(furusawa_oer$is_outlier),
            nrow(furusawa_oer),
            100 * sum(furusawa_oer$is_outlier) / nrow(furusawa_oer)))


# Prepare Final Data

furusawa_final <- furusawa_oer %>%
  select(ion, LET, Z_star2_beta2, cell_line, 
         OER_survival, OER_retention, OER_retention_se, 
         weight, source, is_outlier, outlier_flag) %>%
  mutate(
    ion_std = ion,
    LET_log = log10(LET)
  )


# Summary Statistics

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat(sprintf("Total Furusawa data points: %d\n", nrow(furusawa_final)))
cat(sprintf("  Outliers flagged: %d\n", sum(furusawa_final$is_outlier)))
cat(sprintf("  Clean data points: %d\n", sum(!furusawa_final$is_outlier)))

cat(sprintf("\nLET range: %.2f - %.1f keV/μm\n", 
            min(furusawa_final$LET), max(furusawa_final$LET)))
cat(sprintf("OER_survival range: %.2f - %.2f\n", 
            min(furusawa_final$OER_survival, na.rm = TRUE), 
            max(furusawa_final$OER_survival, na.rm = TRUE)))
cat(sprintf("OER_retention range: %.2f - %.2f\n", 
            min(furusawa_final$OER_retention), 
            max(furusawa_final$OER_retention)))

cat("\nBy ion type (all data):\n")
print(furusawa_final %>%
        group_by(ion_std) %>%
        summarise(
          n = n(),
          n_outliers = sum(is_outlier),
          n_clean = sum(!is_outlier),
          LET_min = round(min(LET), 1),
          LET_max = round(max(LET), 1),
          OER_surv_mean = round(mean(OER_survival, na.rm = TRUE), 2),
          .groups = "drop"
        ) %>%
        arrange(match(ion_std, c("He", "C", "Ne"))))

cat("\nBy cell line:\n")
print(furusawa_final %>% 
        group_by(cell_line) %>%
        summarise(
          n = n(), 
          n_outliers = sum(is_outlier),
          n_clean = sum(!is_outlier)
        ) %>%
        arrange(desc(n)))


# Oer Vs Let Trend Check

cat("\n═════════════════════════════════════════���═════════════════════════\n")
let_bins <- furusawa_final %>%
  filter(!is_outlier) %>%
  mutate(LET_bin = cut(LET, breaks = c(0, 30, 60, 100, 200, 400, 1000),
                       labels = c("<30", "30-60", "60-100", "100-200", "200-400", ">400"))) %>%
  group_by(LET_bin) %>%
  summarise(
    n = n(),
    OER_surv_mean = round(mean(OER_survival, na.rm = TRUE), 2),
    OER_surv_sd = round(sd(OER_survival, na.rm = TRUE), 2),
    .groups = "drop"
  )

cat("OER vs LET trend (clean data):\n")
print(let_bins)

cat("\nExpected trend: OER should DECREASE with increasing LET\n")
# Save Data

cat("\n═══════════════════���═══════════════════════════════════════════════\n")
if (!dir.exists("data")) dir.create("data")

# Save all data (including outliers, flagged)
write_csv(furusawa_final, "data/furusawa_oer_complete_corrected.csv")
saveRDS(furusawa_final, "data/furusawa_oer_complete_corrected.rds")

# Save clean data only (outliers removed)
furusawa_clean <- furusawa_final %>% filter(!is_outlier)
write_csv(furusawa_clean, "data/furusawa_oer_clean.csv")
saveRDS(furusawa_clean, "data/furusawa_oer_clean.rds")

cat("Saved:\n")
cat(sprintf("  data/furusawa_oer_complete_corrected.csv (%d points, includes outliers)\n", 
            nrow(furusawa_final)))
cat(sprintf("  data/furusawa_oer_clean.csv (%d points, outliers removed)\n", 
            nrow(furusawa_clean)))


# Final Summary

cat("\n╔══════════════════════════════════════════════════════════════════╗\n")
cat(sprintf("║  Raw paired experiments: %d                                     ║\n", 
            nrow(furusawa_paired)))
cat(sprintf("║  After OER calculation: %d                                      ║\n", 
            nrow(furusawa_oer)))
cat(sprintf("║  Outliers identified: %d                                         ║\n",
            sum(furusawa_final$is_outlier)))
cat(sprintf("║  Clean data points: %d                                          ║\n",
            nrow(furusawa_clean)))
cat("╚══════════════════════════════════════════════════════════════════╝\n")