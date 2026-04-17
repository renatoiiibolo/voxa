# Step 2: Literature OER dataset compilation
#
# Assembles the 233-observation calibration dataset from 29 published sources,
# covering 10 particle types across LET = 0.2–654 keV/µm.
#
# Data are restricted to 10% clonogenic survival endpoint, in vitro conditions,
# and standard cell lines (V79, HSG, T1, CHO). Tinganelli et al. (2015) adds
# 23 intermediate-O2 observations not in the Wenzl & Wilkens (2011) compilation.
#
# OER_retention = OER_survival × 1.20 (Hirayama 2005 conversion factor).
#
# Output: data/literature_oer_data.csv, data/combined_calibration_data.csv

library(tidyverse)

# Load setup (contains CONVERSION_FACTOR)
if (file.exists("data/uvaom_recalibration_setup.RData")) {
  load("data/uvaom_recalibration_setup.RData")
} else {
  CONVERSION_FACTOR <- list(
    mean = 1.20,
    sd = 0.05,
    range = c(1.15, 1.25),
    source = "Hirayama2005"
  )
}

cat("All OER values are at 10% survival (D10-based) unless noted.\n")
cat("OER_retention = OER_survival × 1.20 (Hirayama 2005 conversion)\n\n")

# Create directories if needed
if (!dir.exists("data")) dir.create("data")
if (!dir.exists("results")) dir.create("results")


# Section 1: Sobp Model Predictions (Wenzl 2011A Table 3)

cat("\n--- SECTION 1: SOBP Model Predictions (Wenzl 2011A Table 3) ---\n")
cat("NOTE: These are MODEL-PREDICTED values, not experimental measurements.\n")
wenzl_sobp <- tribble(
  ~ion, ~position, ~LET, ~pO2_mmHg, ~OER_survival, ~source_model,
  
  # Carbon ions SOBP (Kohno et al. 2005 model)
  "C", "proximal", 40, 0.5, 1.88, "Kohno2005",
  "C", "proximal", 40, 1.0, 1.68, "Kohno2005",
  "C", "proximal", 40, 2.5, 1.39, "Kohno2005",
  "C", "proximal", 40, 5.0, 1.21, "Kohno2005",
  "C", "proximal", 40, 20, 1.03, "Kohno2005",
  "C", "distal", 223, 0.5, 1.27, "Kohno2005",
  "C", "distal", 223, 1.0, 1.23, "Kohno2005",
  "C", "distal", 223, 2.5, 1.15, "Kohno2005",
  "C", "distal", 223, 5.0, 1.09, "Kohno2005",
  "C", "distal", 223, 20, 1.01, "Kohno2005",
  "C", "mean", NA, 0.5, 1.76, "Kohno2005",
  "C", "mean", NA, 1.0, 1.59, "Kohno2005",
  "C", "mean", NA, 2.5, 1.34, "Kohno2005",
  "C", "mean", NA, 5.0, 1.19, "Kohno2005",
  "C", "mean", NA, 20, 1.02, "Kohno2005",
  
  # Proton SOBP (Wilkens & Oelfke 2004 model)
  "proton", "proximal", 2, 0.5, 2.11, "WilkensOelfke2004",
  "proton", "proximal", 2, 1.0, 1.83, "WilkensOelfke2004",
  "proton", "proximal", 2, 2.5, 1.46, "WilkensOelfke2004",
  "proton", "proximal", 2, 5.0, 1.24, "WilkensOelfke2004",
  "proton", "proximal", 2, 20, 1.03, "WilkensOelfke2004",
  "proton", "distal", 13, 0.5, 2.04, "WilkensOelfke2004",
  "proton", "distal", 13, 1.0, 1.78, "WilkensOelfke2004",
  "proton", "distal", 13, 2.5, 1.43, "WilkensOelfke2004",
  "proton", "distal", 13, 5.0, 1.23, "WilkensOelfke2004",
  "proton", "distal", 13, 20, 1.03, "WilkensOelfke2004",
  
  # Gamma reference (Alper-Howard-Flanders model)
  "photon", "reference", 0.22, 0.5, 2.13, "AlperHowardFlanders",
  "photon", "reference", 0.22, 1.0, 1.84, "AlperHowardFlanders",
  "photon", "reference", 0.22, 2.5, 1.46, "AlperHowardFlanders",
  "photon", "reference", 0.22, 5.0, 1.25, "AlperHowardFlanders",
  "photon", "reference", 0.22, 20, 1.03, "AlperHowardFlanders"
) %>%
  mutate(
    O2_pct = pO2_mmHg / 7.6,
    OER_retention = OER_survival * CONVERSION_FACTOR$mean,
    dataset = "Wenzl2011A_SOBP_model",
    data_type = "model_prediction"
  )

cat(sprintf("SOBP model data: %d points\n\n", nrow(wenzl_sobp)))


# Section 2: Hirayama 2005 Data

cat("\n--- SECTION 2: Hirayama et al. 2005 Data ---\n")
cat("Source: Hirayama et al. (2005) J. Radiat. Res.\n")
cat("Cell line: CHO (Chinese Hamster Ovary)\n")
cat("Key finding: OER_DSB at 15 min ≈ 1.2 × OER_killing\n\n")

hirayama_2005 <- tribble(
  ~ion, ~LET, ~LET_se, ~cell_line, ~OER_survival, ~OER_survival_se, ~OER_dsb, ~OER_dsb_se, ~source,
  
  # From Table 4 of Hirayama 2005
  "photon", 1.7, NA, "CHO", 2.8, 0.2, 3.4, 0.2, "Hirayama2005",
  "C", 79.6, 0.6, "CHO", 1.8, 0.0, 2.2, 0.1, "Hirayama2005"
) %>%
  mutate(
    OER_retention = OER_survival * CONVERSION_FACTOR$mean,
    OER_retention_se = sqrt((OER_survival_se * CONVERSION_FACTOR$mean)^2 + 
                              (OER_survival * CONVERSION_FACTOR$sd)^2),
    dataset = "Hirayama2005",
    data_type = "experimental",
    notes = "OER_killing from D10 ratio; OER_dsb from DSB measurements"
  )

cat("Hirayama 2005 OER values:\n")
cat("┌───────────┬────────┬───────────┬───────────────┬─────────────────┐\n")
cat("│ Radiation │ LET    │ Cell Line │ OER_survival  │ OER_retention   │\n")
cat("├───────────┼────────┼───────────┼───────────────┼─────────────────┤\n")
for (i in 1:nrow(hirayama_2005)) {
  cat(sprintf("│ %-9s │ %6.1f │ %-9s │ %6.2f ± %.2f │ %6.2f          │\n",
              hirayama_2005$ion[i],
              hirayama_2005$LET[i],
              hirayama_2005$cell_line[i],
              hirayama_2005$OER_survival[i],
              hirayama_2005$OER_survival_se[i],
              hirayama_2005$OER_retention[i]))
}
cat("└───────────┴────────┴───────────┴───────────────┴─────────────────┘\n\n")


# Section 3: Hirayama 2009 - Direct/Indirect Action Data

cat("\n--- SECTION 3: Hirayama et al. 2009 - Direct/Indirect Action ---\n")
cat("Source: Hirayama et al. (2009) - DMSO scavenger experiments\n")
cat("Cell line: V79\n")
cat("Key finding: Indirect action remains ~32% even at 2106 keV/μm\n\n")

hirayama_2009 <- tribble(
  ~radiation, ~LET, ~indirect_fraction, ~indirect_se, ~D10_inv, ~source,
  
  # From Table 1 of Hirayama 2009
  "X-rays", 9.4, 0.76, 0.05, 0.17, "Hirayama2009",
  "Carbon", 20, 0.65, 0.02, 0.14, "Hirayama2009",
  "Iron", 200, 0.50, 0.21, 0.48, "Hirayama2009",
  "Iron", 797, 0.52, 0.04, 0.25, "Hirayama2009",
  "Iron", 1298, 0.39, 0.04, 0.19, "Hirayama2009",
  "Iron", 2106, 0.32, 0.02, 0.15, "Hirayama2009"
) %>%
  mutate(
    direct_fraction = 1 - indirect_fraction,
    dataset = "Hirayama2009",
    cell_line = "V79",
    notes = "Indirect action from DMSO protection experiments"
  )

cat("Direct/Indirect Action Fractions (Hirayama 2009):\n")
cat("┌───────────┬────────────┬─────────────┬─────────────┐\n")
cat("│ Radiation │ LET        │ Direct (%)  │ Indirect (%)│\n")
cat("├───────────┼────────────┼─────────────┼─────────────┤\n")
for (i in 1:nrow(hirayama_2009)) {
  cat(sprintf("│ %-9s │ %10.0f │ %8.0f    │ %8.0f    │\n",
              hirayama_2009$radiation[i],
              hirayama_2009$LET[i],
              hirayama_2009$direct_fraction[i] * 100,
              hirayama_2009$indirect_fraction[i] * 100))
}
cat("└───────────┴────────────┴─────────────┴─────────────┘\n\n")

cat("This justifies p1_high = 0.68 (68% direct at extreme LET)\n\n")


# Section 4: Tinganelli Et Al. (2015) Data - New

cat("\n--- SECTION 4: Tinganelli et al. (2015) Data - NEW ---\n")
cat("Source: Tinganelli et al. (2015) Sci Rep. 5:17016\n")
cat("Cell line: CHO (Chinese Hamster Ovary)\n")
cat("Endpoint: 10% survival (D10)\n")
cat("Key features:\n")
# Tinganelli data with proper O2 levels
# Note: 0% in paper means anoxia, we use 0.001% for calculations
tinganelli_2015 <- tribble(
  ~ion, ~O2_pct_reported, ~LET, ~OER_survival, ~OER_survival_se,
  
  # Photon data (200 kVp X-rays, LET ~ 9.4 keV/μm)
  "photon", 0,    9.4,   2.83, 0.37,
  "photon", 0,    9.4,   2.40, 0.10,
  "photon", 0.5,  9.4,   1.62, 0.30,
  "photon", 0.5,  9.4,   1.50, 0.10,
  
  # Silicon data
  "Si", 2,    308,  1.07, 0.07,
  "Si", 0,    283,  1.35, 0.12,
  "Si", 0.5,  304,  1.08, 0.11,
  
  # Oxygen data
  "O",  0,    140,  1.44, 0.20,
  
  # Nitrogen data
  "N",  0,    160,  1.30, 0.04,
  
  # Carbon data - Physoxia (2% O2)
  "C",  2,    68.1, 1.20, 0.13,
  "C",  2,    51.2, 1.12, 0.09,
  "C",  2,    29.8, 1.18, 0.14,
  
  # Carbon data - Anoxia (0% O2)
  "C",  0,    317,  1.28, 0.09,
  "C",  0,    141,  1.49, 0.19,
  "C",  0,    51.1, 2.46, 0.23,
  "C",  0,    150,  1.33, 0.05,
  "C",  0,    100,  1.81, 0.12,
  
  # Carbon data - Mild hypoxia (0.5% O2)
  "C",  0.5,  124,  1.16, 0.09,
  "C",  0.5,  51.5, 1.40, 0.15,
  "C",  0.5,  100,  1.29, 0.08,
  
  # Carbon data - Acute hypoxia (0.15% O2)
  "C",  0.15, 124,  1.34, 0.12,
  "C",  0.15, 68,   1.64, 0.18,
  "C",  0.15, 29.8, 1.90, 0.16
) %>%
  mutate(
    # Convert 0% to 0.001% for calculations
    O2_pct = ifelse(O2_pct_reported == 0, 0.001, O2_pct_reported),
    O2_hyp = O2_pct,
    O2_ref = 21.0,
    
    # Convert to retention OER
    OER_retention = OER_survival * CONVERSION_FACTOR$mean,
    OER_retention_se = sqrt((OER_survival_se * CONVERSION_FACTOR$mean)^2 + 
                              (OER_survival * CONVERSION_FACTOR$sd)^2),
    
    # Weights based on uncertainty
    weight = 1 / OER_retention_se^2,
    
    # Metadata
    cell_line = "CHO",
    cell_line_std = "CHO",
    source = "Tinganelli2015",
    dataset = "Tinganelli2015",
    data_type = "experimental",
    
    # Condition labels
    condition = case_when(
      O2_pct_reported == 0 ~ "Anoxia",
      O2_pct_reported == 0.15 ~ "Acute hypoxia",
      O2_pct_reported == 0.5 ~ "Mild hypoxia",
      O2_pct_reported == 2 ~ "Physoxia",
      TRUE ~ "Other"
    )
  )

cat("Tinganelli et al. (2015) data summary:\n")
cat(sprintf("  Total observations: %d\n", nrow(tinganelli_2015)))
cat(sprintf("  Particles: %s\n", paste(unique(tinganelli_2015$ion), collapse = ", ")))
cat(sprintf("  O2 levels: %s%%\n", paste(sort(unique(tinganelli_2015$O2_pct_reported)), collapse = ", ")))
cat(sprintf("  LET range: %.1f - %.1f keV/μm\n", 
            min(tinganelli_2015$LET), max(tinganelli_2015$LET)))
cat("By condition:\n")
tinganelli_2015 %>%
  group_by(condition, O2_pct_reported) %>%
  summarise(
    n = n(),
    particles = paste(unique(ion), collapse = ", "),
    LET_range = sprintf("%.0f-%.0f", min(LET), max(LET)),
    OER_range = sprintf("%.2f-%.2f", min(OER_survival), max(OER_survival)),
    .groups = "drop"
  ) %>%
  print()
cat("By particle:\n")
tinganelli_2015 %>%
  group_by(ion) %>%
  summarise(
    n = n(),
    LET_range = sprintf("%.0f-%.0f", min(LET), max(LET)),
    O2_levels = paste(sort(unique(O2_pct_reported)), collapse = ", "),
    OER_mean = sprintf("%.2f", mean(OER_survival)),
    .groups = "drop"
  ) %>%
  print()
# Section 5: Comprehensive Literature Oer Data

cat("\n--- SECTION 5: Comprehensive Literature OER Data ---\n")
cat("All values are OER at 10% survival from original publications.\n\n")

# Comprehensive literature compilation (same as before)
literature_oer <- tribble(
  ~source, ~ion, ~cell_line, ~LET, ~OER_survival,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Prise et al. 1990
  # ═══════════════════════════════════════════════════════════════════════════
  "Prise1990", "photon", "V79", 2, 3.14,
  "Prise1990", "He", "V79", 110, 0.99,
  "Prise1990", "proton", "V79", 17, 2.77,
  "Prise1990", "proton", "V79", 24, 2.02,
  "Prise1990", "proton", "V79", 32, 1.89,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Chapman 1977
  # ═══════════════════════════════════════════════════════════════════════════
  "Chapman1977", "photon", "V79", 2, 2.8,
  "Chapman1977", "He", "V79", 8, 2.3,
  "Chapman1977", "C", "V79", 80, 1.7,
  "Chapman1977", "Ne", "V79", 120, 1.6,
  "Chapman1977", "Ar", "V79", 310, 1.4,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Barendsen 1966
  # ═══════════════��═══════════════════════════════════════════════════════════
  "Barendsen1966", "photon", "T1", 1.3, 2.7,
  "Barendsen1966", "He", "T1", 166, 1.0,
  "Barendsen1966", "He", "T1", 140, 1.1,
  "Barendsen1966", "He", "T1", 110, 1.3,
  "Barendsen1966", "He", "T1", 88, 1.7,
  "Barendsen1966", "He", "T1", 61, 2.05,
  "Barendsen1966", "He", "T1", 26, 2.4,
  "Barendsen1966", "deuteron", "T1", 20, 2.4,
  "Barendsen1966", "deuteron", "T1", 5.3, 2.6,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Staab et al. 2004
  # ═══════════════════════════════════════════════════════════════════════════
  "Staab2004", "photon", "V79", 1.7, 2.87,
  "Staab2004", "C", "V79", 18, 2.76,
  "Staab2004", "C", "V79", 60, 1.40,
  
  # ════════════════════════════��══════════════════════════════════════════════
  # Tenforde et al. 1980
  # ═══════════════════════════════════════════════════════════════════════════
  "Tenforde1980", "photon", "R1", 2, 2.2,
  "Tenforde1980", "C", "R1", 95, 1.9,
  "Tenforde1980", "Ne", "R1", 177, 1.7,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Katz and Sharma 1974
  # ═══════════════════════════════════════════════════════════════════════════
  "Katz1974", "proton", "T1", 5, 2.6,
  "Katz1974", "proton", "p388", 5, 2.3,
  "Katz1974", "proton", "HeLa", 5, 2.8,
  "Katz1974", "He", "T1", 8, 2.5,
  "Katz1974", "He", "p388", 8, 2.3,
  "Katz1974", "He", "HeLa", 8, 2.6,
  "Katz1974", "Ne", "T1", 46, 1.9,
  "Katz1974", "Ne", "p388", 46, 1.2,
  "Katz1974", "Ne", "HeLa", 46, 1.7,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Raju et al. 1978
  # ═══════════════════════════════════════════════════════════════════════════
  "Raju1978", "photon", "V79", 2, 2.9,
  "Raju1978", "proton", "V79", 0.7, 3.1,
  "Raju1978", "proton", "V79", 1.9, 2.8,
  "Raju1978", "He", "V79", 1.8, 3.0,
  "Raju1978", "He", "V79", 3.575, 2.5,
  "Raju1978", "He", "V79", 7.1, 2.6,
  "Raju1978", "C", "V79", 11.1, 3.4,
  "Raju1978", "C", "V79", 20.073, 2.5,
  "Raju1978", "C", "V79", 36.3, 2.3,
  "Raju1978", "Ne", "V79", 31.7, 2.3,
  "Raju1978", "Ne", "V79", 51.541, 2.0,
  "Raju1978", "Ne", "V79", 83.8, 2.0,
  "Raju1978", "Ar", "V79", 93.2, 1.6,
  "Raju1978", "Ar", "V79", 147.362, 1.9,
  "Raju1978", "Ar", "V79", 233, 1.5,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Raju et al. 1972
  # ═══════════════════════════════════════════════════════════════════════════
  "Raju1972", "He", "T1", 3, 2.5,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Raju et al. 1979
  # ═══════════════════════════════════════════════════════════════════════════
  "Raju1979", "photon", "V79", 2, 2.8,
  "Raju1979", "photon", "V79", 2, 2.9,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Raju et al. 1987
  # ═══════════════════════════════════════════════════════════════════════════
  "Raju1987", "photon", "V79", 20, 2.0,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Fu and Phillips 1976
  # ═══════════════════════════════════════════════════════════════════════════
  "FuPhillips1976", "photon", "EMT6", 1, 3.23,
  "FuPhillips1976", "Ne", "EMT6", 31, 2.28,
  "FuPhillips1976", "Ne", "EMT6", 180, 1.35,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Curtis 1976
  # ═══════════════════════════════════════════════════════════════════════════
  "Curtis1976", "He", "R1", 110, 1.15,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Skarsgard and Harrison 1991
  # ═══════════════════════════════════════════════════════════════════════════
  "Skarsgard1991", "photon", "V79", 2, 2.76,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Hill et al. 2002
  # ═══════════════════════════════════════════════════════════════════════════
  "Hill2002", "photon", "V79", 20, 2.38,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Blakely et al. 1979 - Carbon
  # ═══════════════════════════════════════════════════════════════════════════
  "Blakely1979", "C", "T1", 11, 3.09,
  "Blakely1979", "C", "T1", 13, 2.91,
  "Blakely1979", "C", "T1", 16, 2.87,
  "Blakely1979", "C", "T1", 23, 2.68,
  "Blakely1979", "C", "T1", 30, 2.6,
  "Blakely1979", "C", "T1", 40, 2.33,
  "Blakely1979", "C", "T1", 85, 2.2,
  "Blakely1979", "C", "T1", 124, 2.0,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Blakely et al. 1979 - Neon
  # ═══════════════════════════════════════════════════════════════════════════
  "Blakely1979", "Ne", "T1", 32, 3.03,
  "Blakely1979", "Ne", "T1", 38, 2.92,
  "Blakely1979", "Ne", "T1", 54, 2.31,
  "Blakely1979", "Ne", "T1", 71, 2.25,
  "Blakely1979", "Ne", "T1", 100, 2.05,
  "Blakely1979", "Ne", "T1", 139, 1.95,
  "Blakely1979", "Ne", "T1", 234, 1.5,
  "Blakely1979", "Ne", "T1", 419, 1.36,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Blakely et al. 1979 - Argon
  # ═══════════════════════════════════════════════════════════════════════════
  "Blakely1979", "Ar", "T1", 81, 2.6,
  "Blakely1979", "Ar", "T1", 91, 2.15,
  "Blakely1979", "Ar", "T1", 117, 1.95,
  "Blakely1979", "Ar", "T1", 144, 1.6,
  "Blakely1979", "Ar", "T1", 184, 1.46,
  "Blakely1979", "Ar", "T1", 245, 1.4,
  "Blakely1979", "Ar", "T1", 328, 1.3,
  "Blakely1979", "Ar", "T1", 640, 1.28,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Curtis et al. 1982
  # ═══════════════════════════════════════════════════════════════════════════
  "Curtis1982", "C", "R1", 11, 2.34,
  "Curtis1982", "C", "R1", 64, 1.65,
  "Curtis1982", "C", "R1", 95, 1.72,
  "Curtis1982", "Ne", "R1", 31, 2.20,
  "Curtis1982", "Ne", "R1", 120, 1.72,
  "Curtis1982", "Ne", "R1", 177.5, 1.55,
  "Curtis1982", "Ar", "R1", 95, 1.82,
  "Curtis1982", "Ar", "R1", 430, 1.26,
  "Curtis1982", "Ar", "R1", 750, 1.20,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Sprong et al. 2006
  # ═══════════════════════════════════════════════════════════════════════════
  "Sprong2006", "photon", "AA8", 2, 2.7,
  "Sprong2006", "photon", "CHO-9", 2, 2.7,
  "Sprong2006", "photon", "V79", 2, 3.0,
  "Sprong2006", "photon", "V79-B", 2, 2.6,
  "Sprong2006", "photon", "BJhtert", 2, 2.0,
  "Sprong2006", "photon", "HF30", 2, 2.1,
  "Sprong2006", "photon", "HF293", 2, 2.4,
  
  # ═══════════════════════════════════════════════════════════════════��═══════
  # Todd et al. 1974
  # ═══════════════════════════════════════════════════════════════════════════
  "Todd1974", "He", "V79", 3, 3.0,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Freyer et al. 1991
  # ═══════════════════════════════════════════════════════════════════════════
  "Freyer1991", "photon", "CHO-xrs6", 0.22, 2.89,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Hirayama et al. 2005
  # ═══════════════════════════════════════════════════════════════════════════
  "Hirayama2005", "photon", "CHO", 1.7, 2.8,
  "Hirayama2005", "C", "CHO", 79.6, 1.8,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Hall et al. 1966
  # ════════════════════════════════════════════════���══════════════════════════
  "Hall1966", "photon", "HeLa", 0.22, 2.64,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Nias et al. 1967
  # ═══════════════════════════════════════════════════════════════════════════
  "Nias1967", "photon", "HeLa", 2, 2.39,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Guichard et al. 1977
  # ═══════════════════════════════════════════════════════════════════════════
  "Guichard1977", "photon", "EMT6", 0.22, 2.7,
  "Guichard1977", "He", "EMT6", 10, 2.9,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Wheeler et al. 1979
  # ════════════════════════════════════════════════════════════════���══════════
  "Wheeler1979", "photon", "9L", 2, 2.1,
  "Wheeler1979", "C", "9L", 11, 2.0,
  "Wheeler1979", "C", "9L", 64, 1.7,
  
  # ════════════��══════════════════════════════════════════════════════════════
  # Hall et al. 1977
  # ═══════════════════════════════════════════════════════════════════════════
  "Hall1977", "Ar", "V79", 110.9, 2.1,
  "Hall1977", "Ar", "V79", 304.9, 1.45,
  "Hall1977", "Ar", "V79", 409.2, 1.32,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Ando et al. 1999
  # ═══════════════════════════════════════════════════════════════════════════
  "Ando1999", "photon", "NFSa", 1, 3.4,
  "Ando1999", "C", "NFSa", 74, 1.6,
  
  # ═══════════════════════════════════════════════════════════════════════════
  # Feola et al. 1969
  # ═══════════════════════════════════════════════════════════════════════════
  "Feola1969", "photon", "L2", 2, 3.17,
  "Feola1969", "He", "L2", 1.7, 3.30,
  "Feola1969", "He", "L2", 18, 2.10,
  "Feola1969", "He", "L2", 8, 3.48,
  "Feola1969", "He", "L2", 22, 1.82
)

# Process literature data
literature_oer <- literature_oer %>%
  mutate(
    # Default O2 levels for anoxic OER measurements
    O2_hyp = 0.001,
    O2_ref = 21.0,
    
    # Convert to retention OER
    OER_retention = OER_survival * CONVERSION_FACTOR$mean,
    
    # Assign uncertainty (10% for well-documented studies, 15% for others)
    OER_survival_se = OER_survival * 0.12,
    OER_retention_se = OER_survival * 0.12 * CONVERSION_FACTOR$mean,
    
    # Weights
    weight = 1 / OER_retention_se^2,
    
    # Standardize cell line names
    cell_line_std = case_when(
      is.na(cell_line) ~ "Unknown",
      str_detect(toupper(cell_line), "V79") ~ "V79",
      str_detect(toupper(cell_line), "HSG") ~ "HSG",
      str_detect(toupper(cell_line), "CHO") ~ "CHO",
      str_detect(toupper(cell_line), "T1") ~ "T1",
      str_detect(toupper(cell_line), "R1") ~ "Other",
      str_detect(toupper(cell_line), "HELA") ~ "Other",
      str_detect(toupper(cell_line), "EMT6") ~ "Other",
      str_detect(toupper(cell_line), "9L") ~ "Other",
      str_detect(toupper(cell_line), "L2") ~ "Other",
      str_detect(toupper(cell_line), "P388") ~ "Other",
      str_detect(toupper(cell_line), "NFSA") ~ "Other",
      TRUE ~ "Other"
    ),
    
    # Standardize ion names
    ion_std = ion,
    
    # Dataset and data type
    dataset = source,
    data_type = "experimental",
    
    # Log LET for modeling
    LET_log = log10(LET)
  )


# Section 6: Combine All Calibration Data

cat("\n═══════════════════════════════════════════════════════════════════════\n")
# Prepare Tinganelli data for combination
tinganelli_for_combine <- tinganelli_2015 %>%
  select(
    ion, LET, O2_hyp, O2_ref, 
    OER_survival, OER_survival_se, OER_retention, OER_retention_se,
    weight, cell_line, cell_line_std, source, dataset, data_type
  ) %>%
  mutate(ion_std = ion)

# Prepare literature data for combination
literature_for_combine <- literature_oer %>%
  select(
    ion, LET, O2_hyp, O2_ref,
    OER_survival, OER_survival_se, OER_retention, OER_retention_se,
    weight, cell_line, cell_line_std, source, dataset, data_type, ion_std
  )

# Combine
all_calibration_data <- bind_rows(
  literature_for_combine,
  tinganelli_for_combine
)

# Summary statistics
cat(sprintf("Total calibration data points: %d\n", nrow(all_calibration_data)))
cat(sprintf("  Literature (various sources): %d\n", nrow(literature_for_combine)))
cat(sprintf("  Tinganelli 2015 (CHO): %d\n\n", nrow(tinganelli_for_combine)))

cat("By particle type:\n")
print(all_calibration_data %>%
        group_by(ion) %>%
        summarise(
          n = n(),
          LET_min = round(min(LET), 1),
          LET_max = round(max(LET), 1),
          OER_min = round(min(OER_survival), 2),
          OER_max = round(max(OER_survival), 2),
          OER_mean = round(mean(OER_survival), 2),
          .groups = "drop"
        ) %>%
        arrange(match(ion, c("photon", "proton", "deuteron", "He", "C", "N", "O", "Ne", "Si", "Ar"))))

cat("\nBy cell line:\n")
print(all_calibration_data %>%
        count(cell_line_std) %>%
        arrange(desc(n)))

cat("\nBy source:\n")
print(all_calibration_data %>%
        count(source) %>%
        arrange(desc(n)))


# Section 7: Data Quality Checks

cat("\n═══════════════════════════════════════════════════════════════════════\n")
# Check for potential outliers
cat("Checking for potential outliers...\n\n")

outlier_candidates <- all_calibration_data %>%
  mutate(
    potential_outlier = case_when(
      # OER < 1 is physically impossible
      OER_survival < 1.0 ~ "OER < 1.0 (impossible)",
      
      # Very high OER for high LET
      LET > 100 & OER_survival > 2.5 ~ "High LET with OER > 2.5",
      
      # Very low OER for low LET photons
      ion == "photon" & LET < 5 & OER_survival < 2.0 ~ "Low-LET photon with OER < 2.0",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(potential_outlier))

if (nrow(outlier_candidates) > 0) {
  cat("⚠ Potential outliers identified:\n")
  print(outlier_candidates %>% 
          select(source, ion, cell_line, LET, OER_survival, potential_outlier))
} else {
  cat("✓ No obvious outliers detected\n")
}

# Check OER vs LET trend
cat("\nOER vs LET trend (should decrease with increasing LET):\n")
let_trend <- all_calibration_data %>%
  mutate(LET_bin = cut(LET, 
                       breaks = c(0, 10, 30, 60, 100, 200, 400, 1000),
                       labels = c("<10", "10-30", "30-60", "60-100", "100-200", "200-400", ">400"))) %>%
  group_by(LET_bin) %>%
  summarise(
    n = n(),
    OER_mean = round(mean(OER_survival), 2),
    OER_sd = round(sd(OER_survival), 2),
    .groups = "drop"
  )

print(let_trend)


# Section 8: D-Kondo Oxygen Curve (Validation Only)

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("Source: D-Kondo (2005/2009)\n")
cat("NOTE: Reports HRF (Hypoxia Reduction Factor) with radioprotector\n")
dkondo_oxygen <- tribble(
  ~O2_pct, ~HRF, ~is_measured,
  21.0, 1.0, TRUE,
  10.0, 1.05, FALSE,
  5.0, 1.10, FALSE,
  2.1, 1.21, TRUE,
  1.0, 1.45, FALSE,
  0.5, 1.75, FALSE,
  0.21, 2.33, TRUE,
  0.1, 2.55, FALSE,
  0.05, 2.70, FALSE,
  0.021, 2.83, TRUE,
  0.01, 2.88, FALSE,
  0.001, 2.95, FALSE
) %>%
  mutate(
    ion = "photon",
    LET = 1.0,
    O2_ref = 21.0,
    dataset = "D-Kondo",
    data_type = "validation_HRF",
    notes = ifelse(is_measured, "Measured HRF", "Interpolated")
  )

cat("D-Kondo Oxygen Curve (HRF values):\n")
print(dkondo_oxygen %>% select(O2_pct, HRF, is_measured))


# Section 9: Combine And Save All Data

cat("\n═══════════════════════════════════════════════════════════════��═══════\n")
# Save calibration data (literature + Tinganelli OER)
write_csv(all_calibration_data, "data/wenzl_literature_oer.csv")

# Save Tinganelli data separately for reference
write_csv(tinganelli_2015, "data/tinganelli_2015_oer.csv")

# Save SOBP model data (validation only)
write_csv(wenzl_sobp, "data/wenzl_sobp_oxygen.csv")

# Save D-Kondo data (validation only)
write_csv(dkondo_oxygen, "data/dkondo_oxygen_curve.csv")

# Save Hirayama 2009 direct/indirect data
write_csv(hirayama_2009, "data/hirayama_2009_direct_indirect.csv")

# Save complete documentation
complete_data <- list(
  # Calibration data
  all_calibration = all_calibration_data,
  literature_oer = literature_oer,
  tinganelli_2015 = tinganelli_2015,
  
  # Validation data
  sobp_model = wenzl_sobp,
  dkondo_oxygen = dkondo_oxygen,
  hirayama_2005 = hirayama_2005,
  hirayama_2009 = hirayama_2009,
  
  # Metadata
  metadata = list(
    created = Sys.Date(),
    n_total_points = nrow(all_calibration_data),
    n_literature_points = nrow(literature_oer),
    n_tinganelli_points = nrow(tinganelli_2015),
    n_sources = length(unique(all_calibration_data$source)),
    particles = unique(all_calibration_data$ion),
    cell_lines = unique(all_calibration_data$cell_line),
    
    notes = c(
      "All OER values at 10% survival (D10-based)",
      "OER_retention = OER_survival × 1.20 (Hirayama 2005)",
      "Tinganelli 2015: CHO cells, multiple O2 levels",
      "SOBP data is MODEL-PREDICTED, not experimental",
      "D-Kondo reports HRF with radioprotector, not direct OER",
      "New particles from Tinganelli: N, O, Si"
    )
  )
)

saveRDS(complete_data, "data/step2_complete_data.rds")

cat("Saved:\n")
# Section 10: Summary

cat(sprintf("║   Total calibration data: %d points from %d sources               ║\n", 
            nrow(all_calibration_data), length(unique(all_calibration_data$source))))
cat(sprintf("║     Literature (various): %d points                               ║\n", 
            nrow(literature_oer)))
cat(sprintf("║     Tinganelli 2015 (CHO): %d points                               ║\n", 
            nrow(tinganelli_2015)))
for (p in c("photon", "proton", "deuteron", "He", "C", "N", "O", "Ne", "Si", "Ar")) {
  n_p <- sum(all_calibration_data$ion == p)
  if (n_p > 0) {
    cat(sprintf("║     %-10s: %3d points                                        ║\n", p, n_p))
  }
}
for (cl in c("V79", "HSG", "CHO", "T1", "Other")) {
  n_cl <- sum(all_calibration_data$cell_line_std == cl)
  if (n_cl > 0) {
    cat(sprintf("║     %-10s: %3d points                                        ║\n", cl, n_cl))
  }
}
cat(sprintf("║     SOBP model predictions: %d points                             ║\n", nrow(wenzl_sobp)))
cat(sprintf("║     D-Kondo oxygen curve: %d points                               ║\n", nrow(dkondo_oxygen)))
cat(sprintf("║     Hirayama 2009 (direct/indirect): %d points                    ║\n", nrow(hirayama_2009)))
cat("╚═════════════════════════════════════════════════���════════════════════╝\n")