# ================================================================================
# Correlated Response to Selection on Methane Emission Traits
# ================================================================================
#
# Description:
#   Computes the correlated response to selection for 27 traits and
#   composites under each of the 12 U.S. lifetime net merit selection indices
#   from 1971 (PD$) to 2025 (NM$).
#   Special focus is placed on enteric methane emission traits (MEP, RMI, RMY).
#
#   The correlated response formula used is:
#       CR = i * (b'G) / sqrt(b'Pb) / L
#   where:
#       i  = selection intensity (0.35, consistent with 2025 NM$ derivation)
#       b  = vector of economic values for traits in the index ($/PTA unit, 2017$)
#       G  = matrix of genetic covariances between index traits and all 27 traits
#       P  = phenotypic covariance matrix of index traits, weighted by PTA reliabilities
#       L  = generation interval (2.2 yr for Holstein bulls; Guinan et al., 2023)
#
#   Economic values are inflation-adjusted to 2017 U.S. dollars using CPI-U-X1
#   correction factors prior to computing each index's correlated response.
#
# ================================================================================

# ---- Libraries ----------------------------------------------------------------
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(patchwork)
library(forcats)

rm(list = ls())

# ---- Working directory --------------------------------------------------------
setwd("/Users/GuillermoMartinez/Documents/Projects/Project_NMICH4")

# ================================================================================
# SECTION 1: PARAMETERS AND MATRIX DEFINITIONS
# ================================================================================

# --- 1.1 Genetic (co)variance matrix -------------------------------------------
# Genetic correlation matrix estimated with Calo's method (Calo et al., 1973)
# Sheet "Calo 2026_05_18" contains the final adjusted matrix used in analyses.
genCORR <- as.matrix(readxl::read_excel("Aim1_Past/Data/genCORR.xlsx", sheet = "Calo 2026_05_18"))

# Genetic SDs (TTA basis) and heritabilities for all 27 traits
val <- readxl::read_excel("Aim1_Past/Data/genSD_h2.xlsx", col_types = c("text", "numeric", "numeric"))

# Build genetic covariance matrix: G = D * R * D, where D = diag(genSD)
genCOV <- genCORR * (val$genSD %*% t(val$genSD))
print(genCOV[, 25:27])   # Preview methane trait columns (MEP, RMI, RMY)

# --- 1.2 Phenotypic (co)variance matrix ----------------------------------------
# Phenotypic correlations computed assuming zero environmental correlations:
#   rP(x,y) = h_x * h_y * rG(x,y),  where h = sqrt(h2)
# See Methods for justification of this assumption.
phenCORR <- as.matrix(readxl::read_excel("Aim1_Past/Data/phenCORR.xlsx", sheet = "Calc_rP"))

# Phenotypic SDs derived from genetic SDs and heritabilities: SD(P) = SD(G) / sqrt(h2)
val$phenSD <- val$genSD / sqrt(val$h2)

# Build phenotypic covariance matrix
phenCOV <- phenCORR * (val$phenSD %*% t(val$phenSD))
print(phenCOV[, 25:27])   # Preview methane trait columns

# --- 1.3 Index matrices (P, G, C) ----------------------------------------------
# The selection index formulation uses three matrices:
#   P (17x17): phenotypic (co)variances of traits included in the full 2025 index
#   G (17x27): genetic covariances between index traits and all 27 criterion traits
#   C (27x27): genetic (co)variances among all criterion traits (= genCOV)
#
# The 17 index traits (in order) are the traits present in NM$ 2025:
#   MY, FAT, PROT, PL, SCS, BWC, UDC, FLC, DPR, CA$, HCR, CCR, LIV [1:13]
#   RFI, HTH$ [rows 15-16 of genCOV → index rows 14-15]
#   EFC, HLIV [rows 23-24 of genCOV → index rows 16-17]
#
# Note: Milk (index row 1) maps to genCOV row 1; HTH$ maps to genCOV row 15, etc.

# P matrix: weighted phenotypic covariance matrix
# Off-diagonal elements are scaled by the geometric mean of PTA reliabilities
# for each trait pair, following Cole (2025). Composite REL approximations:
#   CA$  → REL(SCE);  HTH$ → REL(MAS);  BWC, UDC, FLC → REL(Type)
Pwgt <- as.matrix(readxl::read_excel("Aim1_Past/Data/Cole 2025 Supplemental Table 2 Weighted Phenotypic Covariance Matrix.xlsx"))
P <- Pwgt

# G matrix (17x27): genetic covariances between the 17 index traits and all 27 traits
G <- matrix(0, nrow = 17, ncol = 27)
G[1:13, ] <- genCOV[1:13, ]
G[14:15, ] <- genCOV[15:16, ]
G[16:17, ] <- genCOV[23:24, ]

# C matrix (27x27): full genetic covariance matrix (criterion traits)
C <- genCOV

# Gstar (17x17): genetic covariances among the 17 index traits only
# Used to compute the heritability of each index
Gstar <- matrix(0, nrow = 17, ncol = 17)
Gstar[1:13, 1:13] <- genCOV[1:13, 1:13]
Gstar[1:13, 14:15] <- genCOV[1:13, 15:16]
Gstar[14:15, 1:13] <- genCOV[15:16, 1:13]
Gstar[14, 14] <- genCOV[15, 15]
Gstar[15, 15] <- genCOV[16, 16]
Gstar[14, 15] <- genCOV[15, 16]
Gstar[15, 14] <- genCOV[16, 15]
Gstar[1:13, 16:17] <- genCOV[1:13, 23:24]
Gstar[16:17, 1:13] <- genCOV[23:24, 1:13]
Gstar[14:15, 16:17] <- genCOV[15:16, 23:24]
Gstar[16:17, 14:15] <- genCOV[23:24, 15:16]
Gstar[16:17, 16:17] <- genCOV[23:24, 23:24]

# Trait identifiers (ordered as in genCOV/val)
tidx <- val$tidx

# Clean up intermediate matrices no longer needed
rm(GCORRstar, genCORR, genCOV)
rm(phenCORR, phenCOV, Pwgt)

# ================================================================================
# SECTION 2: CORRELATED RESPONSE CALCULATIONS
# ================================================================================
# For each of the 12 U.S. selection indices, the correlated response is computed as:
#       CR_yr = i * (b'G) / sqrt(b'Pb) / L
#
# Index heritability is computed as:
#       h2_index = (b'G*b) / (b'Pb)
#
# Economic values (b) are expressed in nominal dollars for each index revision year
# and inflation-adjusted to 2017 U.S. dollars using CPI-U-X1 correction factors
# before entry into the formula.
#
# Indices with fewer traits use the appropriate sub-matrices of G, P, and Gstar
# (e.g., G[1:3, ] for indices that only include milk, fat, and protein).

# Selection intensity and generation interval
i  <- 0.35   # Annual selection differential in SD units (VanRaden et al., 2025)
GI <- 2.2    # Generation interval in years (Guinan et al., 2023)

# --- 2.1 NM$ 2025 --------------------------------------------------------------
# Source: VanRaden et al. (2025). Net merit as a measure of lifetime profit:
#         2025 revision. USDA ARS Research Report NM$9 (01-25).
# Formula (nominal 2025$):
#   NM$2025 = 0.022(Milk) + 5.01(Fat) + 3.33(Prot) + 30(PL) - 74(SCS-3)
#             - 57(BWC) + 8(UDC) + 3(FLC) + 6(DPR) + 1(CA$) + 1.5(HCR)
#             + 4.3(CCR) + 14.3(LIV) + 1(HTH$) - 0.35(RFI) + 2(EFC) + 8.2(HLIV)
b_2025 <- matrix(c(0.022, 5.01, 3.33, 30.00, -74.00, -57.00, 8.00, 3.00, 6.00, 1.00, 1.50, 4.30, 14.30, 1.00, -0.35, 2.00, 8.20))
factor <- 1.199
b_2025 <- b_2025 / factor  # Adjust to 2017$

CR_2025 <- (t(b_2025) %*% G) / as.numeric(sqrt(t(b_2025) %*% P %*% b_2025))
CR_2025 <- (CR_2025*i)/GI
print(t(CR_2025))

responses <- data.frame(
  index = "NM$ 2025",
  year = 2025,
  trait = tidx,
  response = as.numeric(CR_2025)
)

# --- 2.2 NM$ 2021 --------------------------------------------------------------
# Source: VanRaden et al. (2021). Net merit as a measure of lifetime profit:
#         2021 revision. USDA ARS Research Report NM$8 (05-21).
# Key change: First inclusion of RFI (b = -0.30) and new traits EFC and HLIV.
# Formula (nominal 2021$):
#   NM$2021 = 0.002(Milk) + 4.18(Fat) + 4.67(Prot) + 34(PL) - 74(SCS-3)
#             - 45(BWC) + 19(UDC) + 3(FLC) + 11(DPR) + 1(CA$) + 1.1(HCR)
#             + 2.2(CCR) + 9.8(LIV) + 1(HTH$) - 0.30(RFI) + 2.1(EFC) + 5.0(HLIV)
b_2021 <- matrix(c(0.002, 4.18, 4.67, 34.00, -74.00, -45.00, 19.00, 3.00, 11.00, 1.00, 1.10, 2.20, 9.80, 1.00, -0.30, 2.10, 5.00))
factor <- 1.092
b_2021 <- b_2021 / factor

CR_2021 <- (t(b_2021) %*% G) / as.numeric(sqrt(t(b_2021) %*% P %*% b_2021))
CR_2021 <- (CR_2021*i)/GI
print(t(CR_2021))

resp2021 <- data.frame(
  index = "NM$ 2021",
  year = 2021,
  trait = tidx,
  response = as.numeric(CR_2021)
)
responses <- rbind(responses, resp2021)

# --- 2.3 NM$ 2018 --------------------------------------------------------------
# Source: VanRaden et al. (2018). Net merit as a measure of lifetime profit:
#         2018 revision. USDA AIPL. https://aipl.arsusda.gov/reference/nmcalc-2018.htm
# Key change: First inclusion of six health traits via HTH$ composite.
# Note: RFI not yet in index (b_RFI = 0). Index has 14 traits.
# Formula (nominal 2018$):
#   NM$2018 = -0.004(Milk) + 4.03(Fat) + 3.53(Prot) + 19(PL) - 72(SCS-3)
#             - 18(BWC) + 31(UDC) + 10(FLC) + 11(DPR) + 1(CA$) + 2.2(HCR)
#             + 2.2(CCR) + 12(LIV) + 1(HTH$)
b_2018 <- matrix(c(-0.004, 4.03, 3.53, 19.00, -72.00, -18.00, 31.00, 10.00, 11.00, 1.00, 2.20, 2.20, 12.00, 1.00))
factor <- 1.020
b_2018 <- b_2018 / factor

CR_2018 <- (t(b_2018) %*% G[1:14, ]) / as.numeric(sqrt(t(b_2018) %*% P[1:14, 1:14] %*% b_2018))
CR_2018 <- (CR_2018*i)/GI
print(t(CR_2018))

resp2018 <- data.frame(
  index = "NM$ 2018",
  year = 2018,
  trait = tidx,
  response = as.numeric(CR_2018)
)
responses <- rbind(responses, resp2018)

# --- 2.4 NM$ 2017 --------------------------------------------------------------
# Source: VanRaden (2017). Net merit as a measure of lifetime profit: 2017 revision.
#         USDA ARS Research Report NM$6 (02-17).
# Key change: Inclusion of LIV (cow livability). 13 traits.
# Formula (nominal 2017$):
#   NM$2017 = -0.004(Milk) + 3.56(Fat) + 3.81(Prot) + 21(PL) - 117(SCS-3)
#             - 20(BWC) + 31(UDC) + 10(FLC) + 11(DPR) + 1(CA$) + 2.2(HCR)
#             + 2.2(CCR) + 12(LIV)
b_2017 <- matrix(c(-0.004, 3.56, 3.81, 21.00, -117.00, -20.00, 31.00, 10.00, 11.00, 1.00, 2.20, 2.20, 12.00))
factor <- 1
b_2017 <- b_2017 / factor

CR_2017 <- (t(b_2017) %*% G[1:13,]) / as.numeric(sqrt(t(b_2017) %*% P[1:13, 1:13] %*% b_2017))
CR_2017 <- (CR_2017*i)/GI
print(t(CR_2017))

resp2017 <- data.frame(
  index = "NM$ 2017",
  year = 2017,
  trait = tidx,
  response = as.numeric(CR_2017)
)
responses <- rbind(responses, resp2017)

# --- 2.5 NM$ 2014 --------------------------------------------------------------
# Source: VanRaden and Cole (2014). Net merit as a measure of lifetime profit:
#         2014 revision. USDA ARS Research Report NM$5 (10-14).
# Key change: Addition of HCR and CCR (heifer and cow conception rates). 12 traits.
# Formula (nominal 2014$):
#   NM$2014 = -0.006(Milk) + 3.22(Fat) + 4.14(Prot) + 29(PL) - 122(SCS-3)
#             - 16(BWC) + 31(UDC) + 10(FLC) + 11(DPR) + 1(CA$) + 2.3(HCR)
#             + 2.2(CCR)
b_2014 <- matrix(c(-0.006, 3.22, 4.14, 29.00, -122.00, -16.00, 31.00, 10.00, 11.00, 1.00, 2.30, 2.20))
factor <- 0.966
b_2014 <- b_2014 / factor  # Adjust to 2017$

CR_2014 <- (t(b_2014) %*% G[1:12, ]) / as.numeric(sqrt(t(b_2014) %*% P[1:12, 1:12] %*% b_2014))
CR_2014 <- (CR_2014*i)/GI
print(t(CR_2014))

resp2014 <- data.frame(
  index = "NM$ 2014",
  year = 2014,
  trait = tidx,
  response = as.numeric(CR_2014)
)
responses <- rbind(responses, resp2014)

# --- 2.6 NM$ 2010 --------------------------------------------------------------
# Source: Cole, VanRaden, and Multi-State Project S-1040 (2010). Net merit as a
#         measure of lifetime profit: 2010 revision. USDA ARS Research Report NM$4.
# Note: No new traits relative to 2006; emphasis on yield traits increased. 10 traits.
# Formula (nominal 2010$):
#   NM$2010 = 0.001(Milk) + 2.89(Fat) + 3.41(Prot) + 35(PL) - 182(SCS-3)
#             - 23(BWC) + 32(UDC) + 15(FLC) + 27(DPR) + 1(CA$)
b_2010 <- matrix(c(0.001, 2.89, 3.41, 35.00, -182.00, -23.00, 32.00, 15.00, 27.00, 1.00))
factor <- 0.890
b_2010 <- b_2010 / factor

CR_2010 <- (t(b_2010) %*% G[1:10, ]) / as.numeric(sqrt(t(b_2010) %*% P[1:10, 1:10] %*% b_2010))
CR_2010 <- (CR_2010*i)/GI
print(t(CR_2010))

resp2010 <- data.frame(
  index = "NM$ 2010",
  year = 2010,
  trait = tidx,
  response = as.numeric(CR_2010)
)
responses <- rbind(responses, resp2010)

# --- 2.7 NM$ 2006 --------------------------------------------------------------
# Source: VanRaden and Multi-State Project S-1008 (2006). Net merit as a measure
#         of lifetime profit: 2006 revision. USDA ARS Research Report NM$3 (07-06).
# Note: Milk removed from index; 9 traits. PL scaling revised.
# Formula (nominal 2006$):
#   NM$2006 = 2.70(Fat) + 3.55(Prot) + 29(PL) - 150(SCS-3) - 14(BWC)
#             + 28(UDC) + 13(FLC) + 21(DPR) + 1(CA$)
b_2006 <- matrix(c(2.70, 3.55, 29.00, -150.00, -14.00, 28.00, 13.00, 21.00, 1.00))
factor <- 0.822
b_2006 <- b_2006 / factor

# Milk (index position 1) is excluded; index starts at Fat (G row 2)
CR_2006 <- (t(b_2006) %*% G[2:10, ]) / as.numeric(sqrt(t(b_2006) %*% P[2:10, 2:10] %*% b_2006))
CR_2006 <- (CR_2006*i)/GI
print(t(CR_2006))

resp2006 <- data.frame(
  index = "NM$ 2006",
  year = 2006,
  trait = tidx,
  response = as.numeric(CR_2006)
)
responses <- rbind(responses, resp2006)

# --- 2.8 NM$ 2000 --------------------------------------------------------------
# Source: VanRaden (2000). Net merit as a measure of lifetime profit — 2000 version.
#         USDA ARS Research Report NM$1 (11-00).
# Key change: First inclusion of type composites (UDC, FLC, BWC). 8 traits.
# Formula (nominal 2000$):
#   NM$2000 = 0.018(Milk) + 2.14(Fat) + 4.76(Prot) + 28(PL) - 154(SCS-3)
#             - 14(BWC) + 29(UDC) + 15(FLC)
b_2000 <- matrix(c(0.018, 2.14, 4.76, 28.00, -154.00, -14.00, 29.00, 15.00))
factor <- 0.703
b_2000 <- b_2000 / factor

CR_2000 <- (t(b_2000) %*% G[1:8, ]) / as.numeric(sqrt(t(b_2000) %*% P[1:8, 1:8] %*% b_2000))
CR_2000 <- (CR_2000*i)/GI
print(t(CR_2000))

resp2000 <- data.frame(
  index = "NM$ 2000",
  year = 2000,
  trait = tidx,
  response = as.numeric(CR_2000)
)
responses <- rbind(responses, resp2000)

# --- 2.9 NM$ 1994 --------------------------------------------------------------
# Source: VanRaden and Wiggans (1995). Productive life evaluations: Calculation,
#         accuracy, and economic value. J. Dairy Sci. 78(3):631-638.
# Key change: First index combining yield traits with productive life and SCS. 5 traits.
# Formula (adjusted to account for 1994 yield pricing structure):
#   NM$1994 = 0.7[0.016(Milk) + 1.50(Fat) + 1.95(Prot)] + 11.67(PL) - 29.13(SCS-3)
#           = 0.0112(Milk) + 1.049(Fat) + 1.365(Prot) + 11.67(PL) - 29.13(SCS-3)
b_1994 <- matrix(c(0.0112, 1.049, 1.365, 11.67, -29.13))
factor <- 0.605
b_1994 <- b_1994 / factor

CR_1994 <- (t(b_1994) %*% G[1:5, ]) / as.numeric(sqrt(t(b_1994) %*% P[1:5, 1:5] %*% b_1994))
CR_1994 <- (CR_1994*i)/GI
print(t(CR_1994))

resp1994 <- data.frame(
  index = "NM$ 1994",
  year = 1994,
  trait = tidx,
  response = as.numeric(CR_1994)
)
responses <- rbind(responses, resp1994)

# --- 2.10 CY$ 1984 -------------------------------------------------------------
# Source: Norman (1986). Sire evaluation procedures for yield traits. NCDHIP Handbook H-1.
# Cheese Yield Index: milk, fat, and protein only. 3 traits.
# Formula (nominal 1984$):
#   CY$1984 = -0.00211(Milk) + 1.899(Fat) + 1.646(Prot)
# Note: Sign of milk coefficient is positive as entered here (see original source).
b_1984 <- matrix(c(0.00211, 1.899, 1.646))
factor <- 0.424
b_1984 <- b_1984 / factor

CR_1984 <- (t(b_1984) %*% G[1:3, ]) / as.numeric(sqrt(t(b_1984) %*% P[1:3, 1:3] %*% b_1984))
CR_1984 <- (CR_1984*i)/GI
print(t(CR_1984))

resp1984 <- data.frame(
  index = "CY$ 1984",
  year = 1984,
  trait = tidx,
  response = as.numeric(CR_1984)
)
responses <- rbind(responses, resp1984)

# --- 2.11 MFP$ 1976 ------------------------------------------------------------
# Source: Norman et al. (2010). Response to alternative genetic-economic indices
#         for Holsteins across 2 generations. J. Dairy Sci. 93(6):2695-2702.
# Milk-Fat-Protein index: milk, fat, and protein only. 3 traits.
# Formula (nominal 1976$):
#   MFP$1976 = 0.016(Milk) + 1.50(Fat) + 1.95(Prot)
b_1976 <- matrix(c(0.016, 1.50, 1.95))
factor <- 0.242
b_1976 <- b_1976 / factor

CR_1976 <- (t(b_1976) %*% G[1:3, ]) / as.numeric(sqrt(t(b_1976) %*% P[1:3, 1:3] %*% b_1976))
CR_1976 <- (CR_1976*i)/GI
print(t(CR_1976))

resp1976 <- data.frame(
  index = "MFP$ 1976",
  year = 1976,
  trait = tidx,
  response = as.numeric(CR_1976)
)
responses <- rbind(responses, resp1976)

# --- 2.12 PD$ 1971 -------------------------------------------------------------
# Source: Norman et al. (2010). J. Dairy Sci. 93(6):2695-2702.
# Predicted Difference Dollars: milk and fat only. 2 traits.
# Note on economic values: values of $0.0745 (milk) and $1.50 (fat) from
#   Norman et al. (2010) are used. Alternative values ($0.02595 and $0.83)
#   appear in the May 1971 USDA-DHIA Sire Summary List (ARS 44 Issue 231)
#   but could not be verified from a digital source.
# Formula (nominal 1971$):
#   PD$1971 = 0.0745(Milk) + 1.50(Fat)
b_1971 <- matrix(c(0.0745, 1.50))
factor <- 0.176
b_1971 <- b_1971 / factor

CR_1971 <- (t(b_1971) %*% G[1:2, ]) / as.numeric(sqrt(t(b_1971) %*% P[1:2, 1:2] %*% b_1971))
CR_1971 <- (CR_1971*i)/GI
print(t(CR_1971))

resp1971 <- data.frame(
  index = "PD$ 1971",
  year = 1971,
  trait = tidx,
  response = as.numeric(CR_1971)
)
responses <- rbind(responses, resp1971)

# ================================================================================
# SECTION 3: FIGURES
# ================================================================================

# --- 3.1 Linear trends across index versions -----------------------------------
# Fit linear regression of correlated response on index year for a subset of traits.
# Used to summarize directional trends over the 54-year period.

# Note: unit conversion from lb to kg (÷ 2.2054) is available but currently
#       commented out; activate before final manuscript submission.
# tall_response <- tall_response %>%
#   dplyr::mutate(response = ifelse(trait %in% c("Milk", "Fat", "Protein", "RFI"), response / 2.2054, response))

subset_vars <- c("Milk", "FAT", "BWC", "MEP", "RFI", "RMI")

results <- data.frame(trait = character(),
                      slope = numeric(),
                      intercept = numeric(),
                      Rsq = numeric(),
                      pval = numeric(),
                      stringsAsFactors = FALSE)

for (var in subset_vars) {
  df_var <- responses %>%
    dplyr::filter(trait == var) %>%
    dplyr::filter(!is.na(year) & !is.na(response))

  if (nrow(df_var) > 1) {
    model <- lm(response ~ year, data = df_var)
    summary_model <- summary(model)

    slope <- coef(summary_model)["year", "Estimate"]
    intercept <- coef(summary_model)["(Intercept)", "Estimate"]
    r_squared <- summary_model$r.squared
    p_value <- coef(summary_model)["year", "Pr(>|t|)"]

    results <- rbind(results, data.frame(trait = var,
                                         slope = slope,
                                         intercept = intercept,
                                         Rsq = r_squared,
                                         pval = p_value,
                                         stringsAsFactors = FALSE))
  } else {
    results <- rbind(results, data.frame(trait = var,
                                         slope = NA,
                                         intercept = NA,
                                         Rsq = NA,
                                         pval = NA,
                                         stringsAsFactors = FALSE))
  }
}

# Harmonize trait label for fat (FAT -> Fat) to match figure labels
responses$trait[responses$trait == "FAT"] = "Fat"

# --- 3.2 Figure 1: Raw correlated responses for six key traits -----------------
# Six-panel line plot (2x3 grid) showing raw correlated responses for
# Milk, Fat, RFI (top row) and MEP, RMI, RMY (bottom row) across all 12 indices.
# X-axis labels shown only on bottom row; Y-axis title shown only on leftmost panels.

split_data <- split(responses, responses$trait)
plot_order <- c("Milk", "Fat", "RFI", "MEP", "RMI", "RMY")

plots <- lapply(
  plot_order,
  function(trait) {
    sub_df <- split_data[[trait]]
    show_y_title <- trait %in% c("Milk", "MEP")  # Show y-axis title only for Milk and MEP
    show_x_labels <- trait %in% c("MEP", "RMI", "RMY")  # Show x-axis text and ticks only for MEP, RMI, RMY

    y_limits <- switch(
      trait,
      "Milk" = c(0, 120),
      "FAT" = c(0, 6),
      "BWC" = c(-2, 0),
      "RFI" = c(-20, 0),
      "MEP" = c(-0.2, 2),
      "RMI" = c(-0.2, 2),
      "RMY" = c(-0.2, 2),
      NULL
    )

    p <- ggplot(sub_df, aes(x = forcats::fct_reorder(index, year), y = response, group = 1)) +
      geom_line(size = 1, color = "#2171B5") +
      labs(
        x = if (show_x_labels) NULL else NULL,
        y = if (show_y_title) "Correlated Response" else NULL
      ) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 8),
        axis.text.x = if (show_x_labels) {
          element_text(angle = 45, hjust = 1)
        } else {
          element_blank()
        },
        axis.ticks.x = if (show_x_labels) element_line(color = "black") else element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = if (show_y_title) element_text(size = 10) else element_blank(),
        strip.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "none"
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", size = 0.4, color = "gray") +
      ggtitle(trait)

    if (!is.null(y_limits)) {
      p <- p + ylim(y_limits[1], y_limits[2])
    }

    return(p)
  })

final_plot <- patchwork::wrap_plots(plots, ncol = 3)
final_plot

