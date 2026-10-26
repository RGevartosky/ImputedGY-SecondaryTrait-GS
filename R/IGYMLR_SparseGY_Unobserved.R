# ============================================================
# GY + IGY (MLR) - Sparse-GY Training Design
# Prediction of Unobserved Candidates
# ============================================================
#
# Manuscript:
# Imputed Grain Yield Used as a Secondary Trait in
# Multi-Trait Genomic Selection Can Improve Selection
# Accuracy Among Imputed Candidates
#
# Analysis:
# Model: Multi-trait, multi-environment  (GY + IGY)
# Imputation method: MLR (multiple linear regression)
# Training design: Sparse-GY
# Prediction scenario: Unobserved candidates
# Factor-analytic rank: 4 (FA4)
# Software: ASReml-R
#
# Description:
# This script fits a multi-trait multi-environment genomic
# prediction model using observed and imputed grain yield.
# Phenotypes of selection candidates are masked for both
# traits. Observed GY is additionally masked according to
# the sparse-GY training design.
#
# ============================================================

# 1. Load packages -------------------------------------------

library(tidyverse)
library(asreml)
library(here)

# 2. Load supporting functions -------------------------------

source(here("R", "Functions_earlyGS.R"))

# 3. Load input data -----------------------------------------

Pheno <- readRDS(
  here("Data", "Pheno_SparseGY_Unobserved.rds")
)

Ginv <- readRDS(
  here("Data", "Ginv.rds")
)

val <- read.csv(
  here("Data", "selection_candidates.csv")
)

# 4. Set model parameters ------------------------------------

k <- 4  # Factor-analytic rank

# 5. Prepare genotype effects --------------------------------

Pheno <- Pheno |>
  mutate(
    Gkeep = factor(
      if_else(Gen %in% colnames(Ginv), Gen, NA_character_)
    ),
    Gdrop = factor(
      if_else(Gen %in% colnames(Ginv), NA_character_, Gen)
    )
  )

# 6. Filter phenotypic data ----------------------------------
# Pheno dataset contains both IGYs, remove IGY_MF
Pheno_f <- Pheno |>
  filter(!str_starts(TraitEnv, "mf")) |>
  droplevels()

# 7. Mask phenotypic observations ----------------------------

df_seed <- Pheno_f |>
  mutate(
    # Mask all traits for validation candidates
    Pheno = if_else(
      Gen %in% val$Gen, NA_real_, Pheno
    ),
    
    # Mask observed GY under the sparse-GY design
    Pheno = if_else(
      Trait == "GY" & mask_zone == TRUE,
      NA_real_,
      Pheno
    )
  ) |>
  group_by(TraitEnv) |>
  mutate(
    Pheno_z = as.vector(
      scale(Pheno, center = TRUE)
    ),
    Pheno_mean = mean(Pheno, na.rm = TRUE),
    Pheno_sd = sd(Pheno, na.rm = TRUE)
  ) |>
  droplevels()

# 8. Fit MTME model ------------------------------------------

IGYMLR_SparseGY_Unobserved_FA4 <- asreml(
  Pheno_z ~ TraitEnv,
  random = ~
    rr(TraitEnv, k):vm(Gkeep, Ginv) +
    diag(TraitEnv):vm(Gkeep, Ginv) +
    diag(TraitEnv):ide(Gkeep) +
    diag(TraitEnv):Block,
  residual = ~ dsum(
    ~ ar1(Col):ar1(Row) | TraitEnv
  ),
  sparse = ~ TraitEnv:Gdrop,
  data = df_seed,
  na.action = na.method(x = "include"),
  maxit = 13,
  workspace = "12gb"
)

# 9. Update model --------------------------------------------

IGYMLR_SparseGY_Unobserved_FA4 <- update_asreml(
  IGYMLR_SparseGY_Unobserved_FA4,
  max_updates = 20
)

# 10. Save fitted model --------------------------------------

dir.create(
  here("Results"),
  showWarnings = FALSE,
  recursive = TRUE
)

saveRDS(
  IGYMLR_SparseGY_Unobserved_FA4,
  here(
    "Results",
    "IGYMLR_SparseGY_Unobserved_FA4.rds"
  )
)
