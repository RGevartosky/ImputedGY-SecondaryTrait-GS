# ============================================================
# GY + IGY (MF) - Sparse-GY Training Design
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
# Imputation method: MF (missForest)
# Training design: Sparse-GY
# Prediction scenario: Unobserved candidates
# Factor-analytic rank: 5 (FA5)
# Software: ASReml-R
#
# Description:
# This script fits a multi-trait multi-environment genomic
# prediction model using observed and imputed grain yield.
# Phenotypes of validation candidates are masked for both
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

k <- 5  # Factor-analytic rank

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

Pheno_f <- Pheno |>
  filter(!str_starts(TraitEnv, "mlr")) |>
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
  ungroup()|>
  droplevels()

observation_summary <- df_seed |>
  mutate(
    candidate = !is.na(Gen) &
      as.character(Gen) %in% as.character(val$Gen)
  ) |>
  group_by(Trait, candidate) |>
  summarise(
    n_total = n(),
    n_available = sum(!is.na(Pheno)),
    .groups = "drop"
  )

print(observation_summary)

# 8. Fit MTME model ------------------------------------------

IGYMF_SparseGY_Unobserved_FA5 <- asreml(
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

IGYMF_SparseGY_Unobserved_FA5 <- update_asreml(
  IGYMF_SparseGY_Unobserved_FA5,
  max_updates = 20
)

# 10. Save fitted model --------------------------------------

dir.create(
  here("Results"),
  showWarnings = FALSE,
  recursive = TRUE
)

saveRDS(
  IGYMF_SparseGY_Unobserved_FA5,
  here(
    "Results",
    "IGYMF_SparseGY_Unobserved_FA5.rds"
  )
)
