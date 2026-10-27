# ============================================================
# GY-only - Complete-GY Training Design
# Prediction of Unobserved Candidates
# ============================================================
#
# Manuscript:
# Imputed Grain Yield Used as a Secondary Trait in
# Multi-Trait Genomic Selection Can Improve Selection
# Accuracy Among Imputed Candidates
#
# Analysis:
# Model: GY-only, multi-environment genomic prediction
# Training design: Complete-GY
# Prediction scenario: Unobserved candidates
# Factor-analytic rank: 2 (FA2)
# Software: ASReml-R
#
# Description:
# Fits a GY-only multi-environment genomic prediction
# model using the Complete-GY training design.
#
# GY observations for selection candidates are masked.
# All other available GY observations are retained.
# No sparse-GY masking is applied.
# ============================================================

# 1. Load packages -------------------------------------------

library(tidyverse)
library(asreml)
library(here)

# 2. Load supporting functions -------------------------------

source(here("R", "Functions_earlyGS.R"))

# 3. Load input data -----------------------------------------

Pheno <- readRDS(
  here("Data", "Pheno_CompleteGY_Unobserved.rds")
)

Ginv <- readRDS(
  here("Data", "Ginv.rds")
)

val <- read.csv(
  here("Data", "selection_candidates.csv")
)

# 4. Set model parameters ------------------------------------

k <- 2  # Factor-analytic rank

# 5. Prepare phenotypic data ---------------------------------

Pheno_f <- Pheno |>
  ungroup() |>
  mutate(
    # Genotypes with genomic relationship information
    Gkeep = factor(
      if_else(
        as.character(Gen) %in% colnames(Ginv),
        as.character(Gen),
        NA_character_
      )
    ),
    
    # Genotypes without genomic relationship information
    Gdrop = factor(
      if_else(
        as.character(Gen) %in% colnames(Ginv),
        NA_character_,
        as.character(Gen)
      )
    ),
    
    # Spatial coordinates as factors
    Col = factor(Col),
    Row = factor(Row)
  ) |>
  # Retain only observed grain yield
  filter(Trait == "GY") |>
  droplevels()

# 6. Mask and standardize phenotypic observations ------------

df_seed <- Pheno_f |>
  mutate(
    # Mask GY for selection candidates
    Pheno = if_else(
      as.character(Gen) %in% as.character(val$Gen),
      NA_real_,
      Pheno
    )
  ) |>
  group_by(TraitEnv) |>
  mutate(
    # Standardize GY within each environment
    Pheno_z = as.vector(
      scale(Pheno, center = TRUE, scale = TRUE)
    ),
    Pheno_mean = mean(Pheno, na.rm = TRUE),
    Pheno_sd = sd(Pheno, na.rm = TRUE)
  ) |>
  ungroup() |>
  droplevels()

# 7. Check phenotypic availability ---------------------------

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

# 8. Fit multi-environment model -----------------------------

GYonly_CompleteGY_Unobserved_FA2 <- asreml(
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
  workspace = "4gb"
)

# 9. Update model --------------------------------------------

GYonly_CompleteGY_Unobserved_FA2 <- update_asreml(
  GYonly_CompleteGY_Unobserved_FA2,
  max_updates = 20
)

# 10. Save results -------------------------------------------

dir.create(
  here("Results"),
  showWarnings = FALSE,
  recursive = TRUE
)

saveRDS(
  GYonly_CompleteGY_Unobserved_FA2,
  here("Results", "GYonly_CompleteGY_Unobserved_FA2.rds")
)

