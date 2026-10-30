# ============================================================
# GY + NDRE - Complete-GY Training Design
# Observed-Candidate Scenario
# ============================================================
#
# Manuscript:
# Imputed Grain Yield Used as a Secondary Trait in
# Multi-Trait Genomic Selection Can Improve Selection
# Accuracy Among Imputed Candidates
#
# Analysis:
# Model: Multi-trait, multi-environment (GY + NDRE)
# Training design: Complete-GY
# Prediction scenario: Observed candidates
# Factor-analytic rank: 5 (FA5)
# Software: ASReml-R
#
# Description:
# Fits a multi-trait multi-environment genomic prediction
# model using grain yield (GY) and NDRE.
#
# All available phenotypic observations are retained.
# No additional masking is applied to training observations
# or selection candidates.
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
  here("Data", "Pheno_ndre.rds")
)

Ginv <- readRDS(
  here("Data", "Ginv.rds")
)

val <- read.csv(
  here("Data", "selection_candidates.csv")
)

# 4. Set model parameters ------------------------------------

k <- 5  # Factor-analytic rank

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
    
    # Ensure spatial coordinates are factors
    Col = factor(Col),
    Row = factor(Row)
  ) |>
  droplevels()

# 6. Standardize phenotypic observations ---------------------

df_seed <- Pheno_f |>
  group_by(TraitEnv) |>
  mutate(
    # Standardize within each trait-environment
    Pheno_z = as.vector(
      scale(Pheno, center = TRUE, scale = TRUE)
    ),
    Pheno_mean = mean(Pheno, na.rm = TRUE),
    Pheno_sd = sd(Pheno, na.rm = TRUE)
  ) |>
  ungroup() |>
  droplevels()

# 7. Check phenotypic availability --------------------------------------

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

NDRE_CompleteGY_Observed_FA5 <- asreml(
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
  na.action = na.method(
    x = "include",
    y = "include"
  ),
  maxit = 13,
  workspace = "8gb"
)

# 9. Update model --------------------------------------------

NDRE_CompleteGY_Observed_FA5 <- update_asreml(
  NDRE_CompleteGY_Observed_FA5,
  max_updates = 20
)

# 10. Save fitted model --------------------------------------

dir.create(
  here("Results"),
  showWarnings = FALSE,
  recursive = TRUE
)

saveRDS(
  NDRE_CompleteGY_Observed_FA5,
  here(
    "Results",
    "NDRE_CompleteGY_Observed_FA5.rds"
  )
)
