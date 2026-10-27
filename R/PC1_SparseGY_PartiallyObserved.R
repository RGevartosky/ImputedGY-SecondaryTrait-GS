# ============================================================
# GY + PC1 - Sparse-GY Training Design
# Prediction of Partially Observed Candidates
# ============================================================
#
# Manuscript:
# Imputed Grain Yield Used as a Secondary Trait in
# Multi-Trait Genomic Selection Can Improve Selection
# Accuracy Among Imputed Candidates
#
# Analysis:
# Model: Multi-trait, multi-environment (GY + PC1)
# Training design: Sparse-GY
# Prediction scenario: Partially observed candidates
# Factor-analytic rank: 5 (FA5)
# Software: ASReml-R
#
# Description:
# Fits a multi-trait multi-environment genomic prediction
# model using grain yield (GY) and PC1.
#
# Observed GY is masked for selection candidates,
# while PC1 observations are retained as a secondary trait.
# GY is additionally masked in the predefined sparse-GY
# training zones.
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
  here("Data", "Pheno_PC1.rds")
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

Pheno <- Pheno |>
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
    )
  ) |>
  droplevels()

# 6. Mask and standardize phenotypic observations ------------

df_seed <- Pheno |>
  mutate(
    # Mask observed GY for selection candidates
    Pheno = if_else(
      Trait == "GY" & Gen %in% val$Gen,
      NA_real_,
      Pheno
    ),
    
    # Mask observed GY in sparse-GY training zones
    Pheno = if_else(
      Trait == "GY" & mask_zone == TRUE,
      NA_real_,
      Pheno
    )
  ) |>
  group_by(TraitEnv) |>
  mutate(
    # Standardize within each trait-environment combination
    Pheno_z = as.vector(
      scale(Pheno, center = TRUE, scale = TRUE)
    ),
    Pheno_mean = mean(Pheno, na.rm = TRUE),
    Pheno_sd = sd(Pheno, na.rm = TRUE)
  ) |>
  ungroup() |>
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

# 7. Fit MTME model ------------------------------------------

PC1_SparseGY_PartiallyObserved_FA5 <- asreml(
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
  workspace = "8gb"
)

# 8. Update model --------------------------------------------

PC1_SparseGY_PartiallyObserved_FA5 <- update_asreml(
  PC1_SparseGY_PartiallyObserved_FA5,
  max_updates = 20
)

# 9. Save fitted model ---------------------------------------

dir.create(
  here("Results"),
  showWarnings = FALSE,
  recursive = TRUE
)

saveRDS(
  PC1_SparseGY_PartiallyObserved_FA5,
  here(
    "Results",
    "PC1_SparseGY_PartiallyObserved_FA5.rds"
  )
)
