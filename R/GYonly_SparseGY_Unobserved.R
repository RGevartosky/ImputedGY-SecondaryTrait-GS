# ============================================================
# GY-only Model - Sparse-GY Training Design
# Prediction of Unobserved Candidates
# ============================================================
#
# Manuscript:
# Imputed Grain Yield Used as a Secondary Trait in
# Multi-Trait Genomic Selection Can Improve Selection
# Accuracy Among Imputed Candidates
#
# Analysis:
# Model: Grain yield only (GY-only)
# Training design: Sparse-GY
# Prediction scenario: Unobserved candidates
# Factor-analytic rank: 2 (FA2)
#
# Description:
# Fits a multi-trait multi-environment model using ASReml-R
# for grain yield across environments.
# Phenotypes of validation candidates are masked.
# Observed GY is additionally masked according to
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

k <- 2  # Factor-analytic loadings

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

# 6. Filter grain yield traits -------------------------------
# ILYT_Pheno contains multiple traits, for this analysis we will fit GY only

Pheno_f <- Pheno |>
  filter(str_starts(Trait, "GY")) |>
  droplevels()

# 7. Mask phenotypic data ------------------------------------

df_seed <- Pheno_f |>
  mutate(
    Pheno = if_else(
      Gen %in% val$Gen, NA_real_, Pheno
    ),
    Pheno = if_else(
      mask_zone == TRUE, NA_real_, Pheno
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

GYonly_SparseGY_Unobserved_FA2 <- asreml(
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

GYonly_SparseGY_Unobserved_FA2 <- update_asreml(
  GYonly_SparseGY_Unobserved_FA2,
  max_updates = 20
)

# 10. Save model ---------------------------------------------

dir.create(
  here("Results"),
  showWarnings = FALSE,
  recursive = TRUE
)

saveRDS(
  GYonly_SparseGY_Unobserved_FA2,
  here("Results", "GYonly_SparseGY_Unobserved_FA2.rds")
)