# ============================================================
# Predictive ability evaluation - all GS models
# ============================================================
#
# Evaluates 28 ASReml-R models:
#   GYonly, NDRE, PC1, IGYMF, IGYMLR
#   CompleteGY and SparseGY
#   Observed, Partially Observed, Unobserved
#
# FA rank is detected INSIDE each fitted model.
#
# Reference: BLUPS_NFA2.csv
# Candidates: selection_candidates.csv
#
# ============================================================

library(tidyverse)
library(asreml)
library(here)

# 1. Input and output directories ----------------------------

# Project-relative paths
# Run the script from the GitHub repository's RStudio project.

model_dir <- here("Results")

reference_path <- here(
  "Data",
  "BLUPS_NFA2.csv"
)

validation_path <- here(
  "Data",
  "selection_candidates.csv"
)

out_dir <- here(
  "Results",
  "PredictiveAbility"
)

# Check required inputs
if (!dir.exists(model_dir)) {
  stop(
    "Model directory not found: ", model_dir,
    "\nPlace fitted ASReml model files in Results/."
  )
}

if (!file.exists(reference_path)) {
  stop("Reference BLUP file not found: ", reference_path)
}

if (!file.exists(validation_path)) {
  stop("Validation file not found: ", validation_path)
}

# Create output directory
dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Project root: ", here())
message("Models: ", model_dir)
message("Reference BLUPs: ", reference_path)
message("Validation candidates: ", validation_path)
message("Results output: ", out_dir)

# 2. Identify the 28 fitted model files ----------------------

# The filename is used ONLY for scenario identification.
# FA rank is determined from the saved ASReml object.

file_pattern <- paste0(
  "^(GYonly|NDRE|PC1|IGYMF|IGYMLR)_",
  "(CompleteGY|SparseGY)_",
  "(Observed|PartiallyObserved|Unobserved)_",
  "FA[0-9]+\\.rds$"
)

files <- list.files(
  model_dir,
  pattern = file_pattern,
  full.names = FALSE
)

if (!length(files)) {
  stop("No matching ASReml models found in: ", model_dir)
}

model_specs <- tibble(filename = sort(files)) |>
  tidyr::extract(
    filename,
    into = c(
      "model_type",
      "training_id",
      "scenario_raw",
      "filename_K"
    ),
    regex = paste0(
      "^(GYonly|NDRE|PC1|IGYMF|IGYMLR)_",
      "(CompleteGY|SparseGY)_",
      "(Observed|PartiallyObserved|Unobserved)_",
      "FA([0-9]+)\\.rds$"
    ),
    remove = FALSE
  ) |>
  mutate(
    path = file.path(model_dir, filename),
    candidate_id = case_when(
      scenario_raw == "Observed" ~ "observed",
      scenario_raw == "PartiallyObserved" ~
        "partially_observed",
      scenario_raw == "Unobserved" ~ "unobserved"
    ),
    secondary_trait = recode(
      model_type,
      GYonly = "None",
      NDRE = "NDRE",
      PC1 = "PC1",
      IGYMF = "IMPGY_MF",
      IGYMLR = "IMPGY_MLR"
    ),
    TrainingSet = recode(
      training_id,
      CompleteGY = "Complete-GY TRN",
      SparseGY = "Sparse-GY TRN"
    ),
    CandidateScenario = recode(
      candidate_id,
      observed = "Observed",
      partially_observed = "Partially observed",
      unobserved = "Unobserved"
    )
  )

# A complete inventory contains 28 unique models.
if (nrow(model_specs) != 28L) {
  stop(
    "Expected 28 fitted models, but found ",
    nrow(model_specs),
    ". Inspect the model filenames."
  )
}

if (anyDuplicated(
  model_specs[c(
    "model_type", "training_id", "candidate_id"
  )]
)) {
  stop("Duplicate model configurations found.")
}

print(
  model_specs |>
    select(
      filename, training_id, candidate_id,
      secondary_trait
    ),
  n = Inf
)

# 3. Load validation genotypes -------------------------------

validation <- read_csv(
  validation_path, show_col_types = FALSE
)

stopifnot("Gen" %in% names(validation))

validation_gen <- unique(
  as.character(validation$Gen[!is.na(validation$Gen)])
)

if (!length(validation_gen)) {
  stop("No validation genotypes found.")
}

# 4. Load reference BLUPs ------------------------------------

reference <- read_csv(
  reference_path, show_col_types = FALSE
)

stopifnot(
  all(c("TraitEnv", "Gen", "blups") %in%
        names(reference))
)

reference <- reference |>
  transmute(
    TraitEnv = as.character(TraitEnv),
    TraitEnv_key = str_to_lower(TraitEnv),
    Gen = as.character(Gen),
    blups = as.numeric(blups)
  ) |>
  filter(
    str_detect(TraitEnv_key, "^gy(?:-|$)"),
    Gen %in% validation_gen
  )

if (!nrow(reference)) {
  stop("No reference BLUPs for validation candidates.")
}

if (anyDuplicated(
  reference[c("TraitEnv_key", "Gen")]
)) {
  stop("Duplicated reference genotype-environment pairs.")
}

# 5. Safe Pearson correlation -------------------------------

safe_cor <- function(x, y, min_n = 3L) {
  
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < min_n) {
    return(NA_real_)
  }
  
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(NA_real_)
  }
  
  cor(x[ok], y[ok], method = "pearson")
}

# 6. Detect FA rank from model internals ---------------------

detect_fa_rank <- function(mod) {
  
  names_vp <- names(mod$vparameters)
  
  # Example:
  # rr(TraitEnv,k):vm(Gkeep,Ginv)!GY-Neo22!fa1
  # rr(TraitEnv,k):vm(Gkeep,Ginv)!GY-Neo22!fa2
  
  loading_names <- grep(
    "^rr\\(.*\\):vm\\(Gkeep, *Ginv\\)!.*!fa[0-9]+$",
    names_vp,
    value = TRUE
  )
  
  if (!length(loading_names)) {
    stop("No genomic FA loadings found in model.")
  }
  
  factor_numbers <- as.integer(
    sub(".*!fa([0-9]+)$", "\\1", loading_names)
  )
  
  K <- max(factor_numbers)
  
  if (!identical(
    sort(unique(factor_numbers)),
    seq_len(K)
  )) {
    stop("Unexpected factor numbering in model.")
  }
  
  K
}

# 7. Extract genetic predictions for GY ----------------------

extract_gy_gebvs <- function(mod) {
  
  K <- detect_fa_rank(mod)
  
  vparams <- mod$vparameters
  coefs <- mod$coefficients$random
  cn <- rownames(coefs)
  
  if (is.null(cn)) {
    stop("ASReml random coefficients lack row names.")
  }
  
  # 7a. FA loadings
  load_i <- grep(
    "^rr\\(.*\\):vm\\(Gkeep, *Ginv\\)!.*!fa[0-9]+$",
    names(vparams)
  )
  
  load_names <- names(vparams)[load_i]
  
  te_load <- sub(
    "^.*!([^!]+)!fa[0-9]+$",
    "\\1",
    load_names
  )
  
  comp_load <- as.integer(
    sub(
      "^.*!fa([0-9]+)$",
      "\\1",
      load_names
    )
  )
  
  te <- unique(te_load)
  
  if (
    length(load_i) != length(te) * K ||
    !identical(
      comp_load,
      rep(seq_len(K), each = length(te))
    ) ||
    !identical(
      te_load,
      rep(te, times = K)
    )
  ) {
    stop("FA loadings have unexpected structure.")
  }
  
  Lam <- matrix(
    as.numeric(vparams[load_i]),
    nrow = length(te),
    ncol = K
  )
  
  # 7b. Genotype scores
  score_i <- grep(
    "^rr\\(.*\\)_Comp[0-9]+:vm\\(Gkeep, *Ginv\\)_",
    cn
  )
  
  score_names <- cn[score_i]
  
  comp_score <- as.integer(
    sub(
      "^.*_Comp([0-9]+):vm.*$",
      "\\1",
      score_names
    )
  )
  
  score_gen <- sub(
    "^.*:vm\\(Gkeep, *Ginv\\)_",
    "",
    score_names
  )
  
  gen <- score_gen[comp_score == 1L]
  
  if (
    !length(gen) ||
    length(score_i) != length(gen) * K ||
    !identical(
      comp_score,
      rep(seq_len(K), each = length(gen))
    ) ||
    !identical(
      score_gen,
      rep(gen, times = K)
    )
  ) {
    stop("Unexpected FA genotype score ordering.")
  }
  
  F <- matrix(
    as.numeric(coefs[score_i]),
    nrow = length(gen),
    ncol = K
  )
  
  # 7c. Environment-specific genomic effects
  delta_i <- grep(
    "^TraitEnv_.*:vm\\(Gkeep, *Ginv\\)_",
    cn
  )
  
  delta_names <- cn[delta_i]
  
  delta_gen <- sub(
    "^.*:vm\\(Gkeep, *Ginv\\)_",
    "",
    delta_names
  )
  
  delta_te <- sub(
    "^TraitEnv_(.*):vm\\(Gkeep, *Ginv\\)_.*$",
    "\\1",
    delta_names
  )
  
  if (
    length(delta_i) != length(te) * length(gen) ||
    !identical(
      delta_gen,
      rep(gen, times = length(te))
    ) ||
    !identical(
      delta_te,
      rep(te, each = length(gen))
    )
  ) {
    stop("Unexpected specific-effect ordering.")
  }
  
  # 7d. Reconstruct genomic breeding values
  gebv <- as.vector(
    F %*% t(Lam)
  ) + as.numeric(coefs[delta_i])
  
  predictions <- tibble(
    TraitEnv = delta_te,
    TraitEnv_key = str_to_lower(delta_te),
    Gen = delta_gen,
    gebv = gebv
  ) |>
    filter(
      str_detect(TraitEnv_key, "^gy(?:-|$)")
    )
  
  list(
    K = K,
    predictions = predictions
  )
}

# 8. Evaluate one fitted model -------------------------------

evaluate_one <- function(
    filename, model_type, training_id, scenario_raw,
    filename_K, path, candidate_id,
    secondary_trait, TrainingSet, CandidateScenario
) {
  
  message("Evaluating ", filename)
  
  mod <- readRDS(path)
  
  if (!inherits(mod, "asreml")) {
    stop("File is not an ASReml model: ", filename)
  }
  
  if (!isTRUE(mod$converge)) {
    stop("Model is not converged: ", filename)
  }
  
  extraction <- extract_gy_gebvs(mod)
  
  K <- extraction$K
  
  # Check, but do not derive K from the filename.
  if (K != as.integer(filename_K)) {
    warning(
      "Model rank differs from filename: ", filename,
      "; internal K = ", K
    )
  }
  
  predictions <- extraction$predictions |>
    filter(Gen %in% validation_gen)
  
  if (anyDuplicated(
    predictions[c("TraitEnv_key", "Gen")]
  )) {
    stop("Duplicated GEBV pairs: ", filename)
  }
  
  # Matched reference BLUPs and GEBVs
  matched <- reference |>
    inner_join(
      predictions |>
        select(TraitEnv_key, Gen, gebv),
      by = c("TraitEnv_key", "Gen")
    ) |>
    filter(
      is.finite(blups),
      is.finite(gebv)
    )
  
  if (!nrow(matched)) {
    stop("No reference/prediction matches: ", filename)
  }
  
  # 8a. Within-environment correlations
  by_env <- matched |>
    group_by(TraitEnv) |>
    summarise(
      Cor = safe_cor(blups, gebv),
      n = n(),
      .groups = "drop"
    )
  
  within_cor <- if (all(is.na(by_env$Cor))) {
    NA_real_
  } else {
    mean(by_env$Cor, na.rm = TRUE)
  }
  
  # 8b. Across-environment genotype means:
  #     original method (independent means)
  ref_means <- reference |>
    group_by(Gen) |>
    summarise(
      reference_mean = mean(blups, na.rm = TRUE),
      .groups = "drop"
    )
  
  pred_means <- predictions |>
    filter(is.finite(gebv)) |>
    group_by(Gen) |>
    summarise(
      predicted_mean = mean(gebv),
      .groups = "drop"
    )
  
  original_means <- inner_join(
    ref_means, pred_means, by = "Gen"
  )
  
  across_original <- safe_cor(
    original_means$reference_mean,
    original_means$predicted_mean
  )
  
  # 8c. Matched-environment genotype means
  matched_means <- matched |>
    group_by(Gen) |>
    summarise(
      reference_mean = mean(blups),
      predicted_mean = mean(gebv),
      .groups = "drop"
    )
  
  across_matched <- safe_cor(
    matched_means$reference_mean,
    matched_means$predicted_mean
  )
  
  # Save model metadata with all outputs
  metadata <- tibble(
    training_id = training_id,
    TrainingSet = TrainingSet,
    candidate_id = candidate_id,
    CandidateScenario = CandidateScenario,
    secondary_trait = secondary_trait,
    model_type = model_type,
    K_selected = K,
    filename = filename,
    result_source = "Fitted model"
  )
  
  summary <- tidyr::crossing(
    metadata,
    Method = c(
      "Across Env",
      "Within Env",
      "Across Env (matched)"
    )
  ) |>
    mutate(
      Cor = case_when(
        Method == "Across Env" ~ across_original,
        Method == "Within Env" ~ within_cor,
        Method == "Across Env (matched)" ~
          across_matched
      ),
      n_environments = n_distinct(
        matched$TraitEnv_key
      ),
      n_matched_pairs = nrow(matched),
      n_genotypes = n_distinct(matched$Gen)
    )
  
  by_env <- bind_cols(
    metadata[rep(1L, nrow(by_env)), ],
    by_env
  )
  
  matched <- bind_cols(
    metadata[rep(1L, nrow(matched)), ],
    matched
  )
  
  rm(mod)
  invisible(gc(verbose = FALSE))
  
  list(
    summary = summary,
    by_env = by_env,
    matched = matched
  )
}

# 9. Evaluate all 28 fitted models ---------------------------

evaluations <- purrr::pmap(
  model_specs,
  evaluate_one
)

results_fitted <- purrr::map_dfr(
  evaluations, "summary"
)

environment_results <- purrr::map_dfr(
  evaluations, "by_env"
)

matched_results <- purrr::map_dfr(
  evaluations, "matched"
)

# 10. Reuse GY-only unobserved baselines ----------------------

# GY-only has no secondary trait. Thus, for partially
# observed candidates, the corresponding unobserved
# GY-only result is the relevant baseline.
#
# These are copied comparisons, not newly fitted models.

gy_partial <- results_fitted |>
  filter(
    model_type == "GYonly",
    candidate_id == "unobserved"
  ) |>
  mutate(
    candidate_id = "partially_observed",
    CandidateScenario = "Partially observed",
    result_source = "GY-only unobserved baseline"
  )

results_all <- bind_rows(
  results_fitted,
  gy_partial
)

# 11. Prepare publication plotting data ----------------------

# This matches the columns required by the existing
# manuscript plotting script (Figures 5, 6, and 7).
#
# missing_pct is retained as a compatibility field,
# but is not calculated from the model files.
# Use TrainingSet to distinguish training designs.

plot_data <- results_all |>
  filter(
    Method %in% c("Across Env", "Within Env")
  ) |>
  mutate(
    missing_pct = NA_real_
  ) |>
  select(
    training_id,
    TrainingSet,
    candidate_id,
    CandidateScenario,
    secondary_trait,
    K_selected,
    Method,
    Cor,
    missing_pct,
    filename,
    result_source
  ) |>
  arrange(
    candidate_id,
    TrainingSet,
    Method,
    secondary_trait
  )

# Validate publication plotting inventory:
# 2 training designs x 3 scenarios x 5 trait types
# x 2 predictive-ability methods = 60 rows

if (nrow(plot_data) != 60L) {
  stop(
    "Expected 60 plotting rows; found ",
    nrow(plot_data)
  )
}

if (anyDuplicated(
  plot_data[c(
    "training_id",
    "candidate_id",
    "secondary_trait",
    "Method"
  )]
)) {
  stop("Duplicated publication plotting combinations.")
}


# 14. Exploratory visualization ------------------------------

# This plot is for checking results only.

plot_check <- plot_data |>
  filter(
    Method %in% c("Across Env", "Within Env")
  ) |>
  mutate(
    secondary_trait = recode(
      secondary_trait,
      "IMPGY_MF" = "IGY[MF]",
      "IMPGY_MLR" = "IGY[MLR]"
    ),
    secondary_trait = factor(
      secondary_trait,
      levels = c(
        "None", "NDRE", "PC1",
        "IGY[MF]", "IGY[MLR]"
      )
    ),
    TrainingSet = factor(
      TrainingSet,
      levels = c(
        "Complete-GY TRN",
        "Sparse-GY TRN"
      )
    ),
    CandidateScenario = factor(
      CandidateScenario,
      levels = c(
        "Observed",
        "Partially observed",
        "Unobserved"
      )
    ),
    Method = factor(
      Method,
      levels = c("Across Env", "Within Env")
    )
  )

p_check <- ggplot(
  plot_check,
  aes(
    x = secondary_trait,
    y = Cor,
    fill = secondary_trait
  )
) +
  geom_col(
    width = 0.75,
    color = "black",
    alpha = 0.9,
    na.rm = TRUE
  ) +
  geom_text(
    aes(
      label = if_else(
        is.na(Cor),
        "NA",
        paste0(sprintf("%.2f", Cor), "\nK=", K_selected)
      )
    ),
    vjust = -0.3,
    size = 2.8,
    na.rm = TRUE
  ) +
  facet_grid(
    CandidateScenario + TrainingSet ~ Method
  ) +
  scale_fill_manual(
    values = c(
      "None" = "#FF5F05",
      "NDRE" = "#FCB316",
      "PC1" = "#007E8E",
      "IGY[MF]" = "#1D58A7",
      "IGY[MLR]" = "#009FD4"
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0.08, 0.16))
  ) +
  coord_cartesian(
    ylim = c(0, 1.12),
    clip = "off"
  ) +
  labs(
    title = "Predictive ability across genomic selection scenarios",
    subtitle = "Exploratory comparison of all fitted models",
    x = "Secondary trait",
    y = "Pearson correlation",
    fill = "Secondary trait"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

# Display in RStudio; do not save
print(p_check)

# 12. Save all results ---------------------------------------

write_csv(
  plot_data,
  file.path(
    out_dir,
    "predictive_ability_plot_data.csv"
  )
)

write_csv(
  results_all,
  file.path(
    out_dir,
    "predictive_ability_all_metrics.csv"
  )
)

saveRDS(
  results_all,
  file.path(
    out_dir,
    "predictive_ability_all_metrics.rds"
  )
)

write_csv(
  environment_results,
  file.path(
    out_dir,
    "predictive_ability_by_environment.csv"
  )
)

saveRDS(
  matched_results,
  file.path(
    out_dir,
    "predictive_ability_matched_pairs.rds"
  )
)

write_csv(
  model_specs,
  file.path(
    out_dir,
    "model_inventory.csv"
  )
)

# 13. Display final results ----------------------------------

print(
  plot_data |>
    select(
      TrainingSet,
      CandidateScenario,
      secondary_trait,
      K_selected,
      Method,
      Cor
    ),
  n = Inf
)


