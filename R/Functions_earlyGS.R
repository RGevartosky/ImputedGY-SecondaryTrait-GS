# Functions

# Update ASReml-R ----
# This function updates ASReml-R until it converges
update_asreml <- function(mod, max_updates = 100) {
  count <- 0
  
  while (!mod$converge && count < max_updates) {
    count <- count + 1
    mod <- update(mod)
    
    # Print model update information
    message("Update ", count, ": Convergence = ", mod$converge)
    
    # Print LogLik value
    loglik <- mod$trace |>
      as.data.frame() |>
      tibble::rownames_to_column("Iteration") |>
      dplyr::filter(Iteration == "LogLik")
    
    print(loglik)
  }
  
  if (mod$converge) {
    message("Model successfully converged!")
  } else {
    message("Maximum updates reached. Model did not converge.")
  }
  
  return(mod)
}

# % Va----
# Percentage of additive genetic variance explained by the model

VaPct <- function(mod,k,TE_fct) {
  # variance parameters
  vparams <- mod$vparameters
  # latent environmental covariates (loadings)
  Lam <- matrix(vparams[grep('^rr.*vm.*fa', names(vparams))], ncol = k)
  # specific variances
  Psi <- diag(vparams[grep(paste0('^',TE_fct,'.*vm'), names(vparams))])
  
  # % var explained:
  # mean of ratios
  meanVaPct <- mean(diag(Lam %*% t(Lam))/diag(Lam %*% t(Lam)+Psi)) * 100
  
  # for each trait x env combo
  TE <- names(vparams)[grep(paste0('^',TE_fct,'.*vm'), names(vparams))] |>
    sub(paste0('.*',TE_fct,'_'), '', x=_)
  
  TraitEnv_VaPct <- data.frame(TE_fct = TE,
                               VaPct = diag(Lam %*% t(Lam))/diag(Lam %*% t(Lam)+Psi) * 100)
  
  return(list(meanVaPct = meanVaPct, TraitEnv_VaPct = TraitEnv_VaPct))
}

# BLUPs ----
# Get blups (no relationship matrix)
blup_asreml <- function(mod,k,TE_fct,G_fct) {
  # variance parameters
  vparams <- mod$vparameters
  # latent environmental covariates (loadings)
  Lam <- matrix(vparams[grep(paste0('^rr.*',G_fct,'.*fa'), names(vparams))], ncol = k)
  # specific variances
  Psi <- diag(vparams[grep(paste0('^',TE_fct,'.*',G_fct), names(vparams))])
  
  # coefficients
  coefs <- mod$coefficients$random
  # genotype scores (slopes)
  f <- coefs[grep('Comp', rownames(coefs))]
  Lamf <- c(matrix(f, ncol = k) %*% t(Lam)) # common GET effects
  # genetic regressions residuals
  delta <- coefs[grep(paste0('^',TE_fct,'.*',G_fct), rownames(coefs))] # specific GET effects
  # total GET effects
  Lamfdelta <- Lamf + delta
  
  # get blup for each genotype, by trait-environment combination
  
  # vector of GET
  effects_all <- rownames(coefs)
  head(effects_all)
  effects_GET <- effects_all[grepl(paste0('^',TE_fct,'.*',G_fct), effects_all)]
  head(effects_GET)
  GET <- sub(paste0(TE_fct,'_'), '', effects_GET)
  GET <- sub(paste0(':', G_fct), '', GET)
  head(GET)
  
  split_GET <- do.call(rbind, strsplit(GET, '_', fixed = TRUE))
  
  blup <- data.frame(TE_fct=split_GET[,1],
                     G_fct=split_GET[,2],
                     blup=Lamfdelta)
  rownames(blup) <- GET
  return(blup)
}


# GEBVs ----
# Get genetic estimated breeding values (GEBVs)

gebvs_asreml <- function(mod,k,TE_fct) {
  # variance parameters
  vparams <- mod$vparameters
  # latent environmental covariates (loadings)
  Lam <- matrix(vparams[grep('^rr.*vm.*fa', names(vparams))], ncol = k)
  # specific variances
  Psi <- diag(vparams[grep(paste0('^',TE_fct,'.*vm'), names(vparams))])
  
  # coefficients
  coefs <- mod$coefficients$random
  # genotype scores (slopes)
  f <- coefs[grep('Comp', rownames(coefs))]
  Lamf <- c(matrix(f, ncol = k) %*% t(Lam)) # common GET effects
  # genetic regressions residuals
  delta <- coefs[grep(paste0('^',TE_fct,'.*vm'), rownames(coefs))] # specific GET effects
  # total GET effects
  Lamfdelta <- Lamf + delta
  
  # get gebv for each genotype, by trait-environment combination
  
  # vector of GET
  effects_all <- rownames(coefs)
  head(effects_all)
  effects_GET <- effects_all[grepl(paste0('^',TE_fct,'.*vm'), effects_all)]
  head(effects_GET)
  GET <- sub(paste0(TE_fct,'_'), '', effects_GET)
  GET <- sub(':vm\\(.*?\\)', '', GET)
  head(GET)
  
  split_GET <- do.call(rbind, strsplit(GET, '_', fixed = TRUE))
  
  blup <- data.frame(TE_fct=split_GET[,1],
                     G_fct=split_GET[,2],
                     blup=Lamfdelta)
  rownames(blup) <- GET
  return(blup)
}

# Genetic value ----
# Get total genetic value

genval_asreml <- function(mod,k,TE_fct) {
  # variance parameters
  vparams <- mod$vparameters
  # latent environmental covariates (loadings)
  Lam <- matrix(vparams[grep('^rr.*vm.*fa', names(vparams))], ncol = k)
  # specific variances
  Psi <- diag(vparams[grep(paste0('^',TE_fct,'.*vm'), names(vparams))])
  
  # coefficients
  coefs <- mod$coefficients$random
  # genotype scores (slopes)
  f <- coefs[grep('Comp', rownames(coefs))]
  Lamf <- c(matrix(f, ncol = k) %*% t(Lam)) # common GET effects
  # genetic regressions residuals
  delta <- coefs[grep(paste0('^',TE_fct,'.*vm'), rownames(coefs))] # specific GET effects
  # non-additive effects
  naeff <- coefs[grep(paste0('^',TE_fct,'.*ide'), rownames(coefs))]
  # total GET effects
  Lamfdelta <- Lamf + delta + naeff
  
  # get gebv for each genotype, by trait-environment combination
  
  # vector of GET
  effects_all <- rownames(coefs)
  head(effects_all)
  effects_GET <- effects_all[grepl(paste0('^',TE_fct,'.*vm'), effects_all)]
  head(effects_GET)
  GET <- sub(paste0(TE_fct,'_'), '', effects_GET)
  GET <- sub(':vm\\(.*?\\)', '', GET)
  head(GET)
  
  split_GET <- do.call(rbind, strsplit(GET, '_', fixed = TRUE))
  
  genval <- data.frame(TE_fct=split_GET[,1],
                     G_fct=split_GET[,2],
                     genval=Lamfdelta)
  rownames(genval) <- GET
  return(genval)
}

# Gen Corr ----
gcorr_asreml <- function(mod,k,data,TE_fct) {
  # variance parameters
  vparams <- mod$vparameters
  # latent environmental covariates (loadings)
  Lam <- matrix(vparams[grep('^rr.*vm.*fa', names(vparams))], ncol = k)
  # Rotate
  # Perform Singular Value Decomposition (SVD) for rotation
  svd <- svd(Lam) # Perform SVD on the loadings matrix
  # Compute rotated estimated loadings
  LamStar <- Lam %*% svd$v
  # specific variances
  Psi <- diag(vparams[grep(paste0('^',TE_fct,'.*vm'), names(vparams))])
  # Variance covariance matrix
  Gvar <- LamStar%*%t(LamStar)+Psi
  rownames(Gvar) <- levels(data$TraitEnv)
  colnames(Gvar) <- levels(data$TraitEnv)
  # Genetic correlation matrix
  Cmat <- cov2cor(Gvar)
  return(list(Gvar = Gvar, gcorr = Cmat))
}

# Predictive ability ----
# Function to compute GS predictive ability
pred_ability <- function(model_info, blup, TST, calc_GTEacc = F, k = 5) {
  GTEacc_list <- list()
  Gacc_list <- list()
  Conv_list <- list()
  TotGacc_list <- list()
  
  for (model_name in names(model_info)) {
    load(model_info[[model_name]])
    model_obj <- get(model_name)
    
    Conv_tmp <- model_obj$converge
    
    gebv <- gebvs_asreml(model_obj, k, 'TraitEnv') |>
      rename(gebv = blup) |>
      separate(TE_fct, c('Trait','Year','Loc'), sep = '-', remove = FALSE) |>
      filter(Loc %in% c('Adv', 'Neo', 'Stp', 'Urb'))
    
    genval <- genval_asreml(model_obj, k, 'TraitEnv') |>
      separate(TE_fct, c('Trait','Year','Loc'), sep = '-', remove = FALSE) |>
      filter(Loc %in% c('Adv', 'Neo', 'Stp', 'Urb'))
    
    if (calc_GTEacc) {
      GTEacc_tmp <- gebv |>
        left_join(blup) |>
        filter(G_fct %in% TST) |>
        group_by(Trait) |>
        summarise(GTEcor = cor(gebv, blup, use = "complete.obs"), .groups = "drop") |>
        mutate(Model = model_name)
      GTEacc_list[[model_name]] <- GTEacc_tmp
    }
    
    Gacc_tmp <- gebv |>
      filter(G_fct %in% TST) |>
      group_by(Trait, G_fct) |>
      summarise(G_gebv = mean(gebv), .groups = "drop") |>
      left_join(
        blup |>
          filter(G_fct %in% TST) |>
          group_by(Trait, G_fct) |>
          summarise(G_blup = mean(blup), .groups = "drop"),
        by = c("Trait", "G_fct")
      ) |>
      group_by(Trait) |>
      summarise(Gcor = cor(G_gebv, G_blup, use = "complete.obs"), .groups = "drop") |>
      mutate(Model = model_name)
    
    TotGacc_tmp <- genval |>
      filter(G_fct %in% TST) |>
      group_by(Trait, G_fct) |>
      summarise(G_genval = mean(genval), .groups = "drop") |>
      left_join(
        blup |>
          filter(G_fct %in% TST) |>
          group_by(Trait, G_fct) |>
          summarise(G_blup = mean(blup), .groups = "drop"),
        by = c("Trait", "G_fct")
      ) |>
      group_by(Trait) |>
      summarise(TotGcor = cor(G_genval, G_blup, use = "complete.obs"), .groups = "drop") |>
      mutate(Model = model_name)
    
    Gacc_list[[model_name]] <- Gacc_tmp
    TotGacc_list[[model_name]] <- TotGacc_tmp
    Conv_list[[model_name]] <- Conv_tmp
  }
  
  Conv_summary <- bind_cols(Conv_list) |> pivot_longer(cols = everything(), names_to = 'model', values_to = 'Conv')
  Gacc_summary <- bind_rows(Gacc_list) |> arrange(Trait)
  TotGacc_summary <- bind_rows(TotGacc_list) |> arrange(Trait)
  
  if (calc_GTEacc) {
    GTEacc_summary <- bind_rows(GTEacc_list) |> arrange(Trait)
  } else {
    GTEacc_summary <- NULL
  }
  
  return(list(
    Conv_summary = Conv_summary,
    GTEacc_summary = GTEacc_summary,
    Gacc_summary = Gacc_summary,
    TotGacc_summary = TotGacc_summary
  ))
}

