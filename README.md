# Imputed Grain Yield as a Secondary Trait in Genomic Selection

Reproducible R code for analyses associated with the manuscript **“Imputed Grain Yield Used as a Secondary Trait in Multi-Trait Genomic Selection Can Improve Selection Accuracy Among Imputed Candidates.”**

This repository contains multi-trait, multi-environment (MTME) genomic prediction models comparing observed grain yield (GY) alone with GY analyzed jointly with imputed GY or remotely sensed secondary traits. The analyses examine the effects of training-data availability and candidate phenotyping on genomic prediction.


## Analysis overview

### Traits and models

| Model | Responses | Factor-analytic rank used in the supplied scripts |
| --- | --- | --- |
| GY-only | Observed grain yield | FA2 |
| GY + IGY (MLR) | GY and multiple-linear-regression-imputed GY | FA4 |
| GY + IGY (MF) | GY and missForest-imputed GY | FA5 |
| GY + NDRE | GY and mean NDRE | FA5 |
| GY + PC1 | GY and the first principal component of vegetation indices (NDVI, GNDVI, and NDRE) | FA5 |


### Training set designs

- **Complete-GY:** Uses all the available GY training phenotypes.
- **Sparse-GY:** Masks GY observations in predefined zones. Within each trial, all GY observations in contiguous columns representing approximately 30% of the trial at each end are set to missing (`mask_zone`).

### Prediction scenarios

| Scenario | GY in selection candidates | Secondary trait in selection candidates |
| --- | --- | --- |
| Observed | Available where originally measured | Available |
| Partially observed | Masked | Retained where available |
| Unobserved | Masked | Masked |

In GY-only models, the secondary trait column does not apply (Figure 1). 

### Experimental design

![Experimental design showing Complete-GY and Sparse-GY training designs across observed, partially observed, and unobserved candidate scenarios](Figures/TRN_design.png)

**Figure 1.** Schematic of the phenotypic information available under the Complete-GY and Sparse-GY training designs and the three candidate-observation scenarios. **Orange** indicates observed grain yield (GY), **dark blue** indicates an available secondary trait, and **light gray** indicates missing or withheld observations. The panels represent observed (top), partially observed (middle), and unobserved (bottom) candidates. The figure is conceptual; the R scripts and input data determine which records are actually used for each analysis.

## Repository organization

```text
ImputedGY-SecondaryTrait-GS/
├── README.md
├── .gitignore
├── Figures/
│   └── TRN_design.png
├── IGY-SecondaryTrait-GS.Rproj
├── R/
│   ├── Functions_earlyGS.R
│   └── *.R                         # Model-fitting scripts
├── Data/
│   ├── Ginv.rds
│   ├── selection_candidates.csv
│   ├── Pheno_CompleteGY_Observed.rds
│   ├── Pheno_CompleteGY_Unobserved.rds
│   ├── Pheno_SparseGY_Observed.rds
│   ├── Pheno_SparseGY_Unobserved.rds
│   ├── Pheno_ndre.rds
│   └── Pheno_PC1.rds
└── Results/                        # Generated outputs 
```


## Running a model

1. Clone or download the repository, preserving the directory structure.
2. Make sure required input files are available under `Data/` and the ASReml-R license works.
3. Open `IGY-SecondaryTrait-GS.Rproj` in RStudio, or run a script from the project root:

```bash
Rscript R/GYonly_CompleteGY_Unobserved.R
```

The model scripts use `here()` for portable paths. Example output:

```text
Results/GYonly_CompleteGY_Unobserved_FA2.rds
```

Load a saved model using:

```r
model <- readRDS(here::here("Results", "GYonly_CompleteGY_Unobserved_FA2.rds"))
```


## Input data and masking

- `Ginv.rds`: inverse genomic relationship matrix used in `vm()` terms.
- `selection_candidates.csv`: genotype identifiers (`Gen`) for candidates whose phenotypes are masked under the appropriate scenarios.
- `Pheno_*.rds`: prepared phenotype tables containing trait, environment, genotype, spatial position, and response columns, as required by each model.
- `Pheno_ndre.rds` and `Pheno_PC1.rds`: prepared tables containing GY and the corresponding secondary trait.

## Citation

If you use this repository, please cite the associated manuscript:

*Gevartosky et al. (2026). Imputed Grain Yield Used as a Secondary Trait in Multi-Trait Genomic Selection Can Improve Selection Accuracy Among Imputed Candidates.*

Full publication details and DOI should be available soon.
