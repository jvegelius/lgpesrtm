# LGPE-SRTM scalar-intercept reproduction workflow

This R package fits the LGPE-SRTM model and includes a script for
reproducing the simulated PET example with the revised fixed-effect
parameterisation.

## Model parameterisation

For the reference group,

`k2a_reference(t) = beta_intercept + beta_reference_deviation(t)`,

where `beta_reference_deviation(0) = 0`. The first component is a scalar and
the remaining `m - 1` coefficients have an anchored GP prior. Each additional
group/covariate effect remains an ordinary `m`-point GP function.

## Requirements

Install the required packages once:

```r
install.packages(c("ggplot2", "dplyr", "Matrix", "mvtnorm", "patchwork"))
```

## Run

From an R session with its working directory set to the repository root,
install the package and run the bundled example:

```r
install.packages(normalizePath("."), repos = NULL, type = "source")
example <- new.env()
sys.source(system.file("reproduce", "figure_ex.R", package = "lgpesrtm",
                       mustWork = TRUE), envir = example)
print(example$p_k2a)
print(example$p_activation_beta)
```

Add `fit$beta_k2a_intercept` to `fit$beta_k2a_mat[, 1]` to reconstruct the
reference-group function. `fit$beta_k2a_mat[, 2]` is the group-effect
function.
