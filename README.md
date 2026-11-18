# LGPE-SRTM scalar-intercept reproduction workflow

This directory is a self-contained set of the four R scripts needed to run
`reproduce_figure_ex.R` with the revised fixed-effect parameterisation.

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

Set the working directory to this directory (the directory containing this
README and `reproduce_figure_ex.R`), then run:

```r
source("reproduce_figure_ex.R")
```

For example, in RStudio open `reproduce_figure_ex.R`, use **Session > Set Working
Directory > To Source File Location**, and then click **Source**.

The script sources these files automatically:

- `R/lgpe_srtm.R`
- `R/simulation.R`
- `R/lpntpet_comparison.R`

It writes the PDF figures to the working directory.


## Validation note

The scripts include explicit checks of the new prior dimensions before model
fitting. Their brackets, parentheses, file paths, and cross-file object names
have also been checked. This environment did not contain an R executable, so a
full numerical run could not be performed here.
