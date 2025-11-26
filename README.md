
<!-- README.md is generated from README.Rmd. Please edit that file -->

# UtilsFaustSV

<!-- badges: start -->

[![R-CMD-check](https://github.com/SATVILab/UtilsFaustSV/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SATVILab/UtilsFaustSV/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/SATVILab/UtilsFaustSV/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/SATVILab/UtilsFaustSV/actions/workflows/test-coverage.yaml)
[![Codecov test
coverage](https://codecov.io/gh/SATVILab/UtilsFaustSV/graph/badge.svg)](https://app.codecov.io/gh/SATVILab/UtilsFaustSV)
<!-- badges: end -->

## Overview

UtilsFaustSV provides utility functions to extract and work with outputs
from [FAUST](https://github.com/RGLab/FAUST) (Functional Annotation of
Unsupervised T-cell clustering) for use in downstream analyses.

The main purpose of this package is to allow flexible extraction of
FAUST-identified cell populations, including:

- **FCS files**: Export cells matching specific FAUST population
  definitions as FCS files for further analysis or visualization
- **Cell counts**: Extract counts of cells in FAUST-identified
  populations
- **Flexible population definitions**: Define populations using any
  subset of the markers used by FAUST, enabling extraction of broader
  cell lineages (e.g., all CD3+CD4+ cells) or combinations of multiple
  FAUST clusters

## Installation

You can install the package from [GitHub](https://github.com/) with:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("SATVILab/UtilsFaustSV")
```

## Key Functions

| Function                   | Description                                     |
|----------------------------|-------------------------------------------------|
| `faust_fcs_write()`        | Write FAUST-identified populations to FCS files |
| `faust_count_get()`        | Get raw count matrix from FAUST output          |
| `faust_count_get_pop()`    | Get counts for specific population definitions  |
| `faust_count_plot()`       | Plot counts of FAUST populations                |
| `faust_marker_get_usage()` | Get markers and levels used by FAUST            |

## Example Usage

``` r
library(UtilsFaustSV)
```

First, set the path to the FAUST project directory:

``` r
# set the FAUST project path
proj_path <- usethis::proj_path("tests/testthat")
#> ✔ Setting active project to "/home/runner/work/UtilsFaustSV/UtilsFaustSV".
```

### Get Markers Used by FAUST

Retrieve the markers FAUST used for clustering and the number of
expression levels for each marker:

``` r
faust_marker_get_usage(project_path = proj_path)
#>         CD33          CD7         CCR7      CD8-IgD HLA-DR-beads         CD14 
#>            2            2            2            2            2            2 
#>         CD27          CD4         CD16         CD20   TCRgd-CD19          CD3 
#>            2            2            2            2            2            2 
#>       CD45RA        CXCR5 
#>            2            2
```

### Write FAUST Populations to FCS Files

Export cells matching a population definition to FCS files. Populations
can be defined using only a subset of markers, enabling extraction of
broader cell types:

``` r
# Export CD3+CD4+ cells (defined by any number of levels of other markers)
faust_fcs_write(
  project_path = proj_path,
  pop = c("CD3" = "+", "CD4" = "+"),
  fr_source = gs  # GatingSet used by FAUST
)

# Export multiple populations at once
faust_fcs_write(
  project_path = proj_path,
  pop = list(
    c("CD3" = "+", "CD4" = "+"),
    c("CD3" = "+", "CD4" = "-", "CD8-IgD" = "+")
  ),
  fr_source = gs
)
```

### Get Population Counts

Extract counts of cells matching specific population definitions:

``` r
# Get counts for all subsets matching an annotation
faust_count_get_pop(
  project_path = proj_path,
  pop = c("CD4" = "-", "CD8-IgD" = "+")
)

# Sum counts across matching subsets for multiple populations
faust_count_get_pop(
  project_path = proj_path,
  pop = list(
    c("CD3" = "+", "CD4" = "+", "CD8-IgD" = "-"),
    c("CD3" = "+", "CD4" = "-", "CD8-IgD" = "+")
  )
)
```

### Plot Population Counts

Create visualizations of population frequencies:

``` r
# Plot all subsets matching an annotation
faust_count_plot(
  project_path = proj_path,
  pop = c("CD3" = "+", "CD4" = "+", "CD8-IgD" = "-")
)

# Plot summed counts for specified populations
faust_count_plot(
  project_path = proj_path,
  pop = list(
    c("CD3" = "+", "CD4" = "+", "CD8-IgD" = "-"),
    c("CD3" = "+", "CD4" = "-", "CD8-IgD" = "+")
  )
)
```

## Population Definition Syntax

Populations can be specified in several ways:

- **Named character vector**: `c("CD3" = "+", "CD4" = "+")`
  - Marker names as names, expression levels as values
  - Use `+`, `-`, `Dim`, `Bright`, etc. for levels
- **List of vectors**:
  `list(c("CD3" = "+", "CD4" = "+"), c("CD3" = "+", "CD8" = "+"))`
  - Combine multiple population definitions
  - Counts are summed across all matching populations
- **Numeric notation**: `c("CD3" = 2, "CD4" = 2)` where `2` = positive,
  `1` = negative
- **Full FAUST annotation**: e.g., `"CD3+CD8+CD45RA+IFNg+IL2-TNF+"`
