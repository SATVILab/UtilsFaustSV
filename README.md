
<!-- README.md is generated from README.Rmd. Please edit that file -->

# UtilsFaustSV

<!-- badges: start -->

[![R-CMD-check](https://github.com/SATVILab/UtilsFaustSV/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SATVILab/UtilsFaustSV/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/SATVILab/UtilsFaustSV/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/SATVILab/UtilsFaustSV/actions/workflows/test-coverage.yaml)
[![codecov](https://codecov.io/gh/SATVILab/UtilsFaustSV/branch/main/graph/badge.svg)](https://codecov.io/gh/SATVILab/UtilsFaustSV)
<!-- badges: end -->

The goal of UtilsFaustSV is to provide convenient functions to
interrogate a FAUST clustering and extract output from it.

## Installation

You can install the package from [GitHub](https://github.com/) with:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("SATVILab/UtilsFaustSV")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(UtilsFaustSV)
```

First save the directory FAUST saved results to:

``` r
# set the FAUST project path
proj_path <- usethis::proj_path('tests/testthat')
#> v Setting active project to 'C:/Users/migue/Work/Packages/UtilsFAUSTSV'
```

You can get the names of the markers FAUST used to cluster and the
number of levels for each (after FAUST has clustered):

``` r
get_faust_markers_and_levels(project_path = proj_path)
#>         CD33          CD7         CCR7      CD8-IgD HLA-DR-beads         CD14 
#>            2            2            2            2            2            2 
#>         CD27          CD4         CD16         CD20   TCRgd-CD19          CD3 
#>            2            2            2            2            2            2 
#>       CD45RA        CXCR5 
#>            2            2
```

You can save FAUST-identified clusters as FCS files using only a subset
of the markers (not run):

``` r
save_faust_pop(project_path = proj_path, pop = list("CD3" = 2, "CD4" = 2), 
               gs = gs) # set gs to be the GatingSet FAUST clustered
```

You can get a table of counts of FAUST-identified clusters defined using
only a subset of the markers used (not run):

``` r
get_pop_counts(project_path = proj_path, 
               pop = c("CD4" = "-", "CD8" = "+"))
```

You can save plots of counts of FAUST-identified clusters, where the
clusters are defined using only a subset of the markers (not run):

``` r
# plot all subsets that match these annotations individually
pop <- c("CD3" = "+", "CD4" = "+", "CD8-IgD" = "-")
plot_faust_count(project_path = proj_path,
                 pop = pop)

# plot total counts of subsets that match these annotations
# within each list element
pop <- list(c("CD3" = "+", "CD4" = "+", "CD8-IgD" = "-"),
           c("CD3" = "+", "CD4" = "-", "CD8-IgD" = "+"))
plot_faust_count(project_path = proj_path,
                 pop = pop)
```
