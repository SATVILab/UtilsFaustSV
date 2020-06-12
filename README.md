
<!-- README.md is generated from README.Rmd. Please edit that file -->

# faustutils

<!-- badges: start -->

<!-- badges: end -->

The goal of faustutils is to provide convenient functions to interrogate
a FAUST clustering and extract output from it.

## Installation

You can install the package from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MiguelRodo/faustutils")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(faustutils)
```

First save the directory FAUST saved results to:

``` r
# set the FAUST project 
proj_path <- usethis::proj_path('tests/testthat')
#> √ Setting active project to 'C:/Users/migue/OneDrive - University of Cape Town/Work/PhD/Code/faustutils'
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
