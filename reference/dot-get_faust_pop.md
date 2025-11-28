# Filter an expression matrix to return only cells matching one or more FAUST annotation combinations

Filter an expression matrix to return only cells matching one or more
FAUST annotation combinations

## Usage

``` r
.get_faust_pop(ex, pop, faust_ann)
```

## Arguments

- ex:

  `matrix`. Matrix containing marker expression values for `sample`.

- pop:

  `list`, `named character vector` or unnamed character vector of length
  one. If an unnamed character vector, then all cells matching that full
  FAUST annotation will be returned (e.g. "CD3+CD8+CD45RA+IFNg+IL2-TNF+"
  or, equivalently (depending on the number of sub-populations FAUST
  detects) "CD3~2~2~CD8~2~2~CD45RA~2~2~IFNg~2~2~IL2~1~2~TNF~2~2~")). If
  a named `character vector`, then all cells matching the set of marker
  levels are returned. (names are markers, and elements are levels, e.g.
  c("CD3" = 2, "CD8" = 2)). If a `list`, then each element must be a
  named `character vector`, and then all cells matching either of these
  specified sets of marker levels are returned.

- faust_ann:

  `data.frame`. FAUST annotation table read from
  \<project_path\>/faustData/sampleData//faustAnnotation.csv.

## Value

Numeric matrix.

## Details

This is effectively an inclusive OR statement across the different
population annotations.
