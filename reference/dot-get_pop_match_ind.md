# Get indices of columns that have a specified annotation

Get indices of columns that have a specified annotation

## Usage

``` r
.get_pop_match_ind(data, pop)
```

## Arguments

- data:

  dataframe. Columns containg FAUST-population counts.

- pop:

  `named character vector`. Names specify marker and values specify
  level, e.g. c("CD4" = "-", "CD8" = "+").

## Value

`Integer vector`.
