# Get a searchable version of the marker levels

Get a searchable version of the marker levels

## Usage

``` r
.collapse_pop(pop, search = FALSE)
```

## Arguments

- pop:

  `list` or `named character vector`. If a `character vector`, then
  counts of all subsets matching the set of marker levels are plotted.
  If a `list`, then each list element must be a character vector. In
  that case, one boxplot of frequencies is plotted for each element of
  the list of where the frequency is the sum of all cells that match the
  population specified. If `NULL`, then all subsets found are plotted.

- search:

  logical. If `TRUE`, then square brackets are concatenated around each
  of the elements in the level strings.
