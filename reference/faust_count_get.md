# Get counts of subsets identified by FAUST

Get counts of subsets identified by FAUST

## Usage

``` r
faust_count_get_pop(
  project_path = NULL,
  data = NULL,
  pop = NULL,
  dem_col = c("sample", "exp_unit", "tot_count", "tot_count_classified", "sampleName",
    "experimentalUnit"),
  exhaustive = FALSE,
  simplify_names = TRUE
)

faust_count_get(project_path, exhaustive = FALSE)
```

## Arguments

- project_path:

  character. FAUST project directory.

- data:

  dataframe. Dataframe where columns specify counts. If `NULL`, then it
  is read in from `project_path/faustData/faustCountMatrix.rds`.

- pop:

  `named character vector` or `list`. If a `character vector`, then
  names specify markers and values specify levels, e.g. c("CD4" = "-",
  "CD8" = "+"). Counts of all annotations will be returned that match
  these annotations separately, i.e. not summed. If a `list`, then each
  element must be a character vector as above. For a given list element,
  instead of all returning all subsets individually that have the
  correct level for the specified markers, we sum over all the subsets
  matching the specified annotation and plot the final count. If `NULL`,
  then all subsets are returned.

- dem_col:

  `character vector`. Specifies names of columns in `data` that we wish
  to keep, regardless of if they match a FAUST annotation.

- exhaustive:

  logical. If `TRUE`, then counts are taken from the exhaustive FAUST
  count matrix rather than the count matrix after excluding subsets that
  don't appear in sufficiently many experimental units. Default is
  `FALSE`.

- simplify_names:

  logical. If `TRUE`, then population names are simplified by removing
  the markers specified in `pop` from the column names. Default is
  `TRUE`.

## Value

A dataframe with columns as specified in `dem_col`, `tot_count` (total
number of cells classified for a sample), `pop` (population name) and
`count` (number of cells in population).

## Examples

``` r
if (FALSE) { # \dontrun{
# get counts for cells matching one annotation
pop <- c("CD4" = "-", "CD8" = "+")
faust_count_get_pop(project_path = "/path/to/project", pop = pop)
# get counts for cells matching either of the two annotations
pop <- list(c("CD4" = "-", "CD8" = "+"), c("CD8" = "-", "CD4" = "+"))
faust_count_get_pop(project_path = "/path/to/project", pop = pop)
} # }
```
