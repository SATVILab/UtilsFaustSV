# Save FAUST subset as an FCS file

Save a specified FAUST-identified population to
`project_path/faustData/fcsData` as an FCS file for all samples gated or
just a specified subset.

## Usage

``` r
faust_fcs_write(
  project_path,
  dir_save = NULL,
  pop,
  fr_source = NULL,
  sample = NULL,
  trans_fn = NULL,
  trans_chnl = NULL
)
```

## Arguments

- project_path:

  character. FAUST project directory.

- dir_save:

  character. Directory to save to. If `NULL` (the default), then FCS
  files are saved to
  `file.path(project_path, "faustData", "fcsData", <pop_defn>)`, where
  `<pop_defn>` is a concatenation of the population definition, e.g
  list(c("CD3" = 1, "CD4" = 2)) becomes "CD3~1~CD4~2~". If `character`,
  then the FCS files are saved directly to this directory. As stated
  before, default is `NULL`.

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

- fr_source:

  GatingSet or character vector. If a `GatingSet`, then the flowFrames
  within the GatingSet are used to create the output fcs files. If a
  character vector, then must specify a directory containing FCS files.
  These FCS files are used to create the output fcs files.

- sample:

  integer vector or character vector. If an integer vector, then it
  specifies the indices of the samples for which to save output. If
  character, then it specifies the names of the sames to save output
  for. If `NULL`, then the output for every sample is saved. Default is
  `NULL`.

- trans_fn:

  function. If supplied, this function is applied to the expression
  data. Useful for back-transformation. If `NULL`, then no
  transformation is applied. Default is `NULL`.

- trans_chnl:

  character vector. If specified, `trans_fn` is applied to only these
  channels. If `NULL` and if `trans_fn` is not `NULL`, then `trans_fn`
  is applied to entire expression matrix. Default is `NULL`.

## Value

`invisible(TRUE)`. Side effect is the saved FCS file.
