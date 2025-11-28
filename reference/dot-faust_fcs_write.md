# Save FAUST subset as an FCS file

Save FAUST subset as an FCS file

## Usage

``` r
.faust_fcs_write(
  project_path,
  fr_source,
  sample_name,
  sel_sample,
  pop,
  dir_save,
  trans_fn = NULL,
  trans_chnl = NULL
)
```

## Arguments

- project_path:

  character. FAUST project directory.

- fr_source:

  GatingSet or character vector. If a `GatingSet`, then the flowFrames
  within the GatingSet are used to create the output fcs files. If a
  character vector, then must specify a directory containing FCS files.
  These FCS files are used to create the output fcs files.

- sample_name:

  `character vector`. Specifies names of samples in order in which they
  are found in `gs`.

- sel_sample:

  `character vector`. Character vector specifying the names of the
  samples (as saved by FAUST as folder names in the analysis map and
  equivalent to `gs[[i]]@name` where i is an index in `gs`.) that are to
  have their FAUST pops saved.

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

- dir_save:

  `character`. Specifies directory in which files are saved.

- trans_fn:

  function. If supplied, this function is applied to the expression
  data. Useful for back-transformation. If `NULL`, then no
  transformation is applied. Default is `NULL`.

- trans_chnl:

  character vector. If specified, `trans_fn` is applied to only these
  channels. If `NULL` and if `trans_fn` is not `NULL`, then `trans_fn`
  is applied to entire expression matrix. Default is `NULL`.

## Value

`invisible(TRUE)`.
