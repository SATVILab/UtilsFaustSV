# Check if FAUST annotation has a given level for a set of markers

Check if FAUST annotation has a given level for a set of markers

## Usage

``` r
.is_faust_ann_a_match_pop(faust_ann, pop)
```

## Arguments

- faust_ann:

  `data.frame`. FAUST annotation table read from
  \<project_path\>/faustData/sampleData//faustAnnotation.csv.

- pop:

  `named character vector`. Names are marker names and values are levels
  of marker. Values must be of the form or "\<num_1\>" or
  "\<num_1\>~\<num_2\>", where \<num_1\> is the level for the marker and
  \<num_2\> is the total number of levels for the marker.

## Details

This is effectively an AND statement across all markers specified in the
single FAUST population annotation.
