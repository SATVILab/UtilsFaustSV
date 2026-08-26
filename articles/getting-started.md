# Getting Started with UtilsFaustSV

## Introduction

UtilsFaustSV provides utility functions to extract and work with outputs
from [FAUST](https://github.com/RGLab/FAUST) (Functional Annotation of
Unsupervised T-cell clustering) for use in downstream analyses.

The main use cases for this package are:

1.  **Extracting FCS files** of cells belonging to FAUST-identified
    populations for further analysis in other tools
2.  **Getting cell counts** for populations defined using a subset of
    the markers FAUST used for clustering
3.  **Visualizing population frequencies** across samples

The key feature of this package is the ability to define populations
using only a subset of the markers that FAUST used. This allows you to:

- Extract broader cell lineages (e.g., all CD3+CD4+ cells regardless of
  other markers)
- Combine multiple FAUST clusters into a single population
- Get counts or FCS files for cell types defined at different levels of
  granularity

## Setup

``` r

library(UtilsFaustSV)

# Set the path to your FAUST project directory
project_path <- "/path/to/faust/project"
```

## Understanding FAUST Output Structure

After running FAUST, the output is stored in a directory structure like:

    project_path/
    └── faustData/
        ├── faustCountMatrix.rds         # Count matrix
        ├── exhaustiveFaustCountMatrix.rds
        ├── metaData/
        │   ├── analysisMap.rds
        │   ├── scampClusterNames.rds
        │   └── activeChannels.rds
        ├── sampleData/
        │   └── <sample_name>/
        │       └── faustAnnotation.csv  # Per-cell annotations
        └── fcsData/                     # Output FCS files (created by this package)

## Getting Marker Information

First, it’s useful to understand which markers FAUST used and how many
expression levels it identified for each:

``` r

# Get markers and their levels
markers <- faust_marker_get_usage(project_path = project_path)
print(markers)
```

This returns a named numeric vector where: - Names are the marker
names - Values are the number of expression levels (typically 2 for
binary +/-)

## Defining Populations

Populations can be defined in several ways:

### Named Character Vector

The most common way to define a population is using a named character
vector where names are markers and values are expression levels:

``` r

# CD3+CD4+ cells
pop_cd4 <- c("CD3" = "+", "CD4" = "+")

# CD3+CD8+ cells
pop_cd8 <- c("CD3" = "+", "CD4" = "-", "CD8-IgD" = "+")
```

### Expression Level Notation

For markers with more than 2 levels, you can use descriptive levels:

- `-`: Negative
- `+`: Positive (for 2-level markers)
- `Dim`: Dim expression (for 3-level markers)
- `Bright`: Bright expression (for 3-level markers)
- `MedLow`, `MedHigh`, `VeryBright`: For 4-level markers

### Numeric Notation

You can also use numeric levels (1 = lowest, 2 = next, etc.):

``` r

# Equivalent to c("CD3" = "+", "CD4" = "+")
pop_cd4_numeric <- c("CD3" = 2, "CD4" = 2)
```

### Multiple Populations

Use a list to define multiple populations:

``` r

# Both CD4+ and CD8+ T cells
pop_both <- list(
  c("CD3" = "+", "CD4" = "+", "CD8-IgD" = "-"),
  c("CD3" = "+", "CD4" = "-", "CD8-IgD" = "+")
)
```

## Writing FCS Files

The
[`faust_fcs_write()`](https://satvilab.github.io/UtilsFaustSV/reference/faust_fcs_write.md)
function exports cells matching a population definition to FCS files:

``` r

# You need the original GatingSet that FAUST used
# gs <- flowWorkspace::load_gs("/path/to/gatingset")

# Export CD3+CD4+ cells
faust_fcs_write(
  project_path = project_path,
  pop = c("CD3" = "+", "CD4" = "+"),
  fr_source = gs
)

# Files are saved to:
# project_path/faustData/fcsData/<pop_definition>/
```

### Custom Output Directory

``` r

faust_fcs_write(
  project_path = project_path,
  pop = c("CD3" = "+", "CD4" = "+"),
  fr_source = gs,
  dir_save = "/path/to/output"
)
```

### Exporting Specific Samples

``` r

# By sample name
faust_fcs_write(
  project_path = project_path,
  pop = c("CD3" = "+", "CD4" = "+"),
  fr_source = gs,
  sample = c("sample1.fcs", "sample2.fcs")
)

# By index
faust_fcs_write(
  project_path = project_path,
  pop = c("CD3" = "+", "CD4" = "+"),
  fr_source = gs,
  sample = 1:5
)
```

### Back-transforming Data

If your data was transformed before FAUST analysis, you can
back-transform when exporting:

``` r

# Back-transform data that was asinh-transformed
# (sinh is the inverse of asinh)
faust_fcs_write(
  project_path = project_path,
  pop = c("CD3" = "+", "CD4" = "+"),
  fr_source = gs,
  trans_fn = sinh
)
```

## Getting Cell Counts

### Raw Count Matrix

Get the full FAUST count matrix:

``` r

counts <- faust_count_get(project_path = project_path)
```

### Counts for Specific Populations

Get counts for populations matching a definition:

``` r

# Get counts for all CD3+CD4+ subsets
cd4_counts <- faust_count_get_pop(
  project_path = project_path,
  pop = c("CD3" = "+", "CD4" = "+")
)
```

The returned data frame includes: - Sample identifiers - Total counts
and classified counts - Counts for each matching population subset

### Summed Counts

When using a list of populations, counts are summed:

``` r

# Get summed counts for T cell subsets
t_cell_counts <- faust_count_get_pop(
  project_path = project_path,
  pop = list(
    "CD4+ T cells" = c("CD3" = "+", "CD4" = "+", "CD8-IgD" = "-"),
    "CD8+ T cells" = c("CD3" = "+", "CD4" = "-", "CD8-IgD" = "+")
  )
)
```

### Exhaustive Count Matrix

By default, FAUST filters out rare populations. To include all
populations:

``` r

all_counts <- faust_count_get_pop(
  project_path = project_path,
  pop = c("CD3" = "+", "CD4" = "+"),
  exhaustive = TRUE
)
```

## Plotting Population Counts

The
[`faust_count_plot()`](https://satvilab.github.io/UtilsFaustSV/reference/faust_count_plot.md)
function creates box plots of population frequencies:

``` r

# Plot all CD4+ T cell subsets
faust_count_plot(
  project_path = project_path,
  pop = c("CD3" = "+", "CD4" = "+", "CD8-IgD" = "-")
)
```

Plots are saved to:
`project_path/faustData/plotData/pop_stats/<population_name>.png`

### Customizing Plots

``` r

faust_count_plot(
  project_path = project_path,
  pop = c("CD3" = "+", "CD4" = "+"),
  font_size = 12,
  point_size_max = 3,
  p_width = 50,
  p_height = 20
)
```

## Workflow Example

Here’s a complete workflow for extracting and analyzing a T cell subset:

``` r

library(UtilsFaustSV)

# Set paths
project_path <- "/path/to/faust/project"

# 1. Check available markers
markers <- faust_marker_get_usage(project_path = project_path)
print(markers)

# 2. Define populations of interest
cd4_memory <- c(
  "CD3" = "+",
  "CD4" = "+",
  "CD8-IgD" = "-",
  "CD45RA" = "-",
  "CCR7" = "+"
)

# 3. Get counts
counts <- faust_count_get_pop(
  project_path = project_path,
  pop = cd4_memory
)

# 4. Create visualization
faust_count_plot(
  project_path = project_path,
  pop = cd4_memory
)

# 5. Export FCS files for further analysis
faust_fcs_write(
  project_path = project_path,
  pop = cd4_memory,
  fr_source = gs  # Your GatingSet
)
```

## Session Info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices datasets  utils     methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39       desc_1.4.3          R6_2.6.1           
#>  [4] fastmap_1.2.0       xfun_0.60           cachem_1.1.0       
#>  [7] knitr_1.51          htmltools_0.5.9     rmarkdown_2.31     
#> [10] lifecycle_1.0.5     cli_3.6.6           sass_0.4.10        
#> [13] pkgdown_2.2.1       textshaping_1.0.5   jquerylib_0.1.4    
#> [16] renv_1.0.3          systemfonts_1.3.2   compiler_4.6.1     
#> [19] tools_4.6.1         ragg_1.5.2          bslib_0.12.0       
#> [22] evaluate_1.0.5      yaml_2.3.12         otel_0.2.0         
#> [25] BiocManager_1.30.27 jsonlite_2.0.0      rlang_1.3.0        
#> [28] fs_2.1.0
```
