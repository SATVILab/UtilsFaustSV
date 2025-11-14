pkgname <- "UtilsFaustSV"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('UtilsFaustSV')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("dot-collapse_pop")
### * dot-collapse_pop

flush(stderr()); flush(stdout())

### Name: .collapse_pop
### Title: Get a searchable version of the marker levels
### Aliases: .collapse_pop

### ** Examples

.collapse_pop(c("CD4" = "-", "CD3" = "-"), search = FALSE)
.collapse_pop(c("CD4" = "-", "CD3" = "-"), search = TRUE)



cleanEx()
nameEx("faust_count_get")
### * faust_count_get

flush(stderr()); flush(stdout())

### Name: faust_count_get_pop
### Title: Get counts of subsets identified by FAUST
### Aliases: faust_count_get_pop faust_count_get

### ** Examples

# get counts for cells matching one annotation
pop <- c("CD4" = "-", "CD8" = "+")
get_pop_counts(pop = pop)
# get counts for cells matching either of the two annotations
pop <- list(c("CD4" = "-", "CD8" = "+"), c("CD8" = "-", "CD4" = "+"))
get_pop_counts(pop = pop)



cleanEx()
nameEx("faust_count_plot")
### * faust_count_plot

flush(stderr()); flush(stdout())

### Name: faust_count_plot
### Title: Plot counts of FAUST pops by sample
### Aliases: faust_count_plot

### ** Examples

project_path <- usethis::proj_path(ext = "/inst/extdata")
# plot all subsets of pop, as it is a character vector
pop <- c(
  "CD3" = "+", "CD4" = "+", "CD8-IgD" = "-",
  "CD20" = "-", "CD33" = "-", "CD14" = "-",
  "TCRgd-CD19" = "-"
)
plot_faust_count(
  project_path = project_path,
  pop = pop
)
# plot counts of cells matching annotation of each
# element in pop, as pop is a list
pop <- list(
  c(
    "CD3" = "+", "CD4" = "+", "CD8-IgD" = "-",
    "CD20" = "-", "CD33" = "-", "CD14" = "-",
    "TCRgd-CD19" = "-"
  ),
  c(
    "CD3" = "+", "CD4" = "-", "CD8-IgD" = "+",
    "CD20" = "-", "CD33" = "-", "CD14" = "-",
    "TCRgd-CD19" = "-"
  )
)
plot_faust_count(
  project_path = project_path,
  pop = pop
)



cleanEx()
nameEx("faust_fcs_write")
### * faust_fcs_write

flush(stderr()); flush(stdout())

### Name: faust_fcs_write
### Title: Save FAUST subset as an FCS file
### Aliases: faust_fcs_write

### ** Examples

faust_fcs_write(
  project_path = "", pop = list("CD3" = 2),
  gs = gs, sample = 1
)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
