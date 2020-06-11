
```{r .subset}
#' @param .channelData dataframe. Must only contain columns of markers to be plotted. Default is \code{channelData} object.
#' @param .quantiles numeric vector. Quantiles specifying the lower and upper quantiles, respectively, of the expression levels to plot. If \code{NULL}, then all data points are plotted. Default is \code{quantiles} object, which is from \code{faust} unless the value has been changed.
#' @param .cellPop character vector or list. Specifies the cell population. If a character vector, then must be of length 1 and
#'corresponding to a cell population determined by the FAUST pipeline. Must be taken from the
#' projectPath/faustData/faustCountMatrix.rds column names. If a list, then the name of each element is a marker and the value specifies the marker expression level, i.e. one of "-", "+", "Dim", "Bright",  "Med-", "Med+" and "++". An example would be \code{list( "CD4" = "+", "CD3" = "+" )}. If \code{NULL}, then no subsetting is performed and the plot reflects the \code{startingCellPop}. Default is \code{NULL}.
#' @param .sAnn character vector. Each element specifies the final FAUST cluster the corresponding cell has been placed into.
.subset <- function(.channelData = channelData,
                    .cellPop = cellPop,
                    .quantiles = quantiles,
                    .sAnn = sAnn){

  # if root population is required, i.e. .subset is NULL
  if(is.null(.cellPop)) return(.selQuant(.channelData, .quantiles))
  print(.cellPop)
  # if desired population is a fully specified FAUST subset
  if(is.character(.cellPop)){
    if(length(.cellPop)>1) stop(".cellPop, if a character vector, must have length 1.")
    if(!(.cellPop %in% .sAnn)) stop(".cellPop not a cluster annotated by faust.")
    return(channelData[.sAnn==.cellPop, , drop=FALSE])
  }

  # if desired population is an aggregate across FAUST clusters
  # that have specified levels of specified markers

  # convert "natural language" names of levels to faust names of levels, if used
  for(i in seq_along(.cellPop)){
    elemVal <- .cellPop[[i]]
    reps <- c( "Med+"="~3~4~", "Med-"="~2~4~", "Dim"="~2~3~",
               "++"="~4~4~", "Bright"="~3~3~", "+"="~2~2~",
               "-"="~1", "13" = "~1~2~", "13" = "~1~3~", "14" = "~1~4~" ) # note that only "~1~ is necessary for later markelLevel matching
    if(!(elemVal %in% names(reps) | elemVal %in% reps)) stop(paste("Unknown marker level annotation:", elemVal))
    for(j in seq_along(reps)) if(elemVal == names(reps)[j]) elemVal <- setNames(reps[j], NULL)
    .cellPop[[i]] <- elemVal
  }

  # select only cells in sAnn that have specified annotation
  selInd <- rep(TRUE, length(.sAnn))
  for(i in seq_along(.cellPop)){
    m <- names(.cellPop)[i]
    level <- .cellPop[[i]]
    mLevel <- paste0(m, level)
    selInd <- selInd & grepl(mLevel, .sAnn)
  }

  # subset .channelData
  .channelData <- .channelData[selInd, , drop=FALSE]
  .selQuant(.channelData, .quantiles)
}
```

```{r saveFaust}

```

```{r }
library(flowCore)
library(flowWorkspace)
library(ncdfFlow)
```
