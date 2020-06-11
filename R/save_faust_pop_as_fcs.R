#' @title Save FAUST subset as an FCS file
#'
#' @description
#'
#' @param projectPath character. Path to directory containing FAUST output.
#' @param experimentalUnit,projectPath,startingCellPop Should be set to same value as
#' arguments of the same name when \code{faust} was run. See \code{faust} documentation for details.
#' @param sample integer vector or character vector. If an integer vector, then it
#' specifies the indices of the samples for which to save output. If character,
#' then it specifies the names of the sames to save output for. If \code{NULL},
#' then the output for every sample is saved. Default is \code{NULL}.
#' @param trans function. If supplied, this function is applied to the expression data.
#' Useful for back-transformation, as the original
#' GatingSet . If \code{NULL}, then no transformation is applied. Default is \code{NULL}.
#' @param gs GatingSet. Original \code{GatingSet} object to which \code{faust} was applied.
#' @param saveFCS Boolean. If \code{TRUE}, then FCS files of the selected cell population(s)
#' and sample(s) are saved to projectPath/faustData/fcsData. Default is \code{FALSE}.
#' population of all the selected samplesare saved to projectPath/faustData/gsData. Default is \code{FALSE}.
#'
#' @return \code{invisible(TRUE)}
save_faust_pop_as_fcs <- function(projectPath,
                                  gs,
                                  cellPop,
                                  trans = NULL,
                                  sample = NULL,
                                  saveFCS = FALSE,
                                  saveGS = FALSE){

  if(!(saveFCS | saveGS)) stop("One of saveFCS and saveGS must be TRUE.")

  if(is.null(cellPop)) stop("If cellPop is NULL, then rather use the original FCS files.")

  gsNameVec <- rep("", length(gs))
  gsNameVec <- sapply(seq_along(gs),function(i) gs[[i]]@name)

  # set of cluster names
  faustClusterNames <- readRDS(paste0(projectPath,"/faustData/metaData/scampClusterNames.rds"))
  faustClusterNames <- append(faustClusterNames,"0_0_0_0_0")

  # select samples to save
  analysisMap <- readRDS(paste0(projectPath,"/faustData/metaData/analysisMap.rds"))
  activeSamples <- analysisMap[,"sampleName"]
  selSamples <- switch(typeof(sample),
                       "integer"=,
                       "double"=activeSamples[sample],
                       "NULL"=activeSamples,
                       "character"=activeSamples[vapply(activeSamples, function(x) x %in% sample, logical())],
                       stop("Incorrect specification of sample parameter.") )

  listDepth <- switch(as.character(is.list(cellPop[[1]])),
                      "TRUE" = 2,
                      "FALSE" = 1,
                      stop("Error in specification of cellPop."))

  baseFCSPath <- switch(as.character(saveFCS),
                        "TRUE" = file.path(projectPath, "faustData", "fcsData"),
                        "FALSE" = file.path(tempdir(), "faustData", "fcsData"),
                        stop("saveFCS incorrectly specified."))
  if(!dir.exists(baseFCSPath)) dir.create(baseFCSPath)
  for(sampleName in selSamples){
    # annotations for each cell for a given sample
    sAnn <- utils::read.table(file=paste0(projectPath,"/faustData/sampleData/",sampleName,"/faustAnnotation.csv"),
                              header=F,sep="`",
                              stringsAsFactors=FALSE)[,1]

    # expression dataframe for a given sample
    #exprsMat <- readRDS(paste0(projectPath,"/faustData/sampleData/",sampleName,"/exprsMat.rds"))
    gh <- gs[[which(gsNameVec==sampleName)]]
    fr <- flowWorkspace::gh_pop_get_data(gh)
    exprsMat <- exprs(fr)
    exprsDf <- as.data.frame(exprsMat)

    # data to output
    fcsList <- list()

    if(listDepth==1){
      fcsList <- append( fcsList, list( .subset(.channelData = exprsDf, .cellPop = cellPop,
                                                .quantiles = 0:1, .sAnn = sAnn) ) )
      subsetString <- ""
      for(i in seq_along(cellPop)) subsetString <- paste0(subsetString, names(cellPop)[i], cellPop[[i]])
      names(fcsList) <- subsetString
    }

    if(listDepth==2){
      for(i in seq_along(cellPop)){
        fcsList <- append(fcsList, list( .subset(.channelData = exprsMat, .cellPop = cellPop[[i]],
                                                 .quantiles = 0:1, .sAnn = sAnn) ) )
        subsetString <- ""
        for(j in seq_along(cellPop[[i]])) subsetString <- paste0(subsetString, names(cellPop[[i]])[j], cellPop[[i]][[j]])
        names(fcsList[[i]]) <- subsetString
      }
    }

    for(i in seq_along(fcsList)){
      popPath <- file.path(baseFCSPath, names(fcsList)[i])
      if(sampleName==selSamples[1] & dir.exists(popPath)) unlink(popPath, recursive=TRUE)
      if(!dir.exists(popPath)) dir.create(popPath)
      ex <- as.matrix(fcsList[[i]])
      if(nrow(ex)==0){
        colNames <- colnames(ex)
        ex <- matrix(rep(NA_integer_, ncol(ex)), ncol = ncol(ex))
        colnames(ex) <- colNames
      } else if(!is.null(trans)){
        ex <- trans(ex)
      }
      exprs(fr) <- ex
      flowCore::write.FCS(x=fr, filename=stringr::str_replace(file.path(popPath, sampleName),
                                                              "singlets_cleaned", "bt"))#,paste0(subsetString,"_",sampleName)))
    }
  }
  invisible(TRUE)
}

#' @inheritParams .subset
.selQuant <- function(.channelData = channelData, .quantiles = quantiles){
  if(!is.data.frame(.channelData)) stop (".channelData must be a dataframe.")

  if((is.null(.quantiles) | identical(.quantiles, 0:1))) return(.channelData)

  .quantiles <- sort(.quantiles)
  if(min(.quantiles)<0 | max(.quantiles)>1) stop(".quantiles must only have two values between 0 and 1 (inclusive).")

  channelLookup <- rep(TRUE, nrow(.channelData))
  for(i in seq_along(.channelData)){
    currChannelData <- .channelData[[i]]
    channelQs <- as.numeric(quantile(currChannelData, probs=.quantiles))
    channelLookup <- channelLookup & (currChannelData >= channelQs[1])
    channelLookup <- channelLookup & (currChannelData <= channelQs[2])
  }
  .channelData[channelLookup,,drop=FALSE]
}
