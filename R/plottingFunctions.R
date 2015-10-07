plotSignificance.DChIPRepResults <- function(object,
                                meanFunction = function(x){smhuber(x)$mu},
                                lfdrThresh = 0.2, ...){

    if(!("lfdr" %in% names(resultsDChIPRep(object)))){
    stop("no local fdr estimates in the object, run \"runTesting first! \" ")
    }


countsLog2 <- as.data.frame(log2(counts(DESeq2Data(object), norm = TRUE)))


sampleTable <- colData(DESeq2Data(object))

pos <-   seq(from =-unique(sampleTable$upstream)
                            , to = unique(sampleTable$downstream), by = 1)

suppressMessages(c.ggplot2 <- data.frame(pos = pos,
                        melt(countsLog2,
                                       variable.name= "sample",
                                       value.name= "PosSignal")))

### get mean per group for ratio plots
m.f <- meanFunction
m.per.group <- t(aggregate.data.frame(t(countsLog2),
                                      by = list(sampleTable$condition), FUN = m.f )[,-1])
colnames(m.per.group) <-  levels(sampleTable$condition)


idxSig <- FDRresults(object)$lfdr < lfdrThresh
m.per.group <- as.data.frame(m.per.group)
m.per.group$pos <- pos
m.per.group$significant <- idxSig

dataGG <- gather(m.per.group,
       key = "experimental_Group_and_significance",
       value = "mean_log2_counts", 1:2)

dataGG$experimental_Group_and_significance <- ifelse(dataGG$significant,
                                                      "significant",
                                                     as.character(dataGG$experimental_Group_and_significance))

dataGG$experimental_Group_and_significance <- factor(dataGG$experimental_Group_and_significance,
                                                     levels = c(levels(sampleTable$condition),
                                                     "significant"))

pl <- (ggplot(data = dataGG, 
              aes_string(x = "pos",
                         y = "mean_log2_counts",
                         color = "experimental_Group_and_significance")) +
      geom_point() +
      labs(x ="Distance from TSS (bp)",
           y="Nucleosome occupancy (mean counts per group , log2)") +
      scale_color_manual(values = c("#3366CC", "#e41a1c", "black")) +
      theme(panel.background = element_blank(), panel.grid.minor=element_blank(),
      axis.line = element_line(colour = "black", size = 0.5)))



  return(pl)
}

#' Produce a plot that colors the positions identified as significant
#'
#' This function plots the positionwise mean of the two conditions after
#' \code{\link{runTesting}}
#' has been run on a \code{DChIPRepResults} object. The points corresponding
#' to significant
#' positions are colored black in both  of the conditions.
#' The function returns the plot as a \code{ggplot2}
#' object that can be modified afterwards.
#'
#'
#' @include AllGenerics.R
#' @name plotSignificance
#' @rdname plotSignificance
#' @aliases plotSignificance plotSignificance,DChIPRepResults-method
#' @export
#' @importFrom reshape2 melt
#' @importFrom DESeq2 counts
#' @importFrom smoothmest smhuber
#' @importFrom tidyr gather spread
#'
#' @param object a \code{DChIPRepResults} object after \code{\link{runTesting}}
#' @param meanFunction a function to compute the positionwise mean
#' per group, defaults to a Huber estimator of the mean.
#' @param lfdrThresh Threshold for the local FDR
#' @param ... additional parameters for plotting (NOT YET IMPLEMENTED)
#' @return a \code{ggplot2} object
#'
#' @examples
#' data(testData)
#' dcr <- DChIPRepResults(testData)
#' dcr <- runTesting(dcr)
#' plotSignificance(dcr)
setMethod("plotSignificance", signature(object="DChIPRepResults"),
          plotSignificance.DChIPRepResults)


plotProfiles.DChIPRepResults <- function(object,
                                        meanFunction = function(x){smhuber(x)$mu},
                                        ...){

# get  normalized counts
ccc <- counts(DESeq2Data(object), normalized = TRUE)


sampleTable <- colData(DESeq2Data(object))

profilePerGroup <- t(aggregate.data.frame(t(ccc),
                                      by = list( sampleTable$condition),
                                      FUN = meanFunction )[,-1])
colnames(profilePerGroup) <-  levels(sampleTable$condition)

# take log2 for plotting and  substract the sample mean
profilePerGroup <- as.data.frame(scale(log2(profilePerGroup), scale = FALSE))

profilePerGroup$Pos <- integer(dim(profilePerGroup)[1])

profilePerGroup$Pos <-  seq(from =-unique(sampleTable$upstream)
                            , to = unique(sampleTable$downstream), by = 1)

dataGG <- gather(profilePerGroup,
           key = "experimental_Group", value = "log2_counts_centered", -3)

pl <- (ggplot(data = dataGG,
              aes_string(x = "Pos",
                         y = "log2_counts_centered",
                         color = "experimental_Group")) +
       geom_smooth(se = FALSE) +
       labs(x = "Distance from TSS (bp)",
            y = "Nucleosome occupancy (centered counts, log2)") +
       scale_color_manual(values = c("#3366CC", "#e41a1c")) +
       theme(panel.background = element_blank(),
             panel.grid.minor=element_blank(),
             axis.line = element_line(colour = "black", size = 0.5)))
 return(pl)

}

#' Produce a TSS plot of the two conditions in the data
#'
#' This function plots the positionwise mean of the log2 of the normalized
#' counts of the  two conditions
#' after \code{\link{runTesting}} has been run on a \code{DChIPRepResults}
#' object.
#'
#'
#'
#' @include AllGenerics.R
#' @name plotProfiles
#' @rdname plotProfiles
#' @aliases plotProfiles plotProfiles,DChIPRepResults-method
#' @param object a \code{DChIPRepResults} object after \code{\link{runTesting}}
#' @param meanFunction a function to compute the positionwise mean per group,
#' defaults  to a Huber estimator of the mean.
#'
#' @param ... additional parametes for plotting (NOT YET IMPLEMENTED)
#' @return a \code{ggplot2} object
#' @export
#' @importFrom  GenomicRanges colData
#' @importFrom smoothmest smhuber
#' @examples
#' data(testData)
#' dcr <- DChIPRepResults(testData)
#' dcr <- runTesting(dcr)
#' plotProfiles(dcr)
setMethod("plotProfiles", signature(object="DChIPRepResults"),
            plotProfiles.DChIPRepResults)


