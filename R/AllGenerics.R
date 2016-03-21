# @rdname DESeq2Data
# @export
setGeneric("DESeq2Data",
    function(object, ...) standardGeneric("DESeq2Data"))

# @rdname DESeq2Data
# @export
setGeneric("DESeq2Data<-",
    function(object, ..., value) standardGeneric("DESeq2Data<-"))


setGeneric("runTesting",
    function(object, lfcThreshold = 0.05, plotFDR = FALSE, ...)
    standardGeneric("runTesting"))

setGeneric("FDRresults",
    function(object, ...) standardGeneric("FDRresults"))


setGeneric("FDRresults<-",
    function(object, ..., value) standardGeneric("FDRresults<-"))

setGeneric("resultsDChIPRep",
    function(object, ...) standardGeneric("resultsDChIPRep"))


setGeneric("resultsDChIPRep<-",
    function(object, ..., value) standardGeneric("resultsDChIPRep<-"))

setGeneric("plotSignificance",
    function(object, ...) standardGeneric("plotSignificance"))

setGeneric("plotProfiles",
    function(object, ...) standardGeneric("plotProfiles"))
