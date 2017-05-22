#' read psl gz files, assuming psl gz files don't have column header
#' @param pslFile character vector of file name(s)
#' @param toNull  character vector of column names to get rid of
#' @return data.frame, data.table of the psl table
#' @example 
readpsl <- function(pslFile, toNull=NULL) {
  stopifnot(require("data.table"))
  cols <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert",
            "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName",
            "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd",
            "blockCount", "blockSizes", "qStarts", "tStarts")
  cols.class <- c(rep("numeric",8), rep("character",2), rep("numeric",3),
                  "character", rep("numeric",4), rep("character",3))
  
  psl <- lapply(pslFile, function(f) {
    message("Reading ",f)
    data.table::fread( paste("zcat", f), sep="\t" )
  })
  psl <- data.table::rbindlist(psl)
  colnames(psl) <- cols
  
  if(length(toNull)>0) psl[, toNull] <- NULL
  
  return(as.data.frame(psl))
}