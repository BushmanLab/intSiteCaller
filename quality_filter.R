#' Quality filter read alignments
#' @param alignments GRanges object of alignment ranges with alignment metrics
#' of "matches", "repMatches", "qSize", "qStart", and "tBaseInsert".
#' @param qStartMax integer Maximum allowable nucleotide position to start read
#' alignment.
#' @param globalIdentityMin numeric Between 0 and 100 denoting the percent of 
#' global identity for the alignment needed to pass the filter.
#' @param baseInsertMax integer Maximum number of allowable inserted nucleotides
#' within the alignment.
#' @return Subset of input GRanges Object passing filter criteria
#' 
qualityFilter <- function(alignments, qStartMax = NULL, 
                          globalIdentityMin = NULL, baseInsertMax = 5){
  stopifnot(class(alignments) == "data.frame")
  statsNeeded <- c("matches", "repMatches", "qSize", "qStart", "tBaseInsert")
  stopifnot(all(statsNeeded %in% names(alignments)))
  
  alignments$percIdent <- 100*(alignments$matches +
                                 alignments$repMatches)/alignments$qSize
  if(length(qStartMax) > 0){
    alignments <- subset(alignments, alignments$percIdent >= globalIdentityMin)}
  if(length(globalIdentityMin) > 0){
    alignments <- subset(alignments, alignments$qStart <= qStartMax)}
  if(length(baseInsertMax) > 0){
    alignments <- subset(alignments, alignments$tBaseInsert <= baseInsertMax)}
  alignments
}