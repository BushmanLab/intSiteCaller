#    This source code file is a component of the larger INSPIIRED genomic analysis software package.
#    Copyright (C) 2016 Frederic Bushman
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

## load hiReadsProcessor.R
libs <- c("plyr", "BiocParallel", "Biostrings", "GenomicAlignments" ,
          "hiAnnotator" ,"sonicLength", "GenomicRanges", "BiocGenerics", 
          "ShortRead", "GenomicRanges", "igraph", "data.table")
null <- suppressMessages(sapply(libs, library, character.only=TRUE))

codeDir <- get(load("codeDir.RData"))

stopifnot(file.exists(file.path(codeDir, "hiReadsProcessor.R")))
source(file.path(codeDir, "hiReadsProcessor.R"))
source(file.path(codeDir, "standardization_based_on_clustering.R"))
source(file.path(codeDir, "read_psl_files.R"))
source(file.path(codeDir, "quality_filter.R"))

#' find reads originating from vector
#' @param vectorSeq vector sequence fasta file
#' @param primerLTR primer and LTR sequence
#' @param reads.p DNAStringSet, reads on primer side
#' @param reads.l DNAStringSet, reads on linker side
#' @return character, qNames for the vector reads
#' @example findVectorReads(vectorSeq, reads.p, reads.l)
findVectorReads <- function(vectorSeq, primerLTR="GAAAATCTCTAGCA",
                            reads.p, reads.l,
                            debug=FALSE) {
  
  Vector <- readDNAStringSet(vectorSeq)
  
  primerInVector <- matchPattern(pattern=primerLTR,
                                 subject=DNAString(as.character(Vector)),
                                 algorithm="auto",
                                 max.mismatch=4,
                                 with.indels=TRUE,
                                 fixed=TRUE)
  
  #print(primerInVector)
  if( length(primerInVector)<1 ) message("--- Cannot locate primer and ltrBit in vector ---")
  
  
  globalIdentity <- 0.75
  blatParameters <- c(minIdentity=70, minScore=15, stepSize=3, 
                      tileSize=8, repMatch=112312, dots=1000, 
                      q="dna", t="dna", out="psl")
  
  
  hits.v.p <- try(readpsl(blatSeqs(query=reads.p, subject=Vector,     
                                   blatParameters=blatParameters, parallel=F)))
  if( class(hits.v.p) == "try-error" ) hits.v.p <- data.frame()
  if ( debug ) save(hits.v.p, file="hits.v.p.RData")    
  
  hits.v.l <- try(readpsl(blatSeqs(query=reads.l, subject=Vector, 
                                   blatParameters=blatParameters, parallel=F)))
  if( class(hits.v.l) == "try-error" ) hits.v.l <- data.frame()
  if ( debug ) save(hits.v.l, file="hits.v.l.RData")    
  
  ## Sometimes the vector files received from collaborators are different from the 
  ## vector put in human host. So, it is not feasible to put a lot of constrains.
  ## Filtering on globalIdentity identities for both R1 and R2 seems to work well. 
  hits.v.p <- dplyr::filter(hits.v.p, ##tStart  > ltrpos &
                            ##tStart  < ltrpos+nchar(primerLTR)+10 &
                            ##strand == "+" &
                            matches > globalIdentity*qSize &
                              qStart  <= 5 ) 
  hits.v.l <- dplyr::filter(hits.v.l, matches>globalIdentity*qSize )
  ##strand=="-") 
  hits.v <- try(merge(hits.v.p[, c("qName", "tStart")],
                      hits.v.l[, c("qName", "tStart")],
                      by="qName")
                ,silent = TRUE)
  if( class(hits.v) == "try-error" ) hits.v <- data.frame()
  
  ##hits.v <- dplyr::filter(hits.v, tStart.y >= tStart.x &
  ##                                tStart.y <= tStart.x+2000)
  
  if ( debug ) {
    save(reads.p, file="reads.p.RData")
    save(reads.l, file="reads.l.RData")
  }
  
  vqName <- unique(hits.v$qName)
  
  message(paste0("Vector sequences found ", length(vqName)))
  
  return(vqName)
}
## vqName <- findVectorReads(vectorSeq, reads.p, reads.l)


#' as tittled make PairwiseAlignmentsSingleSubject easily accessable
#' as needed by other functions
PairwiseAlignmentsSingleSubject2DF <- function(PA, shift=0) {
  stopifnot("PairwiseAlignmentsSingleSubject"  %in% class(PA))
  
  return(data.frame(
    width=width(pattern(PA)),
    score=score(PA),
    mismatch=width(pattern(PA))-score(PA),
    start=start(pattern(PA))+shift,
    end=end(pattern(PA))+shift
  ))
}


#' subset and substring
#' trim primer and ltrbit off of ltr side of read, R2 in protocol
#' both primer and ltrbit are required, otherwise disgard it
#' allow 2 mismatch for either primer or ltrbit
#' runSeq <- sapply(1:10000, function(i)
#'                  paste(sample(c("A","C","G","T"), 8, replace=TRUE),
#'                        collapse=""))
#' runSeq.p <- pairwiseAlignment(pattern=runSeq,
#'                               subject=primer,
#'                               substitutionMatrix=submat1,
#'                               gapOpening = 0,
#'                               gapExtension = 1,
#'                               type="overlap")
#' runSeq.p.df <- PairwiseAlignmentsSingleSubject2DF(runSeq.p)
#' table(runSeq.p.df$score)
#'   1    2    3    4    5    6    7
#' 211 3338 4330 1751  344   25    1
#' false positive rate 0.0025 and thus maxMisMatch=2,
#' (1/4)^(7-maxMisMatch)*choose(7-maxMisMatch) as expected
#' false positive rate combining both primer and ltr is 0.0025*0.0025=6.25E-6
#' @param reads.p DNAStringSet of reads, normally R2
#' @param primer character string of lenth 1, such as "GAAAATC"
#' @param ltrbit character string of lenth 1, such as "TCTAGCA"
#' @return DNAStringSet of reads with primer and ltr removed
#' 
trim_Ltr_side_reads <- function(reads.p, primer, ltrbit, maxMisMatch=0) {
  
  stopifnot(class(reads.p) %in% "DNAStringSet")
  stopifnot(!any(duplicated(names(reads.p))))
  stopifnot(length(primer)==1)
  stopifnot(length(ltrbit)==1)
  
  ## allows gap, and del/ins count as 1 mismatch
  submat1 <- nucleotideSubstitutionMatrix(match=1,
                                          mismatch=0,
                                          baseOnly=TRUE)
  
  ## p for primer
  ## search for primer from the beginning
  aln.p <- pairwiseAlignment(pattern=subseq(reads.p, 1, 1+nchar(primer)),
                             subject=primer,
                             substitutionMatrix=submat1,
                             gapOpening = 0,
                             gapExtension = 1,
                             type="overlap")
  aln.p.df <- PairwiseAlignmentsSingleSubject2DF(aln.p)
  
  ## l for ltrbit
  ## search for ltrbit fellowing primer
  ## note, for SCID trial, there are GGG between primer and ltr bit and hence 5
  ## for extra bases
  aln.l <- pairwiseAlignment(pattern=subseq(reads.p, nchar(primer)+1, nchar(primer)+nchar(ltrbit)+1),
                             subject=ltrbit,
                             substitutionMatrix=submat1,
                             gapOpening = 0,
                             gapExtension = 1,
                             type="overlap")
  aln.l.df <- PairwiseAlignmentsSingleSubject2DF(aln.l, shift=nchar(primer))
  
  goodIdx <- (aln.p.df$score >= nchar(primer)-maxMisMatch &
                aln.l.df$score >= nchar(ltrbit)-maxMisMatch)
  
  reads.p <- subseq(reads.p[goodIdx], aln.l.df$end[goodIdx]+1)
  
  return(reads.p)
}
##trim_Ltr_side_reads(reads.p, primer, ltrbit)


#' subset and substring
#' trim primerID linker side of read, R1 in protocol
#' a primerIDlinker has N's in the middle
#' allow 3 mismatches for either part before and after Ns
#' see reasonning above
#' @param reads.l DNAStringSet of reads, normally R1
#' @param linker character string of lenth 1, such as
#'               "AGCAGGTCCGAAATTCTCGGNNNNNNNNNNNNCTCCGCTTAAGGGACT"
#' @param maxMisMatch=3
#' @return list of read.l and primerID
#' 
trim_primerIDlinker_side_reads <- function(reads.l, linker, maxMisMatch=3) {
  
  stopifnot(class(reads.l) %in% "DNAStringSet")
  stopifnot(!any(duplicated(names(reads.l))))
  stopifnot(length(linker)==1)
  
  pos.N <- unlist(gregexpr("N", linker))
  len.N <- length(pos.N)
  link1 <- substr(linker, 1, min(pos.N)-1)
  link2 <- substr(linker, max(pos.N)+1, nchar(linker))
  
  ## allows gap, and del/ins count as 1 mismatch
  submat1 <- nucleotideSubstitutionMatrix(match=1,
                                          mismatch=0,
                                          baseOnly=TRUE)
  
  ## search at the beginning for 1st part of linker
  aln.1 <- pairwiseAlignment(pattern=subseq(reads.l, 1, 2+nchar(link1)),
                             subject=link1,
                             substitutionMatrix=submat1,
                             gapOpening = 0,
                             gapExtension = 1,
                             type="overlap")
  aln.1.df <- PairwiseAlignmentsSingleSubject2DF(aln.1)
  
  ## search after 1st part of linker for the 2nd part of linker
  aln.2 <- pairwiseAlignment(pattern=subseq(reads.l, max(pos.N)-1, nchar(linker)+1),
                             subject=link2,
                             substitutionMatrix=submat1,
                             gapOpening = 0,
                             gapExtension = 1,
                             type="overlap")
  aln.2.df <- PairwiseAlignmentsSingleSubject2DF(aln.2, max(pos.N)-2)
  
  goodIdx <- (aln.1.df$score >= nchar(link1)-maxMisMatch &
                aln.2.df$score >= nchar(link2)-maxMisMatch)
  
  primerID <- subseq(reads.l[goodIdx],
                     aln.1.df$end[goodIdx]+1,
                     aln.2.df$start[goodIdx]-1)
  
  reads.l <- subseq(reads.l[goodIdx], aln.2.df$end[goodIdx]+1)
  
  stopifnot(all(names(primerID)==names(reads.l)))
  
  return(list("reads.l"=reads.l,
              "primerID"=primerID))
}
##trim_primerIDlinker_side_reads(reads.l, linker)


#' subseqing, trim off reads from where marker start to match
#' when human part of sequence is short, ltr side read will read in to 
#' linker, and linker side reads may read into ltrbit, primer, etc
#' allow 1 mismatch for linker common
#' @param reads DNAStringSet of reads
#' @param marker over reading marker
#' @return DNAStringSet of reads with linker sequences removed
#' 
trim_overreading <- function(reads, marker, maxMisMatch=3) {
  
  stopifnot(class(reads) %in% "DNAStringSet")
  stopifnot(!any(duplicated(names(reads))))
  stopifnot(length(marker)==1)
  
  
  submat1 <- nucleotideSubstitutionMatrix(match=1,
                                          mismatch=0,
                                          baseOnly=TRUE)
  
  ## allows gap, and del/ins count as 1 mismatch
  tmp <- pairwiseAlignment(pattern=reads,
                           subject=marker,
                           substitutionMatrix=submat1,
                           gapOpening = 0,
                           gapExtension = 1,
                           type="overlap")
  
  odf <- PairwiseAlignmentsSingleSubject2DF(tmp)
  
  odf$isgood <- FALSE
  ## overlap in the middle or at right
  odf$isgood <- with(odf, ifelse(mismatch<=maxMisMatch &
                                   start>1,
                                 TRUE, isgood))
  
  ## overlap at left
  odf$isgood <- with(odf, ifelse(mismatch<=maxMisMatch &
                                   start==1 &
                                   width>=nchar(marker)-1,
                                 TRUE, isgood))
  
  ## note with ovelrap alignmment, it only align with a minimum of 1/2 of the shorter one
  odf$cut <- with(odf, ifelse(isgood, odf$start-1, nchar(reads)))
  if( any(odf$cut < nchar(reads)) ) {
    odf$cut <- nchar(reads)-nchar(marker)/2
    odf$cut <- with(odf, ifelse(isgood, odf$start-1, cut))
  }
  
  reads <- subseq(reads, 1, odf$cut)
}
##trim_overreading(reads.p, linker_common)
##trim_overreading(reads.l, largeLTRFrag)



getTrimmedSeqs <- function(qualityThreshold, badQuality, qualityWindow, primer,
                           ltrbit, largeLTRFrag, linker, linker_common, mingDNA,
                           read1, read2, alias, vectorSeq){
  
  ##### Load libraries #####
  ##library("hiReadsProcessor")
  ##library("ShortRead")
  
  stats <- data.frame()
  
  workingDir <- alias
  suppressWarnings(dir.create(workingDir, recursive=TRUE))
  setwd(workingDir)
  
  stats.bore <- data.frame(sample=alias)
  
  reads <- lapply(list(read1, read2), sapply, readFastq)
  
  stats.bore$barcoded <- sum(sapply(reads[[1]], length))
  
  r <- lapply(reads, function(x){
    seqs <- x[[1]]
    if(length(seqs) > 0){
      #remove anything after 5 bases under Q30 in 10bp window
      ##trimmed <- trimTailw(seqs, badQuality, qualityThreshold,
      ##                     round(qualityWindow/2))
      ## this step is not necessary at  all
      ## trim if 5 bases are below '0'(fred score 15) in a window of 10 bases
      ## trimmed <- trimTailw(seqs, 5, '+', 5)
      ## trimmed <- trimTailw(seqs, 5, '#', 5)
      ## this step is necessary because many shortreads functions work on ACGT only
      ##trimmed <- trimmed[width(trimmed) > 65]
      trimmed <- seqs
      trimmed <- trimmed[!grepl('N', sread(trimmed))]
      if(length(trimmed) > 0){
        trimmedSeqs <- sread(trimmed)
        trimmedqSeqs <- quality(quality(trimmed))
        names(trimmedSeqs) <- names(trimmedqSeqs) <- 
          sapply(sub("(.+) .+","\\1",ShortRead::id(trimmed)),
                 function(z){paste0(alias, "%", strsplit(z, "-")[[1]][2])})
      }
    }
    list(trimmedSeqs, trimmedqSeqs)
  })
  
  reads <- sapply(r, "[[", 1)
  qualities <- sapply(r, "[[", 2)
  #this is needed for primerID quality scores later on
  R1Quality <- qualities[[1]]
  rm(r)
  gc()
  
  ##stats.bore$p.qTrimmed <- length(reads[[2]])
  ##stats.bore$l.qTrimmed <- length(reads[[1]])
  print(t(stats.bore), quote=FALSE)
  
  # message("\nFilter and trim primer and ltrbit")
  ## .p suffix signifies the 'primer' side of the amplicon (i.e. read2)
  ## .l suffix indicates the 'liner' side of the amplicon (i.e. read1)
  reads.p <- trim_Ltr_side_reads(reads[[2]], primer, ltrbit)
  stats.bore$LTRed <- length(reads.p)
  
  # message("\nFilter and trim linker")
  readslprimer <- trim_primerIDlinker_side_reads(reads[[1]], linker)
  reads.l <- readslprimer$reads.l
  primerIDs <- readslprimer$primerID
  stats.bore$linkered <- length(reads.l)
  save(primerIDs, file="primerIDData.RData")
  
  ltrlinkeredQname <- intersect(names(reads.p), names(reads.l))
  reads.p <- reads.p[ltrlinkeredQname]
  reads.l <- reads.l[ltrlinkeredQname]
  stats.bore$ltredlinkered <- length(reads.l)
  
  print(t(stats.bore), quote=FALSE) 
  
  ## check if reads were sequenced all the way by checking for opposite adaptor
  # message("\nTrim reads.p over reading into linker")
  reads.p <- trim_overreading(reads.p, linker_common, 3)
  # message("\nTrim reads.l over reading into ltr")
  ## with mismatch=3, the 20 bases can not be found in human genome
  reads.l <- trim_overreading(reads.l, substr(largeLTRFrag, 1, 20), 3)
  
  # message("\nFilter on minimum length of ", mingDNA)
  reads.p <- subset(reads.p, width(reads.p) > mingDNA)
  reads.l <- subset(reads.l, width(reads.l) > mingDNA)
  
  ltrlinkeredQname <- intersect(names(reads.p), names(reads.l))
  reads.p <- reads.p[ltrlinkeredQname]
  reads.l <- reads.l[ltrlinkeredQname]
  stats.bore$lenTrim <- length(reads.p)
  
  # message("\nRemove reads align to vector") 
  vqName <- findVectorReads(file.path("..", vectorSeq),
                            paste0(primer, ltrbit),
                            reads.p, reads.l,
                            debug=TRUE)
  
  toload <- names(reads.p)[!names(reads.p) %in% vqName]
  reads.p <- reads.p[toload]
  reads.l <- reads.l[toload]
  stats.bore$vTrimed <- length(reads.p)
  
  ##dereplicate seqs for faster alignments
  ##this is re-expand at the beginning of callSeqs
  reads.p.u <- unique(reads.p)
  reads.l.u <- unique(reads.l)
  
  reads.p30.u <- unique(subseq(reads.p,1,mingDNA))
  
  stats.bore$uniqL <- length(reads.l.u)  
  stats.bore$uniqP <- length(reads.p.u)  
  stats.bore$uniqP30 <- length(reads.p30.u)  
  
  names(reads.p.u) <- seq_along(reads.p.u)
  names(reads.l.u) <- seq_along(reads.l.u)
  
  keys <- data.frame("R2"=match(reads.p, reads.p.u),
                     "R1"=match(reads.l, reads.l.u),
                     "names"=toload)
  keys$readPairKey <- paste0(keys$R1, "_", keys$R2)
  
  save(keys, file="keys.RData")
  
  stats.bore$lLen <- as.integer(mean(width(reads.l)))  
  stats.bore$pLen <- as.integer(mean(width(reads.p)))
  
  stats <- rbind(stats, stats.bore)
  
  save(stats, file="stats.RData")
  
  if(length(toload) > 0){
    
    chunks.p <- split(seq_along(reads.p.u), ceiling(seq_along(reads.p.u)/config$chunkSize))
    for(i in c(1:length(chunks.p))){
      writeXStringSet(reads.p.u[chunks.p[[i]]],
                      file=paste0("R2-", i, ".fa"),
                      append=FALSE)
    }
    
    chunks.l <- split(seq_along(reads.l.u), ceiling(seq_along(reads.l.u)/config$chunkSize))
    for(i in c(1:length(chunks.l))){    
      writeXStringSet(reads.l.u[chunks.l[[i]]],
                      file=paste0("R1-", i, ".fa"),
                      append=FALSE)
    }
    
    save(stats, file="stats.RData")
    alias #return 'value' which ultimately gets saved as trimStatus.RData
  }else{
    stop("error - no curated reads")
  }
}

processAlignments <- function(workingDir, minPercentIdentity, maxAlignStart, 
                              maxLength, refGenome){
  
  codeDir <- get(load("codeDir.RData"))
  source(paste0(codeDir, "/programFlow.R"))#for get_reference_genome function
  
  setwd(workingDir)
  
  #' clean up alignments and prepare for int site calling
  #'
  #' @param algns df with: 21 columns from BLAT (psl)
  #' @param from is "R1" or "R2" 
  #' @param refGenome character name of reference genome, ie. "hg38 
  #' @return Granges object
  
  processBLATData <- function(algns, from, refGenome){
    stopifnot(from == "R1" | from == "R2")
    algns$from <- from
    algns$qtStart <- ifelse(
      algns$strand == "+",
      (algns$tStart - (algns$qStart)),
      (algns$tStart - (algns$qSize - algns$qEnd - 1)))
    algns$qtEnd <- ifelse(
      algns$strand == "+",
      (algns$tEnd + (algns$qSize - algns$qEnd - 1)),
      (algns$tEnd + (algns$qStart)))    
    
    algns.gr <- GRanges(seqnames=Rle(algns$tName),
                        ranges = IRanges(
                          start = (algns$qtStart + 1), 
                          end = (algns$qtEnd)), #Convert to 1-base
                        strand=Rle(algns$strand),
                        seqinfo=seqinfo(get_reference_genome(refGenome)))
    
    mcols(algns.gr) <- algns[,c("from", "qName", "matches", "repMatches", 
                                "misMatches", "qStart", "qEnd", "qSize", 
                                "tBaseInsert")]
    rm(algns)
    algns.gr
  }
  
  load("keys.RData")
  load("stats.RData")
  
  psl.R2 <- list.files(".", pattern="R2.*.fa.psl.gz")
  message("R2 psl:\n", paste(psl.R2, collapse="\n"),"\n")
  hits.R2 <- readpsl(psl.R2)

  psl.R1 <- list.files(".", pattern="R1.*.fa.psl.gz")
  message("R1 psl:\n", paste(psl.R1, collapse="\n"),"\n")
  hits.R1 <- readpsl(psl.R1)
  
  #' Record the number of reads with R1 and R2 alignments
  readsAligning <- length(
    which(keys$R2 %in% hits.R2$qName & keys$R1 %in% hits.R1$qName))
  stats <- cbind(stats, readsAligning)
  save(stats, file="stats.RData")
  
  #' Quality filter and convert alignments from data.frame to GRanges
  hits.R2 <- qualityFilter(hits.R2, maxAlignStart, minPercentIdentity)
  hits.R2 <- processBLATData(hits.R2, "R2", refGenome)   
  save(hits.R2, file="hits.R2.RData")
  
  hits.R1 <- qualityFilter(hits.R1, maxAlignStart, minPercentIdentity)
  hits.R1 <- processBLATData(hits.R1, "R1", refGenome)
  save(hits.R1, file="hits.R1.RData")
  
  #no more '.p' or '.l' nomenclature here
  #we're combining alignments from both sides of the read

  #' All alignments should be either "+" or "-" strand.  
  stopifnot(all(strand(hits.R1) == "+" | strand(hits.R1) == "-"))
  stopifnot(all(strand(hits.R2) == "+" | strand(hits.R2) == "-"))
  
  #' Record the number of reads passing quality filtering
  readsWithGoodAlgnmts <- length(
    which(keys$R2 %in% hits.R2$qName & keys$R1 %in% hits.R1$qName))
  stats <- cbind(stats, readsWithGoodAlgnmts)
  save(stats, file="stats.RData")
  
  #' Identify all combinations of unique R1 and R2 sequences present in the data
  unique_key_pairs <- unique(keys[,c("R1", "R2", "readPairKey")])

  #' Reduced alignments identify the distinct genomic locations present in the 
  #' data for the R1 sequences (breakpoint positions) and R2 sequences 
  #' (integration site position). 
  #' Levels: Reads --> Unique Sequences --> Alignments --> Unique Genomic Loci
  red.hits.R1 <- reduce(
    flank(hits.R1, -1, start = TRUE), min.gapwidth = 0L, with.revmap = TRUE)

  red.hits.R2 <- reduce(
    flank(hits.R2, -1, start = TRUE), min.gapwidth = 0L, with.revmap = TRUE)
  
  #' The following finds all posible combinations of R1 and R2 loci which meet
  #' criteria for pairing. These include: oneEach (each pairing must come from
  #' one R1 and one R2 loci), opposite strands (paired loci should be present on
  #' opposite strands), and correct downstream orientation (if an R2 loci is on 
  #' the "+" strand, then the start of the R2 loci should be less than the 
  #' paired R1, and vice versa for "-" strand).
  #' (Inherent check for oneEach with findOverlaps())
  pairs <- findOverlaps(
    red.hits.R1,
    red.hits.R2,
    maxgap = maxLength,
    ignore.strand = TRUE
  )
  
  R1.loci <- red.hits.R1[queryHits(pairs)]
  R2.loci <- red.hits.R2[subjectHits(pairs)]
  
  #' Check isDownstream and isOppositeStrand
  R1.loci.starts <- start(R1.loci)
  R2.loci.starts <- start(R2.loci)
  
  R1.loci.strand <- strand(R1.loci)
  R2.loci.strand <- strand(R2.loci)
  
  keep.loci <- ifelse(
      R2.loci.strand == "+", 
      as.vector(R1.loci.starts > R2.loci.starts & 
                  R1.loci.strand != R2.loci.strand), 
      as.vector(R1.loci.starts < R2.loci.starts & 
                  R1.loci.strand != R2.loci.strand))
  
  keep.loci <- as.vector(
    keep.loci & R2.loci.strand != "*" & R1.loci.strand != "*")
  
  R1.loci <- R1.loci[keep.loci]
  R2.loci <- R2.loci[keep.loci]
  
  #' Below, the code constructs a genomic loci key which links genomic loci to
  #' the various R1 and R2 sequences that were aligned.
  loci.key <- data.frame(
    "R1.loci" = queryHits(pairs)[keep.loci],
    "R2.loci" = subjectHits(pairs)[keep.loci])
  loci.key$lociPairKey <- paste0(loci.key$R1.loci, ":", loci.key$R2.loci)
  
  loci.key$R1.qNames <- IntegerList(lapply(R1.loci$revmap, function(x){
    as.integer(hits.R1$qName[x])
  }))
  
  loci.key$R2.qNames <- IntegerList(lapply(R2.loci$revmap, function(x){
    as.integer(hits.R2$qName[x])
  }))
  
  loci.key$R1.readPairs <- IntegerList(lapply(
    loci.key$R1.qNames, function(x){
      which(unique_key_pairs$R1 %in% x)
  }))
  
  loci.key$R2.readPairs <- IntegerList(lapply(
    loci.key$R2.qNames, function(x){
      which(unique_key_pairs$R2 %in% x)
  }))
  
  #' Using the range information from the filtered paired alignments, the code
  #' constructs a GRanges object from the R1.loci and R2.loci. R2.loci are the 
  #' integration site positions while the R1.loci are the various breakpoints.
  #' The strand of the range is set to the same strand as the R2.loci since the
  #' direction of sequencing from the viral or vector genome is from the U5-host
  #' junction found at the 3' end of the integrated element.
  paired.loci <- GRanges(
    seqnames = seqnames(R2.loci), 
    ranges = IRanges(
      start = ifelse(strand(R2.loci) == "+", start(R2.loci), start(R1.loci)),
      end = ifelse(strand(R2.loci) == "+", end(R1.loci), end(R2.loci))),
    strand = strand(R2.loci),
    lociPairKey = loci.key$lociPairKey)
  
  paired.loci$readPairKeys <- CharacterList(lapply(
    1:length(paired.loci), 
    function(i){
      unique_key_pairs[intersect(
        loci.key$R1.readPairs[[i]], 
        loci.key$R2.readPairs[[i]]),
        "readPairKey"]
  }))
  
  #' Remove R1:R2 pairings that do not appear in the sequence data
  paired.loci <- paired.loci[sapply(paired.loci$readPairKeys, length) > 0]

  #' Expand readPairKeys and lociPairKeys to make a single object that maps loci
  #' to unique sequences. This is analogous to a sparse matrix, but in a 
  #' data.frame object. The keys object is still needed to jump from readPairKey
  #' to read name.
  read.loci.mat <- data.frame(
    "lociPairKey" = Rle(
      values = paired.loci$lociPairKey,
      lengths = sapply(paired.loci$readPairKeys, length)),
    "readPairKey" = unlist(paired.loci$readPairKeys)
  )
  
  #' Record the number of alignments that have been properly paired and passed
  #' filtering criteria.
  numProperlyPairedAlignments <- nrow(
    keys[keys$readPairKey %in% read.loci.mat$readPairKey,])
  stats <- cbind(stats, numProperlyPairedAlignments)
  save(stats, file="stats.RData")
  
  #' Templates aligning to single loci are termed unique, while templates
  #' aligning to multiple loci are termed multihits.
  readPairCounts <- table(read.loci.mat$readPairKey)
  uniq.readPairs <- names(readPairCounts[readPairCounts == 1])
  multihit.readPairs <- names(readPairCounts[readPairCounts > 1])

  #' ########## IDENTIFY IMPROPERLY-PAIRED READS (chimeras) ##########
  #' Further, all unique and multihit templates were mapped successfully to 
  #' genomic loci, yet some templates were sequenced but did not make it through
  #' the selection criteria. These template either do not have alignments to the
  #' reference genome (R1 or R2 did not align) or map to two distant genomic
  #' loci. The latter are termed chimeras and are considered to be artifacts of
  #' PCR amplification.
  failedReads <- keys[!keys$readPairKey %in% read.loci.mat$readPairKey,]
  chimera.reads <- failedReads[
    failedReads$R1 %in% hits.R1$qName & failedReads$R2 %in% hits.R2$qName,]
  
  chimera.alignments <- GRangesList(lapply(1:length(chimera.reads), function(i){
    R1 <- hits.R1[hits.R1$qName == chimera.reads[i, "R1"]]
    R2 <- hits.R2[hits.R2$qName == chimera.reads[i, "R2"]]
    names(R1) <- rep(chimera.reads[i, "names"], length(R1))
    names(R2) <- rep(chimera.reads[i, "names"], length(R2))
    c(R2, R1)
  }))
  
  chimeraData <- list(
    "read_info" = chimera.reads, "alignments" = chimera.alignments)
  save(chimeraData, file = "chimeraData.RData")
  
  #' Record chimera metrics
  chimeras <- length(unique(chimera.reads$names))
  stats <- cbind(stats, chimeras)
  save(stats, file="stats.RData")
  
  #' ########## IDENTIFY UNIQUELY-PAIRED READS (real sites) ##########
  #' Below, the paired.loci object is expanded to create the genomic alignments
  #' for each read that mapped to a single genomic loci. This data is then 
  #' recorded in two formats. "allSites" is a GRanges object where each row is a
  #' single read, while "sites.final" is a condensed form of the data where each
  #' row is a unique integration site with the width of the range refering to 
  #' the longest template aligned to the reference genome. 
  uniq.read.loci.mat <- read.loci.mat[
    read.loci.mat$readPairKey %in% uniq.readPairs,]
  
  uniq.templates <- paired.loci[
    match(uniq.read.loci.mat$lociPairKey, paired.loci$lociPairKey)]
  uniq.templates$readPairKeys <- NULL
  uniq.templates$readPairKey <- uniq.read.loci.mat$readPairKey
  
  uniq.keys <- keys[keys$readPairKey %in% uniq.readPairs,]
  uniq.reads <- uniq.templates[
    match(uniq.keys$readPairKey, uniq.templates$readPairKey)]
  names(uniq.reads) <- as.character(uniq.keys$names)
  uniq.reads$sampleName <- sapply(
    strsplit(as.character(uniq.keys$names), "%"), "[[", 1)
  uniq.reads$ID <- sapply(strsplit(as.character(uniq.keys$names), "%"), "[[", 2)
  
  allSites <- uniq.reads
  save(allSites, file="allSites.RData")

  sites.final <- dereplicateSites(allSites)
  if(length(sites.final)>0){
    sites.final$sampleName <- allSites[1]$sampleName
    sites.final$posid <- paste0(as.character(seqnames(sites.final)),
                                as.character(strand(sites.final)),
                                start(flank(sites.final, width=-1, start=TRUE)))
  }
  save(sites.final, file="sites.final.RData")
  
  #' Record metrics about unique alignments to the stats object
  numAllSingleReads <- length(allSites)
  stats <- cbind(stats, numAllSingleReads)
  numAllSingleSonicLengths <- length(unique(granges(allSites)))
  stats <- cbind(stats, numAllSingleSonicLengths)
  numUniqueSites <- length(sites.final)
  stats <- cbind(stats, numUniqueSites)

  #' Clean up environment for expansion and clustering of multihits
  rm(uniq.read.loci.mat, uniq.templates, uniq.keys, 
     uniq.reads, allSites, sites.final)
  gc()
  
  #' ########## IDENTIFY MULTIPLY-PAIRED READS (multihits) ##########
  #' Multihits are reads that align to multiple locations in the reference 
  #' genome. There are bound to always be a certain proportion of reads aligning
  #' to repeated sequence due to the high level degree of repeated DNA elements
  #' within genomes. The final object generated, "multihitData", is a list of 
  #' three objects. "unclusteredMultihits" is a GRanges object where every 
  #' alignment for every multihit read is present in rows. 
  #' "clusteredMultihitPositions" returns all the possible integration site 
  #' positions for the multihit. Lastly, "clusteredMultihitLengths" contains the
  #' length of the templates mapping to the multihit clusters, used for
  #' abundance calculations.
  unclusteredMultihits <- GRanges()
  clusteredMultihitPositions <- GRangesList()
  clusteredMultihitLengths <- list()
  
  if(length(multihit.readPairs) > 0){
    #' Only consider readPairKeys that aligned to multiple genomic loci
    multi.read.loci.mat <- read.loci.mat[
      read.loci.mat$readPairKey %in% multihit.readPairs,]
  
    multihit.templates <- paired.loci[
      paired.loci$lociPairKey %in% multi.read.loci.mat$lociPairKey]
    multihit.expansion.map <- multihit.templates$readPairKeys
    multihit.templates$readPairKeys <- NULL
    multihit.templates <- multihit.templates[Rle(
      values = 1:length(multihit.templates),
      lengths = sapply(multihit.expansion.map, length)
    )]
    multihit.templates$readPairKey <- unlist(multihit.expansion.map)
    
    #' As the loci are expanded from the paired.loci object, unique templates 
    #' and readPairKeys are present in the readPairKeys unlisted from the 
    #' paired.loci object.
    multihit.templates <- multihit.templates[
      multihit.templates$readPairKey %in% multi.read.loci.mat$readPairKey]
    
    multihit.keys <- keys[keys$readPairKey %in% multihit.readPairs,]
    multihit.keys$sampleName <- sapply(strsplit(
      as.character(multihit.keys$names), "%"), "[[", 1)
    multihit.keys$ID <- sapply(strsplit(
      as.character(multihit.keys$names), "%"), "[[", 2)

    #' Medians are based on all the potential sites for a given read, which will
    #' be identical for all reads associated with a readPairKey.
    multihit.medians <- round(
      median(width(split(multihit.templates, multihit.templates$readPairKey))))
    multihit.keys$medians <- multihit.medians[multihit.keys$readPairKey]
    
    multihits.pos <- flank(multihit.templates, -1, start = TRUE)
    multihits.red <- reduce(multihits.pos, min.gapwidth = 5L, with.revmap = TRUE)  #! Should make 5L a option
    revmap <- multihits.red$revmap
    
    axil_nodes <- as.character(Rle(
      values = multihit.templates$readPairKey[sapply(revmap, "[", 1)], 
      lengths = sapply(revmap, length)
    ))
    nodes <- multihit.templates$readPairKey[unlist(revmap)]
    edgelist <- unique(matrix( c(axil_nodes, nodes), ncol = 2 ))
    multihits.clusterData <- igraph::clusters(
      igraph::graph.edgelist(edgelist, directed=F))
    clus.key <- data.frame(
      row.names = unique(as.character(t(edgelist))),
      "clusID" = multihits.clusterData$membership)
    
    multihits.pos$clusID <- clus.key[multihits.pos$readPairKey, "clusID"]
    clusteredMultihitPositions <- split(multihits.pos, multihits.pos$clusID)
    clusteredMultihitNames <- lapply(
      clusteredMultihitPositions, function(x) unique(x$readPairKey))
    clusteredMultihitPositions <- GRangesList(lapply(
      clusteredMultihitPositions, 
      function(x){
        unname(unique(granges(x)))
    })) 
    
    clusteredMultihitLengths <- lapply(clusteredMultihitNames, function(x){
      readIDs <- unique(multihit.keys[multihit.keys$readPairKey %in% x,]$ID)
      data.frame(table(multihit.keys[multihit.keys$ID %in% readIDs,]$medians))
    })
    
    #' Expand the multihit.templates object from readPairKey specific to read
    #' specific.
    multihit.keys <- multihit.keys[order(multihit.keys$readPairKey),]
    multihit.readPair.read.exp <- IntegerList(lapply(
      unique(multihit.keys$readPairKey), 
      function(x){which(multihit.keys$readPairKey == x)}))
    names(multihit.readPair.read.exp) <- unique(multihit.keys$readPairKey)
    unclusteredMultihits <- multihit.templates
    multihit.readPair.read.exp <- as(multihit.readPair.read.exp[
      as.character(unclusteredMultihits$readPairKey)], "SimpleList")
    unclusteredMultihits <- unclusteredMultihits[Rle(
      values = 1:length(unclusteredMultihits),
      lengths = sapply(multihit.readPair.read.exp, length)
    )]
    names(unclusteredMultihits) <- multihit.keys$names[
      unlist(multihit.readPair.read.exp)]
    unclusteredMultihits$ID <- multihit.keys$ID[
      unlist(multihit.readPair.read.exp)]
    unclusteredMultihits$sampleName <- multihit.keys$sampleName[
      unlist(multihit.readPair.read.exp)]
  }
  
  stopifnot(length(clusteredMultihitPositions)==length(clusteredMultihitLengths))
  multihitData <- list(unclusteredMultihits, clusteredMultihitPositions, clusteredMultihitLengths)
  names(multihitData) <- c("unclusteredMultihits", "clusteredMultihitPositions", "clusteredMultihitLengths")
  
  save(multihitData, file="multihitData.RData")
  
  #' Record multihit metrics (reads, clusters, sonicLengths)
  multihitReads <- nrow(keys[keys$readPairKey %in% multihit.readPairs,])
  stats <- cbind(stats, multihitReads)

  multihitSonicLengths <- 0
  if( length(multihitData$clusteredMultihitLengths) > 0 ) {
    multihitSonicLengths <- sum(
      sapply(multihitData$clusteredMultihitLengths, nrow))
  }
  stats <- cbind(stats, multihitSonicLengths) 
  
  multihitClusters <- length(multihitData$clusteredMultihitPositions) 
  stats <- cbind(stats, multihitClusters)
  
  #' Finalize metrics by combining unique alignments and multihit clusters
  totalSonicLengths <- numAllSingleSonicLengths + multihitSonicLengths
  stats <- cbind(stats, totalSonicLengths)

  totalEvents <- numUniqueSites + multihitClusters
  stats <- cbind(stats, totalEvents)
  
  save(stats, file="stats.RData")
  
  ####### END OF PROCESS ALIGNMENTS ########
}