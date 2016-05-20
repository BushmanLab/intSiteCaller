libs <- c("stringr",
          "ShortRead",
          "BSgenome")
null <- suppressMessages(sapply(libs, library, character.only=TRUE))


#' note for openlava version of LSF, done does not work. Use ended instead.
#' 
bsub <- function(queue="normal", cpus=1, maxmem=NULL, wait=NULL, jobName=NULL, logFile=NULL, command=NULL){
    stopifnot(!is.null(maxmem))
    stopifnot(!is.null(command))
    
    cmd <- paste0("bsub -q ", queue, " -n ", as.character(cpus), " -M ", maxmem)
    ##cmd <- sprintf("bsub -q %s -n %s -M %s", queue, cpus, maxmem)
    
    if(!is.null(wait)){
        LSF.VERSION <- system2("bsub", "-V", stdout=TRUE, stderr=TRUE)[1]
        if( grepl("openlava", LSF.VERSION, ignore.case=TRUE) ) {
            wait <- sub("done", "ended", wait)
        }
        cmd <- paste0(cmd, " -w \"", wait, "\"")
    }
    
    if(!is.null(jobName)) cmd <- paste0(cmd, " -J \"", jobName, "\"")
    if(!is.null(logFile)) cmd <- paste0(cmd, " -o ", logFile)
    
    cmd <- paste0(cmd, " ", command)
    message(cmd)
    system(cmd)
}

#takes a textual genome identifier (ie. hg18) and turns it into the correct
#BSgenome object
get_reference_genome <- function(reference_genome) {
  pattern <- paste0("\\.", reference_genome, "$")
  match_index <- which(grepl(pattern, installed.genomes()))
  if (length(match_index) != 1) {
    write("Installed genomes are:", stderr())
    write(installed.genomes(), stderr())
    stop(paste("Cannot find unique genome for", reference_genome))
  }
  BS_genome_full_name <- installed.genomes()[match_index]
  library(BS_genome_full_name, character.only=T)
  get(BS_genome_full_name)
}

#' align sequences
#' Note: it is important not to change the blat parameters.
#' The parameters were optimized after lengthy experimentations.
#' Leave them as they are unless there is a specific reason other than
#' curoisity. Hard coded for a reason.
#' 
#' To try different blat parameters, create a file named blatOverzRide.txt
#' in the root analysis folder with the blat command template such as
#' 
#' [@node063 I1]$ cat blatOverRide.txt
#' blat %s.2bit %s %s.psl -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead
#' [@node063 I1]$
#' 
alignSeqs <- function(){
    
    Sys.sleep(1)
    sampleID <- as.integer(Sys.getenv("LSB_JOBINDEX"))
    message("LSB_JOBINDEX=", sampleID)
  
    toAlign <- get(load("toAlign.RData"))
    alignFile <- toAlign[sampleID]
    
    message("alignFile=", alignFile)
    alias <- strsplit(alignFile, "/")[[1]][1]
    
    completeMetadata <- get(load("completeMetadata.RData"))
    genome <- completeMetadata[completeMetadata$alias==alias,"refGenome"]
    indexPath <- paste0(genome, ".2bit")
    
    blatTemplate <- "blat %s.2bit %s %s.psl -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead"
    if( file.exists("blatOverRide.txt") ) {
        blatTemplate <- readLines("blatOverRide.txt")
        message("Blat parameters were overridden by file blatOverRide.txt")
    }
    cmd <-sprintf(blatTemplate, genome, alignFile, alignFile)
    message(cmd)
    unlink(paste0(alignFile, c(".psl", ".psl.gz")), force=TRUE)
    system(cmd)
    
    system(paste0("gzip ", alignFile, ".psl"))
}

callIntSites <- function(){
  Sys.sleep(1)
  message("LSB_JOBINDEX=", Sys.getenv("LSB_JOBINDEX"))
  
  codeDir <- get(load("codeDir.RData"))
  source(file.path(codeDir, "intSiteLogic.R"))
  
  ##sampleID <- as.integer(system("echo $LSB_JOBINDEX", intern=T))
  sampleID <- as.integer(Sys.getenv("LSB_JOBINDEX"))
  message(sampleID)
  
  completeMetadata <- get(load("completeMetadata.RData"))[sampleID,]
  print(t(completeMetadata), quote=FALSE)  
  
  status <- tryCatch(eval(as.call(append(processAlignments,
                                         unname(as.list(completeMetadata[c("alias", "minPctIdent",
                                                                           "maxAlignStart", "maxFragLength",
                                                                           "refGenome", "uniqRegion")]))))),
                     error=function(e){print(paste0("Caught error: ", e$message))})

  save(status, file="callStatus.RData") #working directory is changed while executing getTrimmedSeqs
}

demultiplex <- function(){
  I1 <- readFasta(list.files("Data", pattern="correctedI1-.", full.names=T))
  
  completeMetadata <- get(load("completeMetadata.RData"))
  
  I1 <- I1[as.vector(sread(I1)) %in% completeMetadata$bcSeq]
  samples <- completeMetadata[match(as.character(sread(I1)), completeMetadata$bcSeq), "alias"]
  
  #only necessary if using native data - can parse out description w/ python
  I1Names <-  sapply(strsplit(as.character(ShortRead::id(I1)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!
  
  unlink("Data/demultiplexedReps", recursive=TRUE,  force=TRUE)
  suppressWarnings(dir.create("Data/demultiplexedReps"))
  
  R1 <- readFastq("Data/Undetermined_S0_L001_R1_001.fastq.gz")
  demultiplex_reads(R1, "R1", I1Names, samples, completeMetadata)
  
  R2 <- readFastq("Data/Undetermined_S0_L001_R2_001.fastq.gz")
  demultiplex_reads(R2, "R2", I1Names, samples, completeMetadata)
}

#' write fastq for each barcode and each sample
#' @param reads fastq reads as parsed by readFastq()
#' @param suffix either "R1" or "R2"
demultiplex_reads <- function(reads, suffix, I1Names, samples, completeMetadata) {
    RNames <- sapply(strsplit(as.character(ShortRead::id(reads)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!
    names(RNames) <- NULL

    reads <- reads[match(I1Names, RNames)]
    reads <- split(reads, samples)
    for (i in 1:length(reads)){
        barcode.i <- completeMetadata$bcSeq[ completeMetadata$alias == names(reads)[i] ]
        stopifnot(length(barcode.i)==1)
        alias_by_barcode <- completeMetadata$alias[ completeMetadata$bcSeq == barcode.i ]
        stopifnot(length(alias_by_barcode)>=1)
        fqFiles <- paste0("Data/demultiplexedReps/", alias_by_barcode, "_", suffix, ".fastq.gz")
        cat(barcode.i, "\t", paste(fqFiles, collapse=" "), "\n" )
        null <- sapply(fqFiles, function(fq) writeFastq(reads[[i]], fq, mode="w") )
    }  
}

#' note bsub job dependency is better to use done which wait for job exit 0;
#' here ended is used because openlava cannot handle done correctly.
errorCorrectBC <- function(){
  library("ShortRead")
  
  codeDir <- get(load("codeDir.RData"))
  completeMetadata <- get(load("completeMetadata.RData"))
  bushmanJobID <- get(load("bushmanJobID.RData"))
  
  I1 <- readFastq("Data/Undetermined_S0_L001_I1_001.fastq.gz")
  I1 <- trimTailw(I1, 2, "0", 12)
  I1 <- I1[width(I1)==max(width(I1))]
  
  I1 <- split(I1, ceiling(seq_along(I1)/500000))
  
  for(chunk in names(I1)){
    writeFasta(I1[[chunk]], file=paste0("Data/trimmedI1-", chunk, ".fasta"))
  }
    
  bsub(jobName=sprintf("BushmanErrorCorrectWorker_%s[1-%s]", bushmanJobID, length(I1)),
       maxmem=1000,
       logFile="logs/errorCorrectWorkerOutput%I.txt",
       command=paste0("python ", codeDir, "/errorCorrectIndices/processGolay.py")
  )
  
  bsub(wait=sprintf("ended(BushmanErrorCorrectWorker_%s)", bushmanJobID),
       jobName=sprintf("BushmanDemultiplex_%s", bushmanJobID),
       maxmem=64000, #just in case
       logFile="logs/demultiplexOutput.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); demultiplex();\"")
  )
  
  #trim seqs
  ##bsub(wait=paste0("done(BushmanDemultiplex_", bushmanJobID, ")"),
  bsub(wait=sprintf("ended(BushmanDemultiplex_%s)", bushmanJobID),
       jobName=sprintf("BushmanTrimReads_%s[1-%s]", bushmanJobID, nrow(completeMetadata)),
       maxmem=16000,
       logFile="logs/trimOutput%I.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); trimReads();\"")
  )
  
  ##post-trim processing, also kicks off alignment and int site calling jobs
  bsub(wait=sprintf("ended(BushmanTrimReads_%s)", bushmanJobID),
       jobName=sprintf("BushmanPostTrimProcessing_%s", bushmanJobID),
       maxmem=8000,
       logFile="logs/postTrimOutput.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); postTrimReads();\"")
  )
}


postTrimReads <- function(){

  Sys.sleep(1)
  library("BSgenome")
  library("rtracklayer") #needed for exporting genome to 2bit
  completeMetadata <- get(load("completeMetadata.RData"))
  codeDir <- get(load("codeDir.RData"))
  bushmanJobID <- get(load("bushmanJobID.RData"))
  
  numAliases <- nrow(completeMetadata)
  
  toAlign <- list.files(".", "R[12]-.*fa$", recursive=TRUE)
  toAlign <- toAlign[order(-file.info(toAlign)$size)]
  save(toAlign, file="toAlign.RData", compress=FALSE)
  numFastaFiles <- length(toAlign)


  
  #make temp genomes
  genomesToMake <- unique(completeMetadata$refGenome)
  
  for(genome in genomesToMake){
    export(get_reference_genome(genome), paste0(genome, ".2bit"))
    ##system(paste0("blat ", genome, ".2bit /dev/null /dev/null -makeOoc=", genome, ".11.ooc"))
  }
    
  #align seqs
  bsub(wait=sprintf("ended(BushmanPostTrimProcessing_%s)", bushmanJobID),
       jobName=sprintf("BushmanAlignSeqs_%s[1-%s]", bushmanJobID, numFastaFiles),
       maxmem=12000,
       logFile="logs/alignOutput%I.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); alignSeqs();\"")
  )
  
  #call int sites (have to find out which ones worked)
  bsub(wait=sprintf("ended(BushmanAlignSeqs_%s)", bushmanJobID),
       jobName=sprintf("BushmanCallIntSites_%s[1-%s]", bushmanJobID, nrow(completeMetadata)),
       maxmem=120000, #multihits suck lots of memory
       logFile="logs/callSitesOutput%I.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); callIntSites();\"")
  )
  
  
  bsub(wait=sprintf("ended(BushmanCallIntSites_%s)", bushmanJobID),
       jobName=sprintf("BushmanErrorCheck_%s", bushmanJobID),
       maxmem=4000,
       logFile="logs/errorCheck.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); check_error();\"")
       )
  
  
}

trimReads <- function(){
    
  Sys.sleep(1)
  
  codeDir <- get(load("codeDir.RData"))
  source(file.path(codeDir, "intSiteLogic.R"))
  
  ##sampleID <- as.integer(system("echo $LSB_JOBINDEX", intern=T))
  sampleID <- as.integer(Sys.getenv("LSB_JOBINDEX"))
  message("$LSB_JOBINDEX=",sampleID)
  
  completeMetadata <- get(load("completeMetadata.RData"))[sampleID,]
  
  alias <- completeMetadata$alias
  print(t(as.data.frame(completeMetadata)), quote=FALSE)
  
  suppressWarnings(dir.create(alias, recursive=TRUE))
  
  status <- tryCatch(eval(as.call(append(getTrimmedSeqs,
                                         unname(as.list(completeMetadata[c("qualityThreshold", "badQualityBases",
                                                                           "qualitySlidingWindow", "primer", "ltrBit",
                                                                           "largeLTRFrag", "linkerSequence", "linkerCommon",
                                                                           "mingDNA", "read1", "read2", "alias", "vectorSeq")]))))),
                     error=function(e){print(paste0("Caught error: ", e$message))})
  
  save(status, file="trimStatus.RData") #working directory is changed while executing getTrimmedSeqs
}

processMetadata <- function(){

  bushmanJobID <- parsedArgs$jobID
  
  #stop if any jobs already exist with the same job ID as this will confuse LSF
  stopifnot(!any(grepl(bushmanJobID, suppressWarnings(system("bjobs -l | grep -o \"Job Name <[^>]*>\"", intern=T)))))
  
  #expand codeDir to absolute path for saving
  codeDir <- normalizePath(parsedArgs$codeDir)

  source(file.path(codeDir, 'linker_common.R'))
  source(file.path(codeDir, 'read_sample_files.R'))

  #setting R's working dir also sets shell location for system calls, thus
  #primaryAnalysisDir is propagated without being saved
  setwd(parsedArgs$primaryAnalysisDir)

  save(bushmanJobID, file=paste0(getwd(), "/bushmanJobID.RData"))
  save(codeDir, file=paste0(getwd(), "/codeDir.RData"))

  sample_file <- 'sampleInfo.tsv'
  proc_file <- "processingParams.tsv"
  if ( ! file.exists(proc_file)) { # have to use default
      default <- "default_processingParams.tsv"
      proc_file <- file.path(codeDir, default)
  }
  completeMetadata <- read_sample_processing_files(sample_file, proc_file)
  
  completeMetadata$read1 <- paste0(getwd(), "/Data/demultiplexedReps/", completeMetadata$alias, "_R1.fastq.gz")
  completeMetadata$read2 <- paste0(getwd(), "/Data/demultiplexedReps/", completeMetadata$alias, "_R2.fastq.gz")

  stopifnot(all(c("qualityThreshold", "badQualityBases", "qualitySlidingWindow",
                  "primer", "ltrBit", "uniqRegion", "largeLTRFrag", "linkerSequence", "linkerCommon",
                  "mingDNA", "read1", "read2", "alias", "vectorSeq", "minPctIdent",
                  "maxAlignStart", "maxFragLength", "gender") %in% names(completeMetadata)))
  
  stopifnot(all( file.exists(completeMetadata$vectorSeq) ))
  
  
  ## check primer, ltrBit, largeLTRFrag consistency
  ## largeLTRFrag should start with RC(primer+ltrBit)
  rc.primer <- as.character(
      reverseComplement(DNAStringSet(completeMetadata$primer)))
  rc.ltrbit <- as.character(
      reverseComplement(DNAStringSet(completeMetadata$ltrBit)))
  
  rc.primerltrbitInlargeLTR <- mapply(function(x,y, z) grepl(y, x) | grepl(z, x),
                                      x=completeMetadata$largeLTRFrag,
                                      y=rc.primer,
                                      z=rc.ltrbit)
  
  if(!all(rc.primerltrbitInlargeLTR)) {
      print(data.frame(PLTR=names(rc.primerltrbitInlargeLTR),
                       rc.primerltrbitInlargeLTR))
      stop()
  }
  
  
  save(completeMetadata, file="completeMetadata.RData")

  suppressWarnings(dir.create("logs"))

  #error-correct barcodes - kicks off subsequent steps
  bsub(jobName=paste0("BushmanErrorCorrect_", bushmanJobID),
       maxmem=20000,
       logFile="logs/errorCorrectOutput.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); errorCorrectBC();\"")
  )
}


check_error <- function(errFile="error.txt") {
    message("Errors if any were written to file ", errFile)
    cmd <- "grep -i \"exit\\|halt\\|huge\" logs/*.txt"
    err <- system(cmd, intern=TRUE)
    cmd <-  "grep -i max logs/*.txt | grep -i memory | awk '{print $1, $(NF-1)}' | sort -k2nr"
    mem <- system(cmd, intern=TRUE)
    if (length(err)==0) err <- "No obvious error found"
    write(c(err, "\nMemory usage", mem), errFile)
}

