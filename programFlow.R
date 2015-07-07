
bsub <- function(cpus=1, maxmem=NULL, wait=NULL, jobName=NULL, logFile=NULL, command=NULL){
  stopifnot(!is.null(maxmem))
  if(Sys.Date()>=as.Date("2015-06-01", "%Y-%m-%d")){
    queue <- "normal"
  }else{
    queue <- "umem"
  }

  cmd <- paste0("bsub -q ", queue, " -n ", as.character(cpus), " -M ", maxmem)
  
  if(!is.null(wait)){
    cmd <- paste0(cmd, " -w \"", wait, "\"")
  }
  
  if(!is.null(jobName)){
    cmd <- paste0(cmd, " -J \"", jobName, "\"")
  }
  
  if(!is.null(logFile)){
    cmd <- paste0(cmd, " -o ", logFile)
  }
  
  cmd <- paste0(cmd, " ", command) #no default, should crash if no command provided
  system(cmd)
}

#takes a textual genome identifier (ie. hg18) and turns it into the correct
#BSgenome object
get_reference_genome <- function(reference_genome) {
  library("BSgenome")
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

alignSeqs <- function(){
  # do alignments  
  toAlign <- system("ls */*.fa", intern=T)
  alignFile <- toAlign[as.integer(system("echo $LSB_JOBINDEX", intern=T))]
  alias <- strsplit(alignFile, "/")[[1]][1]
  
  completeMetadata <- get(load("completeMetadata.RData"))
  genome <- completeMetadata[completeMetadata$alias==alias,"refGenome"]
  indexPath <- paste0(genome, ".2bit")
  
  system(paste0("blat ", indexPath, " ", alignFile, " -ooc=", genome, ".11.ooc ", alignFile, ".psl -tileSize=11 -repMatch=112312 -t=dna -q=dna -minIdentity=85 -minScore=27 -dots=1000 -out=psl -noHead"))
  system(paste0("gzip ", alignFile, ".psl"))
}

callIntSites <- function(){
  codeDir <- get(load("codeDir.RData"))
  source(paste0(codeDir, "/intSiteLogic.R"))
  
  sampleID <- as.integer(system("echo $LSB_JOBINDEX", intern=T))
  
  completeMetadata <- get(load("completeMetadata.RData"))[sampleID,]

  status <- tryCatch(eval(as.call(append(processAlignments,
                                         unname(as.list(completeMetadata[c("alias", "minPctIdent",
                                                                           "maxAlignStart", "maxFragLength",
                                                                           "refGenome")]))))),
                     error=function(e){print(paste0("Caught error: ", e$message))})

  save(status, file="callStatus.RData") #working directory is changed while executing getTrimmedSeqs
}

cleanup <- function(){
  cleanup <- get(load("cleanup.RData"))
  
  if(cleanup){
    system("rm *.2bit", ignore.stderr=T)
    system("rm *.ooc", ignore.stderr=T)
    system("rm *.RData", ignore.stderr=T)
    system("rm Data/*.fasta", ignore.stderr=T)
    system("rm */hits.R*.RData", ignore.stderr=T)
    system("rm */R*.fa*", ignore.stderr=T)
    system("rm */keys.RData", ignore.stderr=T)
    system("rm */*Status.RData", ignore.stderr=T)
    system("rm -r logs", ignore.stderr=T)
    system("rm -r Data/demultiplexedReps", ignore.stderr=T)
  }
}

demultiplex <- function(){
  #Demultiplexing is currently a single-core process - perhaps it could be made
  #more efficient by having each error-correct worker do its own
  #mini-demultiplex with the barcodes that it error corrected, then write their 
  #own shorter fastq files which can be cat'd together after everything is done
  
  library("ShortRead")
  
  I1 <- readFasta(list.files("Data", pattern="correctedI1-.", full.names=T)) #readFasta("Data/correctedI1.fasta")
  
  completeMetadata <- get(load("completeMetadata.RData"))
  
  I1 <- I1[as.vector(sread(I1)) %in% completeMetadata$bcSeq]
  samples <- completeMetadata[match(as.character(sread(I1)), completeMetadata$bcSeq), "alias"]
  
  #only necessary if using native data - can parse out description w/ python
  I1Names <-  sapply(strsplit(as.character(ShortRead::id(I1)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!

  rm(I1)
  
  suppressWarnings(dir.create("Data/demultiplexedReps"))
  
  #R1 
  R1 <- readFastq("Data/Undetermined_S0_L001_R1_001.fastq.gz")
  R1Names <- sapply(strsplit(as.character(ShortRead::id(R1)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!
  names(R1Names) <- NULL
  
  R1 <- R1[match(I1Names, R1Names)]
  R1 <- split(R1, samples)
  for (i in 1:length(R1)){writeFastq(R1[[i]], paste0("Data/demultiplexedReps/", names(R1[i]), "_R1.fastq.gz"), mode="w")}
  
  rm(R1, R1Names)
  
  #R2
  R2 <- readFastq("Data/Undetermined_S0_L001_R2_001.fastq.gz")
  R2Names <- sapply(strsplit(as.character(ShortRead::id(R2)), " "), "[[", 1) #for some reason we can't dynamically set name/id on ShortRead!
  names(R2Names) <- NULL
  
  R2 <- R2[match(I1Names, R2Names)]
  R2 <- split(R2, samples)
  
  for (i in 1:length(R2)){writeFastq(R2[[i]], paste0("Data/demultiplexedReps/", names(R2[i]), "_R2.fastq.gz"), mode="w")}
  
  rm(R2, R2Names, I1Names, samples) 
}

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
    
  bsub(jobName=paste0("BushmanErrorCorrectWorker_", bushmanJobID, "[1-", length(I1),"]"),
       maxmem=1000,
       logFile="logs/errorCorrectWorkerOutput%I.txt",
       command=paste0("python ", codeDir, "/errorCorrectIndices/processGolay.py")
  )
  
  bsub(wait=paste0("done(BushmanErrorCorrectWorker_", bushmanJobID, ")"),
       jobName=paste0("BushmanDemultiplex_", bushmanJobID),
       maxmem=64000, #just in case
       logFile="logs/demultiplexOutput.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); demultiplex();\"")
  )
  
  #trim seqs
  bsub(wait=paste0("done(BushmanDemultiplex_", bushmanJobID, ")"),
       jobName=paste0("BushmanTrimReads_", bushmanJobID, "[1-", nrow(completeMetadata), "]"),
       maxmem=16000,
       logFile="logs/trimOutput%I.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); trimReads();\"")
  )
  
  #post-trim processing, also kicks off alignment and int site calling jobs
  bsub(wait=paste0("done(BushmanTrimReads_", bushmanJobID, ")"),
       jobName=paste0("BushmanPostTrimProcessing_", bushmanJobID),
       maxmem=8000,
       logFile="logs/postTrimOutput.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); postTrimReads();\"")
  )
}

postTrimReads <- function(){
  library("BSgenome")
  library("rtracklayer") #needed for exporting genome to 2bit
  completeMetadata <- get(load("completeMetadata.RData"))
  codeDir <- get(load("codeDir.RData"))
  bushmanJobID <- get(load("bushmanJobID.RData"))
  
  numAliases <- nrow(completeMetadata)
  
  numFastaFiles <- length(system("ls */*.fa", intern=T))
  
  #make temp genomes
  genomesToMake <- unique(completeMetadata$refGenome)
  
  for(genome in genomesToMake){
    export(get_reference_genome(genome), paste0(genome, ".2bit"))
    system(paste0("blat ", genome, ".2bit /dev/null /dev/null -makeOoc=", genome, ".11.ooc"))
  }
    
  #align seqs
  bsub(wait=paste0("done(BushmanPostTrimProcessing_", bushmanJobID, ")"),
       jobName=paste0("BushmanAlignSeqs_", bushmanJobID, "[1", "-", numFastaFiles, "]"),
       maxmem=8000,
       logFile="logs/alignOutput%I.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); alignSeqs();\"")
  )
  
  #call int sites (have to find out which ones worked)
  successfulTrims <- unname(sapply(completeMetadata$alias, function(x){
    get(load(paste0(x, "/trimStatus.RData"))) == x    
  }))
  
  bsub(wait=paste0("done(BushmanAlignSeqs_", bushmanJobID, ")"),
       jobName=paste0("BushmanCallIntSites_", bushmanJobID, "[", paste(which(successfulTrims), collapse=","), "]"),
       maxmem=24000, #multihits suck lots of memory
       logFile="logs/callSitesOutput%I.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); callIntSites();\"")
  )
  
  bsub(wait=paste0("done(BushmanCallIntSites_", bushmanJobID, ")"),
       jobName=paste0("BushmanCleanup_", bushmanJobID),
       maxmem=1000,
       logFile="logs/cleanupOutput.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); cleanup();\"")
  )
}

trimReads <- function(){
  codeDir <- get(load("codeDir.RData"))
  source(paste0(codeDir, "/intSiteLogic.R"))
  
  sampleID <- as.integer(system("echo $LSB_JOBINDEX", intern=T))
  
  completeMetadata <- get(load("completeMetadata.RData"))[sampleID,]
    
  alias <- completeMetadata$alias
  
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
  cleanup <- parsedArgs$cleanup

    source(file.path(codeDir, 'linker_common.R'))
    source(file.path(codeDir, 'read_sample_files.R'))

  #setting R's working dir also sets shell location for system calls, thus
  #primaryAnalysisDir is propagated without being saved
  setwd(parsedArgs$primaryAnalysisDir)

  save(bushmanJobID, file=paste0(getwd(), "/bushmanJobID.RData"))
  save(codeDir, file=paste0(getwd(), "/codeDir.RData"))
  save(cleanup, file=paste0(getwd(), "/cleanup.RData"))

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
                  "primer", "ltrBit", "largeLTRFrag", "linkerSequence", "linkerCommon",
                  "mingDNA", "read1", "read2", "alias", "vectorSeq", "minPctIdent",
                  "maxAlignStart", "maxFragLength", "gender") %in% names(completeMetadata)))
  
  stopifnot(all( file.exists(completeMetadata$vectorSeq) ))
  
  save(completeMetadata, file="completeMetadata.RData")

  suppressWarnings(dir.create("logs"))

  #error-correct barcodes - kicks off subsequent steps
  bsub(jobName=paste0("BushmanErrorCorrect_", bushmanJobID),
       maxmem=6000,
       logFile="logs/errorCorrectOutput.txt",
       command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); errorCorrectBC();\"")
  )
}
