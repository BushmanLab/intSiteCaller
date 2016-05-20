## needs blat and python
commandLinePrograms <- c("blat", "python")
programsPresent <- sapply(commandLinePrograms, function(app) system2("which", app, stderr=NULL, stdout=NULL))==0
if(any(!programsPresent)){
  stop(paste(commandLinePrograms[!programsPresent]), " is not available")
}

## R packages
rPackages <- c("ShortRead", "GenomicRanges",
               "rtracklayer", "BSgenome", "argparse", "igraph")
rPackagesPresent <- is.element(rPackages, installed.packages()[,1])
if(any(!rPackagesPresent)){
    stop(paste(rPackages[!rPackagesPresent]), " is not available")
}

codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))

library("argparse", quietly=T)

## get args
parser <- ArgumentParser(formatter_class='argparse.RawTextHelpFormatter')

parser$add_argument("-j", "--jobID", type="character", nargs=1,
                    default="intSiteCallerJob",
                    help="Unique name by which to identify this intance of intSiteCaller [default: %(default)s]")
parser$add_argument("-c", "--codeDir", type="character", nargs=1,
                    default=codeDir,
                    help="Directory where intSiteCaller code is stored, can be relative or absolute [default: %(default)s]")
parser$add_argument("-p", "--primaryAnalysisDir", type="character",
                    default=".",
                    help="Location of primary analysis directory, can be relative or absolute [default: %(default)s]")

parsedArgs <- parser$parse_args(commandArgs(trailingOnly = TRUE))

if( !any(commandArgs(trailingOnly = TRUE) %in% c("-j", "--jobID")) ) {
    parsedArgs$jobID <- basename(normalizePath(parsedArgs$primaryAnalysisDir))
}

source(file.path(parsedArgs$codeDir, "programFlow.R")) 

processMetadata()