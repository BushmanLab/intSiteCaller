library(stringr, quietly = TRUE)
library(plyr, quietly = TRUE)
options(stringsAsFactors = FALSE)

### make sure the folder is clean ####
if( length(list.files(".", pattern="*.RData", recursive=TRUE))>0 |
   length(list.files(".", pattern="Data/*.fasta", recursive=TRUE))>0 ) stop(
   "Previous run trace left, first run\ngit clean -df\n" )


### run intSiteValidation test data and wait until finish ####
testrunlog <- "testrun.log"
cmd <- sprintf('Rscript ../../intSiteCaller.R > %s 2>&1', testrunlog)
message(cmd)
if( system(cmd)!=0 ) stop(cmd, " not executed")

jobid_lines <- grep("Job.*<\\d+>", readLines(testrunlog), value=TRUE)
stopifnot( length(jobid_lines)==1 )
message(jobid_lines)
message("This test should finish in 10 minutes if the queues are not busy.")

start_jobid <- as.integer(stringr:::str_match(jobid_lines, "<(\\d+.)>")[2])
stopifnot( length(start_jobid)==1 )

still_running <- TRUE
minutes <- 0
while( any(still_running) ) {
    mybjobsid <- as.integer(system("bjobs | grep -v JOBID | awk '{print $1}'", intern=TRUE))
    still_running <- mybjobsid>=start_jobid
    message("Running_minutes: ", minutes, "\tjobs: ", sum(still_running))
    if ( any(still_running) ) {
        Sys.sleep(60)
        minutes <- minutes + 1
    }
}
message("Run stopped after: ", minutes, " minutes")

### 1. check md5 for RData objects ####
message("\nChecking md5 digest for RData files")
source("../../check_rdata_md5.R")

### 2. check attriton table ####
message("\nChecking attrition tables")
cmd <- "Rscript ../../check_stats.R > testrun.attr"
##message(cmd)
system(cmd)

attr.old <- read.table("intSiteValidation.attr", header=TRUE)
attr.old$workdir <- NULL
attr.new <- read.table("testrun.attr", header=TRUE)
attr.new$workdir <- NULL
if( !identical(attr.old, attr.new) ) {
    message("FAIL")
    q()
} else { message("PASS") }


q(save="no")

