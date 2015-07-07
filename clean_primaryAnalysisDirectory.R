keptpattern <- c("*.txt$",
                 "*.csv$",
                 "*.tsv$",
                 "logs/",
                 "html$",
                 "completeMetadata.RData",
                 "sites.final.RData",
                 "allSites.RData",
                 "rawSites.RData",
                 "multihitData.RData",
                 "chimeraData.RData",
                 "stats.RData" )

allfiles <- list.files(recursive=TRUE)

filelist <- sapply(keptpattern, function(pattern) grep(pattern, allfiles, value=TRUE))

filekeep <- sort(unique(unlist(unname(filelist))))

filedel <- setdiff(allfiles, filekeep)

message("The following files will be deleted\n",
        paste(filedel, collapse="\n") )

file.remove(filedel)

