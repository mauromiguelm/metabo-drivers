#### save global environment objects

globalEnv <- ls()

#### Processing of Echo exception log, from source plate

fileNames <- list.files(path = 'data/exceptions')

fileNames <- 
fileNames[grepl(pattern = "E5XX-18041", x = fileNames, ignore.case = T) ]

exceptions_source <- 
  lapply(fileNames, function(x) read.table(paste0('data/exceptions/',x), sep = "\t", header = T, stringsAsFactors = F))

exceptions_source <- do.call(rbind, exceptions_source)

# Processing of Echo exception log, from target plate

fileNames <- list.files('data/exceptions')

fileNames <- 
  fileNames[grepl(pattern = "exceptions_report", x = fileNames, ignore.case = T) ]


exceptions_target <- 
  lapply(fileNames, function(x) read.table(paste0('data/exceptions/',x), sep = ",", header = T, stringsAsFactors = F))

exceptions_target <- 
  do.call(rbind, exceptions_target)

exceptions_target$cellPlateBC <- gsub(".mat", "", exceptions_target$cellPlateBC)

exceptions <- list(exceptions_source, exceptions_target)

rm(list = ls()[!ls()%in% c(globalEnv,"exceptions")])
