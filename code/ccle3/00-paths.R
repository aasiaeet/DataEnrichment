###############################
### Set the root directory, which is where we can find the CyTOF
### data on the local machine. Need to do this before 'source'ing
### the 'overlay' code

home <- Sys.getenv("HOME", unset = NA)
if (is.na(home)) stop("Cannot find 'HOME' from environment variable s.")

jinfo <- file.path(home, "Paths", "dataenrichment.json")
if (!file.exists(jinfo)) stop("Cannot locate file: '", jinfo, "'.\n", sep='')
library(rjson)
temp <- fromJSON(file = jinfo)
paths <- temp$paths

rm(home, jinfo, temp)
