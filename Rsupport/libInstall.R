args<-commandArgs(TRUE)
instLib = args[1]
ascatPackage = args[2]
battenbergPackage = args[3]

#Use previous version because current version in cran is not compatible with R < 3.2 (dir.exists).
stringi_for_legacyR = "1.1.2"

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    biocLite(new.pkg, ask=FALSE, lib=instLib)
  sapply(pkg, library, character.only = TRUE)
}

r = getOption("repos") # hard code the UK repo for CRAN
r["CRAN"] = "http://cran.uk.r-project.org"
options(repos = r)
rm(r)
source("http://bioconductor.org/biocLite.R")

tmp <- c("devtools")
ipak(tmp)
library(devtools)
options(download.file.method = "auto")

if ( version$major >= 3 && version$minor >= 2 ) {
  biocPackages <- c("stringi", "readr", "doParallel", "ggplot2", "RColorBrewer", "gridExtra", "gtools")
} else {
  install_version("stringi", stringi_for_legacyR)
  biocPackages <- c("readr", "doParallel", "ggplot2", "RColorBrewer", "gridExtra", "gtools")
}

ipak(biocPackages)
install_file(ascatPackage)
install_file(battenbergPackage)
