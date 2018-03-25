args<-commandArgs(TRUE)
instLib = args[1]
ascatPackage = args[2]
battenbergPackage = args[3]

#Use earier version because current version in cran is not compatible with R < 3.2 (dir.exists).
stringi_for_legacyR = "http://cran.uk.r-project.org/src/contrib/Archive/stringi/stringi_1.1.6.tar.gz"

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    biocLite(new.pkg, ask=FALSE, lib=instLib, lib.loc=instLib)
  sapply(pkg, library, character.only = TRUE)
}

r = getOption("repos") # hard code the UK repo for CRAN
r["CRAN"] = "http://cran.uk.r-project.org"
options(repos = r)
rm(r)
source("http://bioconductor.org/biocLite.R")

biocPackages <- c("devtools")
ipak(biocPackages)
library(devtools)
options(download.file.method = "auto")

if ( version$major > 3 || ( version$major == 3 && version$minor >= 2 )) {
  ipak(c("stringi"))
} else {
  install.packages(stringi_for_legacyR, repos=NULL, type="source")
}

ipak("readr")
ipak("ggplot2")
ipak("doParallel")
ipak("RColorBrewer")
ipak("gridExtra")
ipak("gtools")

install.packages(ascatPackage)
install.packages(battenbergPackage)
