args<-commandArgs(TRUE)
instLib = args[1]
ascatPackage = args[2]
battenbergPackage = args[3]

r = getOption("repos") # hard code the UK repo for CRAN
r["CRAN"] = "http://cran.uk.r-project.org"
options(repos = r)
rm(r)
#Use previous version because current version in cran is not compatible with R < 3.2 (dir.exists).
install.packages("http://cran.uk.r-project.org/src/contrib/Archive/stringi/stringi_1.1.2.tar.gz", repos=NULL, type="source")
install.packages("RColorBrewer")
install.packages("ggplot2")
install.packages("readr")
install.packages(ascatPackage)
install.packages(battenbergPackage)
