instLib = commandArgs(T)[1]

r = getOption("repos") # hard code the UK repo for CRAN
r["CRAN"] = "http://cran.uk.r-project.org"
options(repos = r)
rm(r)

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    biocLite(new.pkg, ask=FALSE, lib=instLib, lib.loc=instLib)
  sapply(pkg, library, character.only = TRUE)
}

ipak_bioc <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    BiocManager::install(new.pkg, ask=FALSE, lib=instLib, lib.loc=instLib)
  sapply(pkg, library, character.only = TRUE)
}

if( (version$major == 3 && version$minor >=5) || version$major > 3) {
  # biocmanager versions of R
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(ask=FALSE, lib=instLib, lib.loc=instLib)
  ipak_bioc(c("readr"))
  ipak_bioc(c("ggplot2"))
  ipak_bioc(c("doParallel"))
  ipak_bioc(c("gridExtra"))
  ipak_bioc(c("gtools"))
  ipak_bioc(c("RColorBrewer"))
} else {
  # OLD versions of R
  source("http://bioconductor.org/biocLite.R")
  ipak(c("readr"))
  ipak(c("ggplot2"))
  ipak(c("doParallel"))
  ipak(c("gridExtra"))
  ipak(c("gtools"))
  ipak(c("RColorBrewer"))
}

# works on old and new
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", lib=instLib)
library(devtools)
options(download.file.method = "auto")
install_github("Irrationone/copynumber", ref="87d2663fe6b11c03cf6006b4ee9ed70450eacb5a")
