args=(commandArgs(TRUE))
imputeInfoFile<-toString(args[1])
inputFile<-toString(args[2])
normalInputFile<-toString(args[3])
outputDir<-toString(args[4])
minCount = as.numeric(args[5])
samplename = toString(args[6])

source("GetBAFLogR.R")
source("ascat.R")

getBAFsAndLogRs(imputeInfoFile,inputFile,normalInputFile,paste(outputDir,"normalBAF.tab",sep="",collapse=""),paste(outputDir,"mutantBAF.tab",sep="",collapse=""),paste(outputDir,"normalLogR.tab",sep="",collapse=""),paste(outputDir,"mutantLogR.tab",sep="",collapse=""),minCount=minCount,samplename=samplename)

ascat.bc = ascat.loadData(paste(outputDir,"mutantLogR.tab",sep=""),paste(outputDir,"mutantBAF.tab",sep=""),paste(outputDir,"normalLogR.tab",sep=""),paste(outputDir,"normalBAF.tab",sep=""))
ascat.plotRawData(ascat.bc, parentDir=outputDir)

q(save="no")
