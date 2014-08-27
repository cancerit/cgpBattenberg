args=(commandArgs(TRUE))
lib_path<-toString(args[1])
imputeInfoFile<-toString(args[2])
inputFile<-toString(args[3])
normalInputFile<-toString(args[4])
outputDir<-toString(args[5])
minCount = as.numeric(args[6])
samplename = toString(args[7])
thougenloc<-toString(args[8])

source(paste(lib_path,"GetBAFLogR.R",sep="/"))
source(paste(lib_path,"ascat.R",sep="/"))

getBAFsAndLogRs(imputeInfoFile,inputFile,normalInputFile,paste(outputDir,"normalBAF.tab",sep="",collapse=""),paste(outputDir,"mutantBAF.tab",sep="",collapse=""),paste(outputDir,"normalLogR.tab",sep="",collapse=""),paste(outputDir,"mutantLogR.tab",sep="",collapse=""),minCount=minCount,samplename=samplename,thougenloc=thougenloc)

ascat.bc = ascat.loadData(paste(outputDir,"mutantLogR.tab",sep=""),paste(outputDir,"mutantBAF.tab",sep=""),paste(outputDir,"normalLogR.tab",sep=""),paste(outputDir,"normalBAF.tab",sep=""))
ascat.plotRawData(ascat.bc, parentDir=outputDir)

q(save="no")
