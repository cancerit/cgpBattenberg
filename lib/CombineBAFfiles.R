args=(commandArgs(TRUE)) 
impute_info_file<-toString(args[1])
is.male<-as.logical(args[2])
inputStart<-toString(args[3])
inputEnd<-toString(args[4])
outputFile<-toString(args[5])
source("concatenateFiles.R")

impute.info = read.table(impute_info_file,header=F,row.names=NULL,sep="\t",stringsAsFactors=F)
if(is.male){
	impute.info = impute.info[impute.info[,7]==1,]
}
no.chrs = length(unique(impute.info[,1]))
concatenateBAFfiles(inputStart,inputEnd,outputFile,no.chrs)
q(save="no")
