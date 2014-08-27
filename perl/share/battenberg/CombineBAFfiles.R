args=(commandArgs(TRUE))
lib_path<-toString(args[1])
impute_info_file<-toString(args[2])
is.male<-as.logical(args[3])
inputStart<-toString(args[4])
inputEnd<-toString(args[5])
outputFile<-toString(args[6])
source(paste(lib_path,"concatenateFiles.R",sep="/"))

impute.info = read.table(impute_info_file,header=F,row.names=NULL,sep="\t",stringsAsFactors=F)
if(is.male){
	impute.info = impute.info[impute.info[,7]==1,]
}
no.chrs = length(unique(impute.info[,1]))
concatenateBAFfiles(inputStart,inputEnd,outputFile,no.chrs)
q(save="no")
