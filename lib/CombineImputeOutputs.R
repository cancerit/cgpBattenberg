args=(commandArgs(TRUE)) 
impute_input_file<-toString(args[1])
is.male<-as.logical(args[2])
inputStart<-toString(args[3])
chr<-as.numeric(args[4])

impute.info = read.table(impute_input_file,header=F,row.names=NULL,sep="\t",stringsAsFactors=F)
if(is.male){
	impute.info = impute.info[impute.info[,7]==1,]
}
chr_names=unique(impute.info[,1])
impute.info = impute.info[impute.info[,1]==chr_names[chr],]

all.boundaries = array(0,c(0,2))
for(r in 1:nrow(impute.info)){
	boundaries = seq(as.numeric(impute.info[r,5]),as.numeric(impute.info[r,6]),5000000)
	if(boundaries[length(boundaries)] != impute.info[r,6]){
		boundaries = c(boundaries,impute.info[r,6])
	}
	all.boundaries = rbind(all.boundaries,cbind(boundaries[-(length(boundaries))],boundaries[-1]))
}

source("concatenateFiles.R")
concatenateImputeFiles(paste(inputStart,chr,sep=""),paste(inputStart,chr,"_allHaplotypeInfo.txt",sep=""),all.boundaries)
q(save="no")