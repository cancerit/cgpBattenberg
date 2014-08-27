args=(commandArgs(TRUE))
lib_path<-toString(args[1])
impute_input_file<-toString(args[2])
is.male<-as.logical(args[3])
inputStart<-toString(args[4])
chr<-as.numeric(args[5])

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
source(paste(lib_path,"concatenateFiles.R",sep="/"))
concatenateImputeFiles(paste(inputStart,chr,sep=""),paste(inputStart,chr,"_allHaplotypeInfo.txt",sep=""),all.boundaries)
q(save="no")
