concatenateImputeFiles<-function(inputStart,outputFile,boundaries)
{
	all_data<-NULL
	for(i in 1:nrow(boundaries))
	{
		filename = paste(inputStart,"_",boundaries[i,1]/1000,"K_",boundaries[i,2]/1000,"K.txt_haps",sep="")
		if(file.exists(filename) && file.info(filename)$size>0)
		{
			data<-read.table(filename, sep = " ")
			all_data<-rbind(all_data,data)
		}
	}
	write.table(all_data,outputFile, row.names=F, col.names=F, quote=F, sep=" ")
}

concatenateBAFfiles<-function(inputStart,inputEnd,outputFile,no.chrs)
{
	all_data<-NULL
	colNames<-NULL
	for(i in 1:no.chrs)
	{
		filename = paste(inputStart,i,inputEnd,sep="")
		if(file.exists(filename) && file.info(filename)$size>0)
		{
			data<-read.table(filename, sep = "\t", header=T, row.names=1)
			if(nrow(data) > 0){
				all_data<-rbind(all_data,data)
				colNames<-names(data)
			}
		}
	}
	rnames=paste("snp",1:nrow(all_data),sep="")
	write.table(all_data,outputFile, row.names=rnames, col.names=colNames, quote=F, sep="\t")
}
