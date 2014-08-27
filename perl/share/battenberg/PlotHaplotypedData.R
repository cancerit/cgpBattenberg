args=(commandArgs(TRUE))
lib_path<-toString(args[1])
imputeInfoFile<-toString(args[2])
chr<-as.numeric(args[3])
mutFile<-toString(args[4])
samplename<-toString(args[5])

impute.info = read.table(imputeInfoFile,header=F,row.names=NULL,sep="\t",stringsAsFactors=F)
chr_names=unique(impute.info[,1])
chr_name = chr_names[chr]

imageFileName = paste(samplename,"_chr",chr_name,"_heterozygousData.png",sep="")

mut_data<-read.table(mutFile,sep="\t",header=T)

png(filename = imageFileName, width = 10000, height = 2500, res = 500)
par(pch = ".", cex=1, cex.main=0.8, cex.axis = 0.6, cex.lab=0.7,yaxp=c(-0.05,1.05,6))

	if(nrow(mut_data) > 0){
		plot(c(min(mut_data$Position,na.rm=T),max(mut_data$Position,na.rm=T)),c(0,1),type="n",,main=paste(samplename,", chromosome",mut_data[1,1],sep=" "),xlab="pos", ylab="BAF")
		points(mut_data$Position,mut_data[,3],col="blue")
		points(mut_data$Position,1-mut_data[,3],col="red")
	}
dev.off()

q(save="no")
