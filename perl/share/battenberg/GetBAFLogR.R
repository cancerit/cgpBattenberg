##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
#
# This file is part of cgpBattenberg.
#
# cgpBattenberg is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##########

#getBAFsAndLogRs<-function(imputeInfoFile,inputFile, normalInputFile,  BAFnormalFile, BAFmutantFile, logRnormalFile, logRmutantFile,minCounts=NA,samplename="sample1")
getBAFsAndLogRs<-function(imputeInfoFile,inputFile, normalInputFile,  BAFnormalFile, BAFmutantFile, logRnormalFile, logRmutantFile,minCounts=NA,samplename="sample1",thougenloc)
{
	impute.info = read.table(imputeInfoFile,header=F,row.names=NULL,sep="\t",stringsAsFactors=F)
	chr_names=unique(impute.info[,1])

	input_data=NULL
	normal_input_data=NULL
	allele_data=NULL
	for(chr in 1:length(chr_names)){
		chr_file = gsub(".txt",paste("_chr",chr,".txt",sep=""),inputFile)
		normal_chr_file = gsub(".txt",paste("_chr",chr,".txt",sep=""),normalInputFile)

		input_data<-rbind(input_data,read.table(chr_file,header=TRUE,sep="\t",comment.char=""))
		normal_input_data<-rbind(normal_input_data,read.table(normal_chr_file,header=TRUE,sep="\t",comment.char=""))

		#allele_data<-rbind(allele_data,read.table(paste("1000genomesloci/1000genomesAlleles2012_chr",chr,".txt",sep=""),sep="\t",header=T))
		allele_data<-rbind(allele_data,read.table(paste(thougenloc,"/1000genomesAlleles2012_chr",chr,".txt",sep=""),sep="\t",header=T))
	}

	print("data dims:")
	print(dim(input_data))
	print(dim(normal_input_data))

	head(input_data)
	head(normal_input_data)
	head(allele_data)

	names(input_data)[1] = "CHR"
	names(input_data)[2] = "POS"
	names(normal_input_data)[1] = "CHR"
	names(normal_input_data)[2] = "POS"

	normal_data = normal_input_data[,3:6]
	mutant_data = input_data[,3:6]

	len <- nrow(normal_data)
	normCount1 <-normal_data[cbind(1:len,allele_data[,2])]
	normCount2 <-normal_data[cbind(1:len,allele_data[,3])]
	totalNormal <- normCount1 + normCount2
	mutCount1 <-mutant_data[cbind(1:len,allele_data[,2])]
	mutCount2 <-mutant_data[cbind(1:len,allele_data[,3])]
	totalMutant <- mutCount1 + mutCount2

	indices=1:nrow(input_data)
	if(!is.na(minCounts)){
		print(paste("minCount=",minCount,sep=""))
		#indices=which(totalNormal>=minCounts & totalMutant>=minCounts)
		#only normal has to have min coverage, mutant must have at least 1 read to prevent division by zero
		indices=which(totalNormal>=minCounts & totalMutant>=1)
		#only tumour has to have min coverage, normal must have at least 1 read to prevent division by zero
		#indices=which(totalMutant>=minCounts & totalNormal>=1)

		totalNormal = totalNormal[indices]
		totalMutant = totalMutant[indices]
		normal_data = normal_data[indices,]
		mutant_data = mutant_data[indices,]

		normCount1 = normCount1[indices]
		normCount2 = normCount2[indices]
		mutCount1 = mutCount1[indices]
		mutCount2 = mutCount2[indices]
	}
	n<-length(indices)

	normalBAF=vector(length=n,mode="numeric")
	mutantBAF=vector(length=n,mode="numeric")
	normalLogR=vector(length=n,mode="numeric")
	mutantLogR=vector(length=n,mode="numeric")

	#randomise A and B alleles
	selector<-round(runif(n))

	normalBAF[which(selector==0)] = normCount1[which(selector==0)] / totalNormal[which(selector==0)]
	normalBAF[which(selector==1)] = normCount2[which(selector==1)] / totalNormal[which(selector==1)]

	mutantBAF[which(selector==0)] = mutCount1[which(selector==0)] / totalMutant[which(selector==0)]
	mutantBAF[which(selector==1)] = mutCount2[which(selector==1)] / totalMutant[which(selector==1)]

	normalLogR = vector(length=n,mode="integer") #assume that normallogR is 0, and normalise mutantLogR to normalLogR
	mutantLogR = totalMutant/totalNormal

	Germline_BAF = data.frame(CHR=input_data$CHR[indices],POS=input_data$POS[indices],baf=normalBAF)
	colnames(Germline_BAF)=c("Chromosome","Position",samplename)
	rownames(Germline_BAF)=paste("snp",1:(dim(Germline_BAF)[1]),sep="")
	Germline_LogR = data.frame(CHR=input_data$CHR[indices],POS=input_data$POS[indices],samplename=normalLogR)
	colnames(Germline_LogR)=c("Chromosome","Position",samplename)
	rownames(Germline_LogR)=paste("snp",1:(dim(Germline_BAF)[1]),sep="")
	Tumor_BAF = data.frame(CHR=input_data$CHR[indices],POS=input_data$POS[indices],baf=mutantBAF)
	colnames(Tumor_BAF)=c("Chromosome","Position",samplename)
	rownames(Tumor_BAF)=paste("snp",1:(dim(Germline_BAF)[1]),sep="")
	Tumor_LogR = data.frame(CHR=input_data$CHR[indices],POS=input_data$POS[indices],samplename=log2(mutantLogR/mean(mutantLogR,na.rm=T)))
	colnames(Tumor_LogR)=c("Chromosome","Position",samplename)
	rownames(Tumor_LogR)=paste("snp",1:(dim(Germline_BAF)[1]),sep="")

	write.table(Germline_BAF,file=BAFnormalFile,append=FALSE,quote=FALSE,sep="\t",col.names=c("Chromosome","Position",samplename))
	write.table(Tumor_BAF,file=BAFmutantFile,append=FALSE,quote=FALSE,sep="\t",col.names=c("Chromosome","Position",samplename))
	write.table(Germline_LogR,file=logRnormalFile,append=FALSE,quote=FALSE,sep="\t",col.names=c("Chromosome","Position",samplename))
	write.table(Tumor_LogR,file=logRmutantFile,append=FALSE,quote=FALSE,sep="\t",col.names=c("Chromosome","Position",samplename))
}
