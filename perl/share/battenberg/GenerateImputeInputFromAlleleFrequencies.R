args=commandArgs(TRUE)
lib_path<-toString(args[1])
impute_input_file<-toString(args[2])
is.male<-as.logical(args[3])
SNP_file<-toString(args[4])
Normal_SNP_file<-toString(args[5])
outFileStart<-toString(args[6])
chr<-as.numeric(args[7])
problemLociFile<-toString(args[8])
useLociFile = NA
heterozygousFilter=0.1
if(length(args)>=9){
	problemLociFile = toString(args[9])
	if(length(args)>=10){
		useLociFile = toString(args[10])
		if(length(args)>=11){
			heterozygousFilter=as.numeric(args[11])
		}
	}
}

#260213 - read allele frequency file for just one chromosome - requires less memory
SNP_file = gsub(".txt",paste("_chr",chr,".txt",sep=""),SNP_file)
Normal_SNP_file = gsub(".txt",paste("_chr",chr,".txt",sep=""),Normal_SNP_file)

outfile=paste(outFileStart,chr,".txt",sep="")

impute.info = read.table(impute_input_file,header=F,row.names=NULL,sep="\t",stringsAsFactors=F)
if(is.male){
	impute.info = impute.info[impute.info[,7]==1,]
}
chr_names=unique(impute.info[,1])
impute.info = impute.info[impute.info[,1]==chr_names[chr],]

print(paste("GenerateImputeInput is.male? ",is.male,sep=""))
print(paste("GenerateImputeInput #impute files? ",nrow(impute.info),sep=""))

known_SNPs<-read.table(impute.info[1,2],sep=" ",header=T)
if(nrow(impute.info)>1){
	for(r in 2:nrow(impute.info)){
		known_SNPs<-rbind(known_SNPs,read.table(impute.info[r,2],sep=" ",header=T))
	}
}

#filter out bad SNPs (streaks in BAF)
if((problemLociFile !="NA") & (!is.na(problemLociFile)))
{
	problemSNPs=read.table(problemLociFile,header=T,sep="\t")
	problemSNPs=problemSNPs$Pos[problemSNPs$Chr==chr_names[chr]]
	badIndices=match(known_SNPs[,2],problemSNPs)
	known_SNPs = known_SNPs[is.na(badIndices),]
}

#filter 'good' SNPs (e.g. SNP6 positions)
if((useLociFile !="NA") & (!is.na(useLociFile)))
{
	goodSNPs=read.table(useLociFile,header=T,sep="\t",row.names=NULL)
	goodSNPs=goodSNPs$pos[goodSNPs$chr==chr_names[chr]]
	len=length(goodSNPs)
	goodIndices=match(known_SNPs[,2],goodSNPs)
	known_SNPs = known_SNPs[!is.na(goodIndices),]
}

snp_data<-read.table(SNP_file,comment.char="#",sep="\t",header=F,stringsAsFactors=F)
normal_snp_data<-read.table(Normal_SNP_file,comment.char="#",sep="\t",header=F,stringsAsFactors=F)
snp_data<-cbind(snp_data,normal_snp_data)
indices<-match(known_SNPs[,2],snp_data[,2])
found_snp_data=snp_data[indices[!is.na(indices)],]

normalGenotype = substr(found_snp_data[1:nrow(found_snp_data),6],1,2)
refLoci = which(normalGenotype == paste(known_SNPs[!is.na(indices),3],known_SNPs[!is.na(indices),3],sep=""))
hetLoci = which(normalGenotype == paste(known_SNPs[!is.na(indices),3],known_SNPs[!is.na(indices),4],sep=""))
hetLoci = c(hetLoci,which(normalGenotype == paste(known_SNPs[!is.na(indices),4],known_SNPs[!is.na(indices),3],sep="")))
altLoci = which(normalGenotype == paste(known_SNPs[!is.na(indices),4],known_SNPs[!is.na(indices),4],sep=""))

genotypes = array(0,c(sum(!is.na(indices)),3))

genotypes[refLoci,1]=1
genotypes[hetLoci,2]=1
genotypes[altLoci,3]=1

minBaf=min(heterozygousFilter,1.0-heterozygousFilter)
maxBaf=max(heterozygousFilter,1.0-heterozygousFilter)
nucleotides=c("A","C","G","T")
print(paste("nrow knownsnps=",nrow(known_SNPs),sep=""))
print(paste("ncol knownsnps=",ncol(known_SNPs),sep=""))
print(paste("nrow found_snp_data=",nrow(found_snp_data),sep=""))
print(paste("ncol found_snp_data=",ncol(found_snp_data),sep=""))

ref_indices=match(known_SNPs[!is.na(indices),3],nucleotides)+ncol(normal_snp_data)+2
alt_indices=match(known_SNPs[!is.na(indices),4],nucleotides)+ncol(normal_snp_data)+2

BAFs<-as.numeric(found_snp_data[cbind(1:nrow(found_snp_data),alt_indices)])/(as.numeric(found_snp_data[cbind(1:nrow(found_snp_data),alt_indices)])+as.numeric(found_snp_data[cbind(1:nrow(found_snp_data),ref_indices)]))

genotypes[BAFs<=minBaf,1]=1
genotypes[BAFs>minBaf & BAFs<maxBaf,2]=1
genotypes[BAFs>=maxBaf,3]=1
snp.names=paste("snp",1:sum(!is.na(indices)),sep="")
out.data<-cbind(snp.names,known_SNPs[!is.na(indices),1:4],genotypes)

write.table(out.data,file=outfile,row.names=F,col.names=F,quote=F)
if(is.na(as.numeric(chr_names[chr]))){
	sample.g.file=paste(outFileStart,"sample_g.txt",sep="")
	#not sure this is necessary, because only the PAR regions are used for males
	#if(is.male){
	#	sample_g_data=data.frame(ID_1=c(0,"INDIVI1"),ID_2=c(0,"INDIVI1"),missing=c(0,0),sex=c("D",1))
	#}else{
		sample_g_data=data.frame(ID_1=c(0,"INDIVI1"),ID_2=c(0,"INDIVI1"),missing=c(0,0),sex=c("D",2))
	#}
	write.table(sample_g_data,file=sample.g.file,row.names=F,col.names=T,quote=F)
}

q(save="no")
