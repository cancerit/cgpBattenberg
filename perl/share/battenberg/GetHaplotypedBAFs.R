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

#gets b-allele frequencies for one chromosome and writes them in the format required by ASCAT
GetChromosomeBAFs<-function(chr, SNP_file, startHaplotypeFile,endHaplotypeFile, samplename, outfile, chr_names, minCounts=10)
{
	hetMutBAFs<-NULL

	snp_data<-read.table(SNP_file,comment.char="",sep="\t",header=T,stringsAsFactors=F)

	print(snp_data[1:3,])
	print(chr_names)
	print(chr)

	variant_data<-read.table(paste(startHaplotypeFile,chr,endHaplotypeFile,sep=""),header=F)
	#just select heterozygous SNPs
	het_variant_data=variant_data[variant_data[,6] != variant_data[,7],]

	chr_name = chr_names[as.numeric(chr)]
	print(chr_name)

	indices<-match(het_variant_data[,3],snp_data[,2])
	het_variant_data=het_variant_data[!is.na(indices),]
	snp_indices<-indices[!is.na(indices)]
	filtered_snp_data = snp_data[snp_indices,]
	if(nrow(het_variant_data)==0){
		write.table(array(NA,c(0,3)),outfile,sep="\t",col.names=c("Chromosome","Position",samplename),quote=F)
		return()
	}
	nucleotides=c("A","C","G","T")
	ref_indices=match(het_variant_data[cbind(1:nrow(het_variant_data),4+het_variant_data[,6])],nucleotides)
	alt_indices=match(het_variant_data[cbind(1:nrow(het_variant_data),4+het_variant_data[,7])],nucleotides)

	print(filtered_snp_data[1:3,])

	min_indices=NA
	ref.count = as.numeric(filtered_snp_data[cbind(1:nrow(filtered_snp_data),alt_indices+2)])
	alt.count = as.numeric(filtered_snp_data[cbind(1:nrow(filtered_snp_data),ref_indices+2)])
	denom<- ref.count + alt.count

	min_indices<-denom>=minCounts

	filtered_snp_data=filtered_snp_data[min_indices,]
	denom = denom[min_indices]
	alt.count = alt.count[min_indices]

	hetMutBAFs<-cbind(chr_name,filtered_snp_data[,2],alt.count/denom)
	write.table(hetMutBAFs,outfile,sep="\t",row.names=paste("snp",1:nrow(hetMutBAFs),sep=""),col.names=c("Chromosome","Position",samplename),quote=F)
}
