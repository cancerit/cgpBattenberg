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

args=(commandArgs(TRUE))

lib_path<-toString(args[1])
impute_info_file<-toString(args[2])
is.male<-as.logical(args[3])
SNP_file<-toString(args[4])
Normal_SNP_file<-toString(args[5])
hapFileStart<-toString(args[6])
samplename<-toString(args[7])
BAFfileStart<-toString(args[8])
chr<-toString(args[9])
minCount = 10
if(length(args)>=10){
	minCount<-as.integer(args[10])
}

source(paste(lib_path,"GetHaplotypedBAFs.R",sep="/"))

impute.info = read.table(impute_info_file,header=F,row.names=NULL,sep="\t",stringsAsFactors=F)
if(is.male){
	impute.info = impute.info[impute.info[,7]==1,]
}
chr_names=unique(impute.info[,1])

#260213 - read allele frequency file for just one chromosome - requires less memory
SNP_file = gsub(".txt",paste("_chr",chr,".txt",sep=""),SNP_file)
#Normal_SNP_file = gsub(".txt",paste("_chr",chr,".txt",sep=""),Normal_SNP_file)

GetChromosomeBAFs(chr,SNP_file, hapFileStart,"_allHaplotypeInfo.txt",samplename,paste(BAFfileStart,"chr",chr,"_heterozygousMutBAFs_haplotyped.txt",sep=""),chr_names,minCount)

q(save="no")
