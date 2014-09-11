##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
# 
# Author: Cancer Genome Project cgpit@sanger.ac.uk
# 
# This file is part of battenberg.
# 
# battenberg is free software: you can redistribute it and/or modify it under
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

args = ( commandArgs(TRUE) )
lib_path<-toString(args[1])
samplename<-toString(args[2])
start.file<-toString(args[3])
dist_choice <- as.integer(args[4])
ascat_dist_choice <- as.integer(args[5])
uninformative_BAF_threshold = 0.51
gamma_param=1
use_preset_rho_psi=F
if( length( args ) >= 6)
{
	uninformative_BAF_threshold <- as.numeric(args[6] )
	if( length( args ) >= 7 )
	{
		gamma_param=as.numeric(args[7])
		if( length( args ) >= 8 ){
			min.ploidy = as.numeric(args[8])
			max.ploidy = as.numeric(args[9])
			min.rho = as.numeric(args[10])
			min.goodness.of.fit = as.numeric(args[11])
			if( length( args ) >= 12 ){
				preset_rho = as.numeric(args[12])
				preset_psi = as.numeric(args[13])
				use_preset_rho_psi=T
			}
		}
	}
}

read_depth=30

source(paste(lib_path,"ascat.R",sep="/"))

segmented.BAF.data = read.table(paste(samplename,".BAFsegmented.txt",sep=""),sep="\t",header=F,stringsAsFactors=F,skip=1,row.names=1)
raw.BAF.data = read.table(paste(start.file,"mutantBAF.tab",sep=""),sep="\t",header=T,stringsAsFactors=F)
raw.logR.data = read.table(paste(start.file,"mutantLogR.tab",sep=""),sep="\t",header=T,stringsAsFactors=F)

#300713
raw.BAF.data = raw.BAF.data[!is.na(raw.BAF.data[,3]),]
raw.logR.data = raw.logR.data[!is.na(raw.logR.data[,3]),]

#chromosome names are sometimes 'chr1', etc.
if(length(grep("chr",raw.BAF.data[1,1]))>0){
	raw.BAF.data[,1] = gsub("chr","",raw.BAF.data[,1])
}
if(length(grep("chr",raw.logR.data[1,1]))>0){
	raw.logR.data[,1] = gsub("chr","",raw.logR.data[,1])
}

BAF.data = NULL
logR.data = NULL
segmented.logR.data = NULL
chr.names = unique(segmented.BAF.data[,1])
matched.segmented.BAF.data = NULL
for(chr in chr.names){
	chr.BAF.data = raw.BAF.data[raw.BAF.data$Chromosome==chr,]
	if(nrow(chr.BAF.data)==0){next}
	chr.segmented.BAF.data = segmented.BAF.data[segmented.BAF.data[,1]==chr,]
	print(chr.segmented.BAF.data[1:3,])
	indices = match(chr.segmented.BAF.data[,2],chr.BAF.data$Position )

	#130313
	chr.segmented.BAF.data = chr.segmented.BAF.data[!is.na(indices),]
	matched.segmented.BAF.data = rbind(matched.segmented.BAF.data, chr.segmented.BAF.data)

	BAF.data = rbind(BAF.data, chr.BAF.data[indices[!is.na(indices)],])
	chr.logR.data = raw.logR.data[raw.logR.data$Chromosome==chr,]
	indices = match(chr.segmented.BAF.data[,2],chr.logR.data$Position)
	logR.data = rbind(logR.data, chr.logR.data[indices[!is.na(indices)],])
	chr.segmented.logR.data = chr.logR.data[indices[!is.na(indices)],]
	segs = make_seg_lr(chr.segmented.BAF.data[,5])
	cum.segs = c(0,cumsum(segs))
	for(s in 1:length(segs)){
		chr.segmented.logR.data[(cum.segs[s]+1):cum.segs[s+1],3] = mean(chr.segmented.logR.data[(cum.segs[s]+1):cum.segs[s+1],3])
	}
	segmented.logR.data = rbind(segmented.logR.data,chr.segmented.logR.data)
	#print(paste(nrow(matched.segmented.BAF.data),nrow(segmented.logR.data),nrow(logR.data),sep=","))
}
names(matched.segmented.BAF.data)[5] = samplename
row.names(segmented.logR.data) = row.names(matched.segmented.BAF.data)
row.names(logR.data) = row.names(matched.segmented.BAF.data)

write.table(segmented.logR.data,paste(samplename,".logRsegmented.txt",sep=""),sep="\t",quote=F,col.names=F)

segBAF = 1-matched.segmented.BAF.data[,5]
segLogR = segmented.logR.data[,3]
logR = logR.data[,3]
names(segBAF) = rownames(matched.segmented.BAF.data)
names(segLogR) = rownames(matched.segmented.BAF.data)
names(logR) = rownames(matched.segmented.BAF.data)
print(unique(logR.data[,1]))
chr.segs = NULL
for(ch in 1:length(chr.names)){
	chr.segs[[ch]] = which(logR.data[,1]==chr.names[ch])
}

if(use_preset_rho_psi){
	ascat_optimum_pair = list(rho=preset_rho, psi = preset_psi)
}else{
	distance.outfile=paste(start.file,"distance.png",sep="",collapse="") # kjd 20-2-2014
	copynumberprofile.outfile=paste(start.file,"copynumberprofile.png",sep="",collapse="") # kjd 20-2-2014
	nonroundedprofile.outfile=paste(start.file,"nonroundedprofile.png",sep="",collapse="") # kjd 20-2-2014

	ascat_optimum_pair=runASCAT(logR, 1-BAF.data[,3], segLogR, segBAF, chr.segs, ascat_dist_choice,distance.outfile, copynumberprofile.outfile, nonroundedprofile.outfile, gamma=gamma_param, allow100percent=T, reliabilityFile=NA, min.ploidy, max.ploidy, min.rho, min.goodness.of.fit) # kjd 4-2-2014
}

distance.outfile=paste(start.file,"second_distance.png",sep="",collapse="") # kjd 20-2-2014
copynumberprofile.outfile=paste(start.file,"second_copynumberprofile.png",sep="",collapse="") # kjd 20-2-2014
nonroundedprofile.outfile=paste(start.file,"second_nonroundedprofile.png",sep="",collapse="") # kjd 20-2-2014

out=run_clonal_ASCAT( logR, 1-BAF.data[,3], segLogR, segBAF, chr.segs, matched.segmented.BAF.data, ascat_optimum_pair, dist_choice, distance.outfile, copynumberprofile.outfile, nonroundedprofile.outfile, gamma_param=gamma_param, read_depth, uninformative_BAF_threshold, allow100percent=T, reliabilityFile=NA) # kjd 21-2-2014

ascat_optimum_pair_fraction_of_genome = out$output_optimum_pair_without_ref
ascat_optimum_pair_ref_seg = out$output_optimum_pair
is.ref.better = out$is.ref.better

rho_psi_output = data.frame(rho = c(ascat_optimum_pair$rho,ascat_optimum_pair_fraction_of_genome$rho,ascat_optimum_pair_ref_seg$rho),psi = c(ascat_optimum_pair$psi,ascat_optimum_pair_fraction_of_genome$psi,ascat_optimum_pair_ref_seg$psi), distance = c(NA,out$distance_without_ref,out$distance), is.best = c(NA,!is.ref.better,is.ref.better),row.names=c("ASCAT","FRAC_GENOME","REF_SEG"))
write.table(rho_psi_output,paste(start.file,"rho_and_psi.txt",sep=""),quote=F,sep="\t")

q(save="no")
