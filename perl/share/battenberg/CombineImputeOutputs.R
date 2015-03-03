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
