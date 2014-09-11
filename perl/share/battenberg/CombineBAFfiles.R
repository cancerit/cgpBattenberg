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

args=(commandArgs(TRUE))
lib_path<-toString(args[1])
impute_info_file<-toString(args[2])
is.male<-as.logical(args[3])
inputStart<-toString(args[4])
inputEnd<-toString(args[5])
outputFile<-toString(args[6])
source(paste(lib_path,"concatenateFiles.R",sep="/"))

impute.info = read.table(impute_info_file,header=F,row.names=NULL,sep="\t",stringsAsFactors=F)
if(is.male){
	impute.info = impute.info[impute.info[,7]==1,]
}
no.chrs = length(unique(impute.info[,1]))
concatenateBAFfiles(inputStart,inputEnd,outputFile,no.chrs)
q(save="no")
