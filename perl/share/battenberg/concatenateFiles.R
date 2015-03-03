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
