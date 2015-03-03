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

args=commandArgs(TRUE)
lib_path<-toString(args[1])
imputeInfoFile<-toString(args[2])
impute.executable<-toString(args[3])
is.male<-as.logical(args[4])
inFileStart<-toString(args[5])
outFileStart<-toString(args[6])
chr<-as.numeric(args[7])

inpute.info = read.table(imputeInfoFile,header=F,row.names=NULL,sep="\t",stringsAsFactors=F)
if(is.male){
	inpute.info = inpute.info[inpute.info[,7]==1,]
}
chr_names=unique(inpute.info[,1])
inpute.info = inpute.info[inpute.info[,1]==chr_names[chr],]

for(r in 1:nrow(inpute.info)){
	boundaries = seq(as.numeric(inpute.info[r,5]),as.numeric(inpute.info[r,6]),5000000)
	if(boundaries[length(boundaries)] != inpute.info[r,6]){
		boundaries = c(boundaries,inpute.info[r,6])
	}
	boundaries=boundaries
	for(b in 1:(length(boundaries)-1)){
		cmd = paste(impute.executable," -m ",inpute.info[r,3]," -h ",inpute.info[r,4]," -l ",inpute.info[r,2]," -g ",inFileStart,chr,".txt -int ",boundaries[b]," ",boundaries[b+1]," -Ne 20000 -o ",outFileStart,chr,"_",boundaries[b]/1000,"K_",boundaries[b+1]/1000,"K.txt -phase",sep="")
		system(cmd, wait=T)
	}
}
#allow time for all output files to be written
#Sys.sleep(10)

q(save="no")
