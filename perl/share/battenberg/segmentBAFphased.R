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
sample = toString(args[2])

source(paste(lib_path,"fastPCF.R",sep="/"))

gamma = 10
phasegamma = 3
kmin=3
phasekmin=3
if(length(args)>=3){
	gamma = as.integer(args[3])
	if(kmin>gamma){
		kmin=gamma
	}
	if(length(args)>=4){
		phasegamma = as.numeric(args[4])
		if(phasekmin>phasegamma){
			phasekmin=phasegamma
		}
		if(length(args)>=5){
			kmin = as.numeric(args[5])
			if(length(args)>=6){
				phasekmin = as.numeric(args[6])
			}
		}
	}
}

BAFraw = read.table(paste(sample,"_allChromosomes_heterozygousMutBAFs_haplotyped.txt",sep=""),sep="\t",header=T)

BAFoutput = NULL

#for (chr in sort(unique(BAFraw[,1]))) {
for (chr in unique(BAFraw[,1])) {
  BAFrawchr = BAFraw[BAFraw[,1]==chr,c(2,3)]
  BAFrawchr = BAFrawchr[!is.na(BAFrawchr[,2]),]

  BAF = BAFrawchr[,2]
  pos = BAFrawchr[,1]
  names(BAF) = rownames(BAFrawchr)
  names(pos) = rownames(BAFrawchr)

  sdev <- getMad(ifelse(BAF<0.5,BAF,1-BAF),k=25)
  #DCW 250314
  #for cell lines, sdev goes to zero in regions of LOH, which causes problems.
  #0.09 is around the value expected for a binomial distribution around 0.5 with depth 30
  if(sdev<0.09){
	  sdev = 0.09
  }

  if(length(BAF)<50){
    BAFsegm = rep(mean(BAF),length(BAF))
  }else{
    res= selectFastPcf(BAF,phasekmin,phasegamma*sdev,T)
    BAFsegm = res$yhat
  }

  png(filename = paste(sample,"_RAFseg_",ifelse(grepl("chr", chr), "", "chr"),
                       chr,".png",sep=""), width = 2000, height = 1000, res = 200)
  par(mar = c(5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2, cex.lab = 2)
  plot(c(min(pos)/1000000,max(pos)/1000000),c(0,1),pch=".",type = "n",
	   main = paste(sample,", chromosome ", chr, sep=""), xlab = "Position (Mb)", ylab = "BAF (phased)")
  points(pos/1000000,BAF,pch=".",col="red",cex=2)
  points(pos/1000000,BAFsegm,pch=19,cex=0.5,col="green")
  dev.off()

  BAFphased = ifelse(BAFsegm>0.5,BAF,1-BAF)

  if(length(BAFphased)<50){
	BAFphseg = rep(mean(BAFphased),length(BAFphased))
  }else{
  	res= selectFastPcf(BAFphased,kmin,gamma*sdev,T)
  	BAFphseg = res$yhat
  }

  png(filename = paste(sample,"_",ifelse(grepl("chr", chr), "", "chr"),
                       chr,".png",sep=""), width = 2000, height = 1000, res = 200)
  par(mar = c(5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2, cex.lab = 2)
  plot(c(min(pos)/1000000,max(pos)/1000000),c(0,1),pch=".",type = "n",
	   main = paste(sample,", chromosome ", chr, sep=""), xlab = "Position (Mb)", ylab = "BAF (phased)")
  points(pos/1000000,BAF,pch=".",col=ifelse(BAFsegm>0.5,"red","blue"),cex=2)
  points(pos/1000000,BAFphseg,pch=19,cex=0.5,col="darkred")
  points(pos/1000000,1-BAFphseg,pch=19,cex=0.5,col="darkblue")
  dev.off()

  BAFphased = ifelse(BAFsegm>0.5,BAF,1-BAF)
  BAFoutputchr = cbind(rep(chr,length(BAFphseg)),pos,BAF,BAFphased,BAFphseg)
  BAFoutput=rbind(BAFoutput,BAFoutputchr)
}
colnames(BAFoutput) = c("Chromosome","Position","BAF","BAFphased","BAFseg")
write.table(BAFoutput,paste(sample,".BAFsegmented.txt",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)

q(save="no")
