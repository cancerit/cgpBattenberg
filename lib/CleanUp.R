CleanUp<-function(chr,startOutput,startInput){
	cmd = paste("rm ",startOutput,chr,"_*_*.*",sep="")
	system(cmd)
	cmd = paste("rm ",startInput,chr,".txt",sep="")
	system(cmd)
}