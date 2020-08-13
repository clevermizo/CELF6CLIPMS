#R code
#Make GTF annotation file for genomic bins for counting respecting paired end
#fragment & strand alignment
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

windowFile=paste(chrom,".windows",sep="")
z=read.table(windowFile,header=FALSE,sep="\t",as.is=TRUE)
z2=cbind(z[,1],c("binned_data"),c("bin"),z[,2:3],c("."),c("+"),
	c("."))
colnames(z2)=c("chr","origin","feature","start","end","score","strand","frame")
z3=z2
z3$strand="-"
z = rbind(z2,z3)

annotationColumn = paste(z$chr,z$start,z$end,z$strand,sep="_")
annotationColumn = paste("binName \"",annotationColumn,"\"",sep="")

z=cbind(z,annotationColumn)

write.table(z,file=paste(windowFile,".gtf",sep=""),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
