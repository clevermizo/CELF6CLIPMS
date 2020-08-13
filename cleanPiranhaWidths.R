#R code - adjust peak widths in Piranha output
library(stringr)
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}
# input variables from command line:
#  gtf: GTF of Piranha peaks meeting pval criterion
#  counts: tmpintersect.bed of gene's counts
#  peakwidth: desired truncation width

peakwidth=as.numeric(peakwidth)
#Read GTF
gtfdata=read.table(gtf,header=FALSE,sep="\t",as.is=TRUE)
colnames(gtfdata)=c("chr","V2","V3","start","end","V6","strand","V8","annot")

gtfdata.new=gtfdata
gtfdata.new$start=NA
gtfdata.new$end=NA

for(i in 1:nrow(gtfdata)){
	tmp.start=gtfdata$start[i]
	tmp.end=gtfdata$end[i]
	tmp.chr=gtfdata$chr[i]
	tmp.strand=gtfdata$strand[i]
	
  # use grep to quickly find the count lines that contain the peak boundaries
	searchStart = paste("grep -nP \"",
			    tmp.chr,"\t",tmp.start,"\t[0-9]*\t\\",tmp.strand,"\t\" ",counts," | cut -f1 -d\":\"",sep="")
	startLine=as.numeric(system(searchStart,intern=TRUE))

	searchEnd = paste("grep -nP \"",tmp.chr,"\t[0-9]*\t",tmp.end,"\t\\",tmp.strand,"\t\" ", counts," | cut -f1 -d\":\"",sep="")
	endLine=as.numeric(system(searchEnd,intern=TRUE))

	con=file(counts)
	open(con)
	
	tmp.lines=read.table(con,skip=startLine-1,nrow=(endLine-startLine+1))
	close(con)

	tmp.lines$midpoint=apply(cbind(tmp.lines$V2,tmp.lines$V3),1,median)
	k=which(tmp.lines$V5==max(tmp.lines$V5))[1] # It is possible that two bins in a row have the same count value, so pick the first corresponding to the maximum
	pkCenter=median(tmp.lines$midpoint[k])
	
	gtfdata.new$start[i]=pkCenter-peakwidth/2
	gtfdata.new$end[i]=pkCenter-1+peakwidth/2

	annot.string.tmp = gtfdata$annot[i]
	newpkname = paste("peakName ",tmp.chr,"_",gtfdata.new$start[i],"_",
			  gtfdata.new$end[i],"_",tmp.strand,sep="")

	annot.split=strsplit(annot.string.tmp,split=";",fixed=TRUE)[[1]]
	annot.split[1]=newpkname
		
	gtfdata.new$annot[i]=paste(annot.split,collapse=";")
}

newGTF=str_replace(gtf,".gtf",".widthclean.gtf")
write.table(gtfdata.new,
	    file=newGTF,
	    sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
