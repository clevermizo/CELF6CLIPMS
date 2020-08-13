#R code
#Make GTF annotation file for genomic bins for counting respecting paired end
#fragment & strand alignment
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

peakFile=paste(gtfprefix,".bed",sep="")
gtfOut=paste(gtfprefix,"_","p",piranhacutoff,".gtf",sep="")
comments=readLines(comments)
if(file.size(peakFile) > 0){
  z=read.table(peakFile,header=FALSE,sep="\t",as.is=TRUE)
  z$padj = p.adjust(z$V7, "BH")
  # piranhacutoff can either be "FDR" (benjamini Hochberg) or a p-value
  if(piranhacutoff=="FDR"){
  zsig = z[z$padj < 0.05,]
  zsig = cbind(
  zsig[,1],
   c("called_peaks"),
   c("peak"),
   zsig[,2:3],
   c("."),
   zsig[,4],
   c("."),
   zsig[,c(5,7,8)])
   colnames(zsig)=
   c("chr","origin","feature","start","end","something","strand","something2",
   "pkheight","pval","padj")
  annotationColumn = paste(zsig$chr,zsig$start,zsig$end,zsig$strand,sep="_")
  annotationColumn = paste("peakName \"",annotationColumn,"\"; ",
                         "peakHeight \"",zsig[,"pkheight"],"\"; ",
                         "peakFDR \"",zsig[,"padj"],"\";",sep="")
  annotationColumn = paste(annotationColumn,comments,sep=" ")
  }else{
  piranhacutoff = as.numeric(piranhacutoff)
  zsig = z[z$V7 <= piranhacutoff,]
  zsig = cbind(
        zsig[,1],
        c("called_peaks"),
        c("peak"),
        zsig[,2:3],
        c("."),
        zsig[,4],
        c("."),
        zsig[,c(5,7)])
  colnames(zsig)=
   c("chr","origin","feature","start","end","something","strand","something2",
   "pkheight","pval")
  annotationColumn = paste(zsig$chr,zsig$start,zsig$end,zsig$strand,sep="_")
  annotationColumn = paste("peakName \"",annotationColumn,"\"; ",
                         "peakHeight \"",zsig[,"pkheight"],"\"; ",
                         "peakPvalue \"",zsig[,"pval"],"\";",sep="")
  annotationColumn = paste(annotationColumn,comments,sep=" ")

  }

  zsig = cbind(zsig[,1:8],annotationColumn)
  zsig_strand = zsig[zsig$strand %in% c("+","-"),]

  write.table(zsig_strand,
  file=gtfOut,
  sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

  quit(save="no",status=0)
}else {
  quit(save="no",status=1)
}
