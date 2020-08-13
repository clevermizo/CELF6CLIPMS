# Assemble data for motif enrichment analysis

############################## Initialize #############################################
# Read in target data for 3'UTRs
library(stringr)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)

source("~/Dropbox/alwaysLoad.R")
targetfile='~/tmp/CELF6MS/allmethod_targetsandcontrols.rds'
#source("I:/home/clevermizo/Dropbox/alwaysLoad.R")
#targetfile='I:/home/clevermizo/tmp/CELF6MS/allmethod_targetsandcontrols.rds'

maxpk=10
sequencecols=c("logFC.inp","logFC.wt","chr","strand","pkloc","analysismethod","ensID","geneName","rankcriterion","geneFeature")
memelength=50
outputprefix="~/tmp/CELF6MS/motifenrichment/"
ameenrichmentthreshold=0.05
minmatrixsize=3 # At least three elements for clustering
dendrocutdist=1.8 #From previous analysis, Euclidean distance threshold for cutting
fimomatchp=0.005 # From previous analysis, threshold to call a FIMO match. 


############################## Prepare Data for Sequence Pulls ##############################
# logFC is CLIP enrichment vs. control samples, used to rank order sequences
# then need chromosome, strand, and peak locations. memelength is the size needed for analysis
# 50 bp total was used for the original analysis.

#### Filter out "whacky" numbers of peak calls. 
x=readRDS(targetfile)
x=x[x$npk<=maxpk,]
#### Expand based on peak ####
x = x[,sequencecols]
# listify columns for expansion
x$pkloc = str_split(x$pkloc,";")
# Sometimes "chr" & "strand" are delimited, and sometimes not: they got delimited when there were exon boundaries or alternate UTRs etc. in the
# original annotation.
# We will make the assumption that no elements span multiple chromosomes, or at least we will limit ourselves to those that
# are on single chromosomes or strands. 
x$chr=str_split(x$chr,";")
x$strand=str_split(x$strand,";")
chr.test=sapply(x$chr,function(y){return(length(unique(y)))})
strand.test=sapply(x$strand,function(y){return(length(unique(y)))})
# there was only 1 instance in which there was more than one chr listed, an intron in Gm1821,
# which is annotated at chr11 & chr14. Weird! Only 1 instance, so subset for chr.test==1 & strand.test==1
# This was called as "target" in the "input" rank criterion category. It has a biotype of "pseudogene"
# Apparently lives on multiple chromosomes!

x=x[chr.test==1 & strand.test==1,]
x$chr=sapply(x$chr,unique)
x$strand=sapply(x$strand,unique)
# Now that that is cleaned up, expand based on pk location coordinate
xx= as.data.frame(expand(DataFrame(x),c("pkloc")))
xx$pkloc = as.integer(xx$pkloc)
xx2=x
xx2$analysismethod=str_split(xx2$analysismethod,";")
xx2=as.data.frame(expand(DataFrame(xx2),c("analysismethod")))
xx$analysismethod=xx2$analysismethod


sequencelists=expand.grid(rankcriterion=c("both","input","wt"),geneFeature=c("utr3","utr5","intron","CDS"))
rownames(sequencelists) = paste(sequencelists$rankcriterion,sequencelists$geneFeature,sep="-")
# some instances of "cds" lowercase, most CDS.
xx$geneFeature =str_replace_all(xx$geneFeature,"cds","CDS")


############################## PULL SEQUENCES AND RUN AME ##############################
for(i in rownames(sequencelists)){
  print(paste0("Pulling sequence data for ",i))
  #Subset data
  tmp=xx[xx$rankcriterion == sequencelists[i,"rankcriterion"] &
         grepl(sequencelists[i,"geneFeature"],xx$geneFeature),]
  tmp=tmp[!is.na(tmp$pkloc),] # Remove any that fail to find a peak. Otherwise it pulls the whole chromosome
  #Rank Targets according to CLIP enrichment, decreasing
  tmp=tmp[order(-tmp$logFC.inp,-tmp$logFC.wt),]
  
  #Measure for subsampling controls
  N=nrow(tmp)
  
  tmp2=xx[xx$rankcriterion == "control" &
         grepl(sequencelists[i,"geneFeature"],xx$geneFeature),]
  tmp2=tmp2[!is.na(tmp2$pkloc),]
  N2=nrow(tmp2)
  
  #If Controls outnumber Targets,randomly subsample controls
  if(N2>N){
    tmp2=tmp2[
    sample(rownames(tmp2),size=N,replace=FALSE),
    ]
    }
  
  #If Targets outnumber Controls,take the top N2
  if(N>N2){
    tmp=tmp[1:N2,]
    }
  
  # If N==N2 do nothing
    
  
  outfile1=paste0(outputprefix,i,".fa")
  outfile2=paste0(outputprefix,i,"-control.fa")
  
  #Pull Target Sequences
  seq=RNAStringSet(getSeq(Mmusculus,names=tmp$chr,start=tmp$pkloc-(memelength/2-1),end=tmp$pkloc+memelength/2,
             strand=tmp$strand,as.character=FALSE))
  names(seq)=rownames(tmp)
  export(seq,outfile1,format="fasta")
  
  
  #Pull Control Sequences
  seq=RNAStringSet(getSeq(Mmusculus,names=tmp2$chr,start=tmp2$pkloc-(memelength/2-1),end=tmp2$pkloc+memelength/2,
             strand=tmp2$strand,as.character=FALSE))
  names(seq)=rownames(tmp2)
  export(seq,outfile2,format="fasta")
  
  
  outputdir=paste0(outputprefix,i)
  #AME makes an unnecessary folder with an HTML version of the output - cleanup and just move over the tab delimited file. 
  #E value is #motifs * pvalue, to return everything do E=#motifs. Here the Mus_musculus.meme has ~700 or so motifs. 
  # By the way, the website says there are 700 something, but when you look at the .meme file, there are only 93
  #Method uses ranksum/Mann Whitney on the max score. 
  amecmd=paste("ame --o",
                outputdir,
                "--oc",
                outputdir,
                "--evalue-report-threshold 800 ", 
                "--control",outfile2,
                "--method ranksum --scoring max",
                outfile1," ~/tmp/CELF6MS/code/motif_databases/CISBP-RNA/Mus_musculus.meme")
  outfile3=paste0(outputprefix,i,".ameresults.tsv")
  cpcmd=paste("cp",
             paste0(outputdir,"/ame.tsv"),
             outfile3)
  rmcmd=paste("rm -rf",outputdir)
  
  print("Running AME...")
  system(amecmd)
  print("Copying results...")
  system(cpcmd)
  print("Removing temporary folders...")
  system(rmcmd)
  }
    
############################## CLUSTER AME RESULTS ##############################
# Functions
pullPWM = function(startpoint,endpoint,dblines){
  motif.data = dblines[startpoint:endpoint]
  motif.data = motif.data[(grep("letter-probability\\smatrix",motif.data)+1):
                            (grep("URL",motif.data)-2)]
  motifConn = textConnection(motif.data)
  pwm.data = as.matrix(read.table(motifConn,header=FALSE,sep="\t",as.is=TRUE)[,1:4])
  close(motifConn)
  colnames(pwm.data)=c("A","C","G","T")
  return(pwm.data)
}
makePWMlines = function(pwm.data,motifname,altname){
  motif.data = c(paste("MOTIF",motifname,altname),"",
                 paste("letter-probability matrix: alength= 4 w= ",nrow(pwm.data)," nsites= 1 E= 0",sep=""))
  newPWMtext = apply(matrix(paste(as.character(pwm.data),"\t",sep=""),nrow=nrow(pwm.data),ncol=4),1,
                     function(x){str_c(x,collapse="")})
  motif.data = c(motif.data,newPWMtext,"","URL http://nourl.custom.com","")
  return(motif.data)
}
# Read in Mouse Motif PWM Data
cisbp.db = file("~/tmp/CELF6MS/code/motif_databases/CISBP-RNA/Mus_musculus.meme")
cisbp.dblines = readLines(cisbp.db)
close(cisbp.db)
motif.dbtable = data.frame(start = grep("^MOTIF",cisbp.dblines),stringsAsFactors = FALSE)
motif.dbtable$end = c((motif.dbtable$start[2:length(motif.dbtable$start)]-1),length(cisbp.dblines))
motif.dbheader = cisbp.dblines[1:(motif.dbtable$start[1]-1)]
rownames(motif.dbtable) = sapply(str_split(cisbp.dblines[motif.dbtable$start]," "),function(x){return(x[2])})
# motif.dbtable is each motifs PWM starting & ending lines in the cisbp.db file.
motif.dbtable$PWM=vector('list',length=nrow(motif.dbtable))
motif.dbtable$motiflength=NA
for(i in 1:nrow(motif.dbtable)){
  motif.dbtable$PWM[[i]] = pullPWM(motif.dbtable$start[i],motif.dbtable$end[i],cisbp.dblines)
  motif.dbtable$motiflength[i] = nrow(motif.dbtable$PWM[[i]])
}
names(motif.dbtable$PWM)=rownames(motif.dbtable)
for(i in rownames(sequencelists)){
  # Read in the dataset's ame scores and subset to ameenrichmentthreshold p value
  amedata=read.table(paste0(outputprefix,i,".ameresults.tsv"),header=TRUE,sep="\t")
  # Add rownames indices
  rownames(amedata)=amedata$motif_ID
  # Subset to threshold
  # Note, use the "pright" value, not the "two-tailed value" - we do not want to include the
  # hypothesis that a motif is enriched in the controls & not the targets, the hypothesis is directional
  amedata=subset(amedata,amedata$p.value <= ameenrichmentthreshold)
  if(nrow(amedata) >= minmatrixsize){
  print(paste0("Reading PWM for motifs enriched in the dataset ... ", i))
  pwmdata=motif.dbtable$PWM[rownames(amedata)]
  maxLength=max(motif.dbtable[rownames(amedata),"motiflength"])
  # Add 0s to pwmdata to normalize length
  for(j in 1:length(pwmdata)){
    zpad=rep(c(0,0,0,0),maxLength-nrow(pwmdata[[j]]))
    if(length(zpad)>0){
    zpad=matrix(zpad,ncol=4,nrow=maxLength-nrow(pwmdata[[j]]))
    pwmdata[[j]]=rbind(pwmdata[[j]],zpad)
    }
  }
  PWMvectorize=lapply(pwmdata,function(x){y=as.numeric(x)})
  #PWMvectorize elements can be re-matrixed as matrix(PWMvectorize[[elementname]],nrow=maxLength,ncol=4)
  PWMvectorize=do.call(rbind,PWMvectorize)
  print("...running clustering analysis....")
  #First, compute distances (uses default, Euclidean distance)
  PWMz=dist(PWMvectorize)
  #Hierarchical clustering
  PWMh=hclust(PWMz)
  #Make better labels
  #Str_replace code is ugly, but MOTIF_ID:AltNAME:Consensus is the output
  newlabels=paste(PWMh$labels,
                  str_replace_all(str_split(amedata[PWMh$labels,"motif_alt_ID"],"_",simplify=TRUE)[,1],"\\(|\\)",""),
                  amedata[PWMh$labels,"consensus"],
                  sep=":")
  PWMh.forplot=PWMh
  PWMh.forplot$labels=newlabels
  PWMdnd=as.dendrogram(PWMh.forplot)
  #Make Dendrogram
  setEPS()
  postscript(paste0(outputprefix,i,".dendrogram.eps"))
  plot(PWMdnd)
  dev.off()
    
  # Cut the tree at dendrocutdist
  clusters=cutree(PWMh,h=dendrocutdist)
  amedata$cluster=clusters[rownames(amedata)]
  amedata=amedata[order(amedata$cluster),]
  
  # Generate Average PWMs based on clusters
    
  clusterPWMs=vector('list',max(amedata$cluster))
  clusterDB=motif.dbheader
  clusterAltNames=rep(NA,max(amedata$cluster))
  clusterMotifNames=paste0(i,"_cluster_",unique(amedata$cluster))
  
  for(j in 1:length(clusterPWMs)){
    cluster.motifs = PWMvectorize[amedata[amedata$cluster==j,"motif_ID"],]
    clusterAltNames[j]=str_c(str_replace_all(str_split(amedata[amedata$cluster==j,"motif_alt_ID"],"_",simplify=TRUE)[,1],"\\(|\\)",""),collapse=",")
    if(!is.null(dim(cluster.motifs))){
    cluster.avg = colMeans(cluster.motifs)
    }else{
      cluster.avg=cluster.motifs
    }
    cluster.avg = matrix(cluster.avg,nrow=maxLength,ncol=4)
    #Remove all 0 rows (pads from earlier) - These are always at the end.
    cluster.avg = cluster.avg[rowSums(cluster.avg==0)<4,]
    #Renormalize to 1.0 per row
    clusterPWMs[[j]]=cluster.avg/matrix(rep(rowSums(cluster.avg),4),nrow=nrow(cluster.avg),ncol=4)
    clusterDB=c(clusterDB,makePWMlines(clusterPWMs[[j]],motifname=clusterMotifNames[j],altname=clusterAltNames[j]))
    }
  print("...writing output...")
  clusterDB.file = file(paste0(outputprefix,i,".clusters.meme"))
  writeLines(clusterDB,clusterDB.file)
  close(clusterDB.file) 
  }else{
    print('Not enough data to perform clustering analysis')
  }
}

############################## Annotate Targets For Presence of Clusters ##############################
for(i in rownames(sequencelists)){
  print(paste0("Beginning Annotation for AME Clusters on peak sequences for the dataset ",i))
  #Subset data
  tmp=xx[xx$rankcriterion == sequencelists[i,"rankcriterion"] &
         grepl(sequencelists[i,"geneFeature"],xx$geneFeature),]
  tmp=tmp[!is.na(tmp$pkloc),] # Remove any that fail to find a peak. Otherwise it pulls the whole chromosome
  #Rank Targets according to CLIP enrichment, decreasing
  tmp=tmp[order(-tmp$logFC.inp,-tmp$logFC.wt),]
  #Run FIMO against the appropriate cluster meme file
  clusterDB.file=paste0(outputprefix,i,".clusters.meme")
  seq.file=paste0(outputprefix,i,".fa")
  
  if(file.exists(clusterDB.file)){
  print("....running match search on cluster averages")
  FIMOcmd=paste(
    "fimo --oc",paste0(outputprefix,i,"-fimo"),"--verbosity 4 --thresh 1 --norc",
    clusterDB.file,seq.file
    )
  system(FIMOcmd)
  print("....copying results")
  system(paste(
    "cp",
    paste0(outputprefix,i,"-fimo","/fimo.tsv"),
    paste0(outputprefix,i,".clusters.fimo")
    )
  )
  print("....removing temporarily files")
  system(paste(
    "rm -rf",
     paste0(outputprefix,i,"-fimo")))
  
  FIMOdata=read.table(paste0(outputprefix,i,".clusters.fimo"),header=TRUE,sep="\t")
  
  clusterMatch=matrix(NA,nrow=nrow(tmp),ncol=length(unique(FIMOdata$motif_id)))
  colnames(clusterMatch)=unique(FIMOdata$motif_id)
  rownames(clusterMatch)=rownames(tmp)
  print("....analyzing presence or absence of matches for each cluster ....")
  for(j in rownames(clusterMatch)){
    #IF sequence made it in
    if(sum(FIMOdata$sequence_name==j)>0){
    tmp2=FIMOdata[FIMOdata$p.value <= fimomatchp & FIMOdata$sequence_name==j,]
    tmp2=table(tmp2$motif_id)
    clusterMatch[j,]=0
    clusterMatch[j,names(tmp2)]=tmp2
    }
    }
  
  #Sort columns by "top ranking" clusters, the clusters containing the most representation across the sequences.
  clustertable=FIMOdata[,c("motif_id","motif_alt_id")]
  clustertable=clustertable[!duplicated(clustertable),]
  clustertable$rank=colSums(clusterMatch>0,na.rm=TRUE)
  clustertable=clustertable[order(clustertable$rank,decreasing=TRUE),]
  
  tmp=cbind(tmp,as.data.frame(clusterMatch[rownames(tmp),clustertable$motif_id]))
  print('.... writing output annotation & relative cluster rankings ....')
  write.table(tmp,paste0(outputprefix,i,".clusterannotation.tab"),sep="\t",row.names=TRUE,col.names=NA,quote=TRUE)
  write.table(clustertable,paste0(outputprefix,i,".clusterranking.tab"),sep="\t",row.names=FALSE,col.names=TRUE,quote=TRUE)
  }else{
  print("...not enough data to perform cluster analsyis, skipping.")
 }
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
