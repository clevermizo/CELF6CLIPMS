#make_replicate_scatterplots.R
if(!file.exists(paste0("replicatescatterplots/",
                       configuration$experiment_name,
                       "pairwisesamplecorrelations.tab"))){
if(sum(ls()=="d")==0){
  if(!file.exists(paste0(configuration$experiment_name,"_DGE.Rdata",sep=""))){
    print("Go back and run setup_ptreseq_dge.R before proceeding. DGE is missing")
  }else{
    print("DGE is not yet loaded. Loading...")
    load(paste0(configuration$experiment_name,"_DGE.Rdata",sep=""))
    print("DGE is now loaded. Proceeding with replicate scatter plot generation.")
  }
}else{
  print("DGE is already loaded. Proceeding with replicate scatter plot generation.")
}

C=cor(d$log)
CSummary = corSummaryTable(C)
CSummary$sample1class=d$samples$fraction[d$samples$sample %in% CSummary$sample1]
CSummary$sample2class=d$samples$fraction[d$samples$sample %in% CSummary$sample2]
CSummary$comparison=paste(CSummary$sample1class,CSummary$sample2class,sep="X")
if(!file.exists("replicatescatterplots")){
  dir.create("replicatescatterplots")
}
write.table(CSummary,
            file=paste0("replicatescatterplots/",
                        configuration$experiment_name,
                        "pairwisesamplecorrelations.tab"),
                    sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

for(i in 1:nrow(CSummary)){
print(paste0("Printing scatter plot ",
             CSummary$sample1[i],"X",CSummary$sample2[i]))
setEPS()
postscript(file=paste0(
  "replicatescatterplots/",
  configuration$experiment_name,
  "_",
  CSummary$sample1[i],
  "_",
  CSummary$sample2[i],
  ".eps"))
plot(d$log[,as.character(CSummary[i,c(1:2)])],pch=20,
     main=round(CSummary[i,3],4),
     ylim=c(floor(min(d$log)),ceiling(max(d$log))),
     xlim=c(floor(min(d$log)),ceiling(max(d$log))),
     xaxt="none",
     yaxt="none")
axis(side = 1,at = seq(floor(min(d$log)),ceiling(max(d$log)),by=1))
axis(side = 2,at = seq(floor(min(d$log)),ceiling(max(d$log)),by=1))
dev.off()
}
}else{
  print("It looks like replicate scatterplots have already been generated. Check under replicatescatterplots/")
  print("If necessary, delete this folder's contents and re-run, or run source(make_replicate_scatterplots.R) manually")
}
