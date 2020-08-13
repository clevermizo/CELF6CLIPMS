################ This code was used to generate motif mutations used in the PTRE-Seq Library ##################################
library(stringr)
library(Biostrings)
#####Read AME output#####
ame.data = read.table("./AME/ame_50nt_maxscore_ranksum.tab",header=TRUE,sep="\t",as.is=TRUE)
ame.data$motif.name = make.unique(sapply(str_split(ame.data$motif_alt,"_"),function(x){return(x[1])}))
rownames(ame.data)=ame.data$motif
ame.motifs = file("./FIMO/ame50nt.list")
writeLines(ame.data$motif,ame.motifs)
close(ame.motifs)

#####Read cloning sequences and MEME sequences#####
cloning.sequences = read.table("3utr_targets.cloning.tab",as.is=TRUE,sep="\t",header=TRUE,row.names = 1)
meme.sequences = read.table("3utr_targets.memeanalysis.tab",as.is=TRUE,sep="\t",header=TRUE,row.names = 1)

#####Read lines from mus musculues MEME file#####
cisbp.db = file("~/jdlab/Mike/MEME_Motifs/motif_databases/CISBP-RNA/Mus_musculus.meme")
cisbp.dblines = readLines(cisbp.db)
close(cisbp.db)

#####Make a small CISBP for running FIMO#####
motif.dbtable = data.frame(start = grep("^MOTIF",cisbp.dblines))
motif.dbtable$end = c((motif.dbtable$start[2:length(motif.dbtable$start)]-1),length(cisbp.dblines))
motif.dbheader = cisbp.dblines[1:(motif.dbtable$start[1]-1)]
rownames(motif.dbtable) = sapply(str_split(cisbp.dblines[motif.dbtable$start]," "),function(x){return(x[2])})
motif.dbtable.sub = motif.dbtable[ame.data$motif,]
# Make a small database
small.cisbpdb = motif.dbheader
for(i in 1:nrow(motif.dbtable.sub)){
  linesToPull = as.integer(motif.dbtable.sub[i,"start"]:motif.dbtable.sub[i,"end"])
  small.cisbpdb = c(small.cisbpdb,cisbp.dblines[linesToPull])
}
# Write output
small.cisbpfile = file("FIMO/ameout.meme")
writeLines(small.cisbpdb,small.cisbpfile)
close(small.cisbpfile)


##### Run FIMO and read output #####
system(command = "echo Running FIMO search of AME motifs")
system(command = "/usr/local/bin/meme/bin/fimo --oc ./FIMO/output --verbosity 4 --thresh 1 --norc ./FIMO/ameout.meme ./3utr_targets.memeanalysis.fa")
system(command = "cp ./FIMO/output/fimo.txt ./FIMO/ameout.fimo")
system(command = "sed 's/#//g' ./FIMO/ameout.fimo > ./FIMO/ameout.fimo.tab")
fimo.data = read.table("FIMO/ameout.tab",header=TRUE,sep="\t",as.is=TRUE)
fimo.thresh = 5E-3 # threshold for a sig match
pwm.thresh = 0.8 # threshold for mutating 

meme.sequences$ame.top.match = ""
meme.sequences$ame.top.alt = ""
meme.sequences$ame.top.alt.short = ""
meme.sequences$ame.top.start = NA
meme.sequences$ame.top.end = NA
meme.sequences$ame.adjpval = NA
meme.sequences$fimo.score = NA
meme.sequences$fimo.pval = NA
meme.sequences$fimo.context = ""
meme.sequences$matched.seqs = ""
meme.sequences$mut.seqs=""
meme.sequences$mut.context =""
##### Explore FIMO results by sequence #####
for(i in rownames(meme.sequences)){
  k=which(fimo.data$sequence.name==i)
  tmp=fimo.data[k,]
  tmp=tmp[tmp$p.value<fimo.thresh,]
  tmp=tmp[order(tmp$start,decreasing=FALSE),]
  if(nrow(tmp)>0){
  # Consider sites < fimo.thresh
  rownames(tmp)=NULL
  tmp$context = as.character(meme.sequences[i,"sequences"])
  context.tmp = as.character(meme.sequences[i,"sequences"])
  context.output=context.tmp
  for(j in 1:nrow(tmp)){
    tmp$matched.sequence[j] = str_replace_all(tmp$matched.sequence[j],"U","T")
    tmp$context[j] = paste(
      substr(tmp$context[j],1,tmp$start[j]-1),
      str_to_lower(substr(tmp$context[j],tmp$start[j],tmp$stop[j])),
      substr(tmp$context[j],tmp$stop[j]+1,str_length(tmp$context[j])),sep="")
    context.output = paste(
          substr(context.output,1,tmp$start[j]-1),
          str_to_lower(substr(context.output,tmp$start[j],tmp$stop[j])),
          substr(context.output,tmp$stop[j]+1,str_length(context.output)),sep="")
  }
  meme.sequences[i,"ame.top.match"] = str_c(tmp$pattern.name,collapse=",")
  meme.sequences[i,"ame.top.alt"] = str_c(as.character(ame.data[tmp$pattern.name,"motif_alt"]),collapse=",")
  meme.sequences[i,"ame.top.alt.short"] = str_c(as.character(ame.data[tmp$pattern.name,"motif.name"]),collapse=",")
  meme.sequences[i,"ame.top.start"] = str_c(tmp$start,collapse=",")
  meme.sequences[i,"ame.top.end"] = str_c(tmp$stop,collapse=",")
  meme.sequences[i,"ame.adjpval"] = str_c(as.numeric(ame.data[tmp$pattern.name,"adj_two.tailed"]),collapse=",")
  meme.sequences[i,"fimo.score"] = str_c(tmp$score,collapse=",")
  meme.sequences[i,"fimo.pval"] = str_c(tmp$p.value,collapse=",")
  meme.sequences[i,"fimo.context"] = context.output
  meme.sequences[i,"matched.seqs"] = str_c(tmp$matched.sequence,collapse=",")
  # Determine Mutations to Make
  tmp$matched.mutant = ""
  tmp$context.mutant = ""
  for(j in 1:nrow(tmp)){
   motif.data = cisbp.dblines[motif.dbtable.sub[tmp$pattern.name[j],"start"]:motif.dbtable.sub[tmp$pattern.name[j],"end"]]
   motif.data = motif.data[(grep("letter-probability\\smatrix",motif.data)+1):
                            (grep("URL",motif.data)-2)]
   motifConn = textConnection(motif.data)
   pwm.data = as.matrix(read.table(motifConn,header=FALSE,sep="\t",as.is=TRUE)[,1:4])
   close(motifConn)
   colnames(pwm.data)=c("A","C","G","T")
   pwm.test = which(rowSums(pwm.data >= pwm.thresh)!=0)
   mut.tmp = str_split(str_to_lower(tmp$matched.sequence[j]),pattern="")[[1]]
   for(p in pwm.test){
     if(str_to_upper(mut.tmp[p])==colnames(pwm.data)[as.numeric(pwm.data[p,])>=pwm.thresh]){
       mut.tmp[p]=colnames(pwm.data)[which.min(as.numeric(pwm.data[p,]))]
     }
   }
   tmp$matched.mutant[j] = str_c(mut.tmp,collapse="")
   tmp$context.mutant[j] = paste(
     substr(tmp$context[j],1,tmp$start[j]-1),
     tmp$matched.mutant[j],
     substr(tmp$context[j],tmp$stop[j]+1,str_length(tmp$context[j])),sep="")
  }
  # Make a full length mutant, respecting hierarchy of FIMO p-values and overlap
  # keep mutant sequence from previous in list, basically, to concatenate a full mutant when motifs overlap:
  context.mutant=tmp$context.mutant[1]
  if(nrow(tmp)>1){
      context.mutant = data.frame(start=1,stop=tmp$stop[1])
  for(j in 2:nrow(tmp)){
    if(j<nrow(tmp)){
      context.mutant = rbind(context.mutant,
                             data.frame(start=tmp$stop[j-1]+1,
                                        stop=tmp$stop[j]))
    }else{
      context.mutant = rbind(context.mutant,
                             data.frame(start=tmp$stop[j-1]+1,
                                        stop=str_length(tmp$context.mutant[j])))
    }
  }
      context.mutant$keep = (context.mutant$stop - context.mutant$start)>=0
      tmp2 = tmp[context.mutant$keep,]
      context.mutant = context.mutant[context.mutant$keep,]
      context.mutant$start = c(1,context.mutant$stop[1:(nrow(context.mutant)-1)]+1)
      rownames(context.mutant)=NULL
      rownames(tmp2)=NULL
      #Progressive cleaning:
      while(sum((context.mutant$stop - context.mutant$start)<0)>0){
      context.mutant$keep = (context.mutant$stop - context.mutant$start)>=0
      tmp2 = tmp2[context.mutant$keep,]
      context.mutant = context.mutant[context.mutant$keep,]
      context.mutant$start = c(1,context.mutant$stop[1:(nrow(context.mutant)-1)]+1)
      rownames(context.mutant)=NULL
      rownames(tmp2)=NULL
      }
      context.mutant$seq=""
  for(j in 1:nrow(tmp2)){
      context.mutant$seq[j] = substr(tmp2$context.mutant[j],context.mutant$start[j],context.mutant$stop[j])
  }
      context.mutant = str_c(context.mutant$seq,collapse="")
  }
  
  meme.sequences[i,"mut.seqs"] = str_c(tmp$matched.mutant,collapse=",")
  meme.sequences[i,"mut.context"] = context.mutant
  }
}


##### Remove sequences without sig matches to motifs####
meme.sequences.sub = meme.sequences[meme.sequences$mut.context!="",]
write.table(meme.sequences.sub,file="3utr_targets.memeanalysisWithMutations.tab",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)

##### Make Cloning Sequence Set ####
cloning.sequences = cloning.sequences[rownames(meme.sequences.sub),]
cloning.sequences$mutant=""
tackOn = (cloning.sequences$retrievedLen[1] - meme.sequences$retrievedLen[1])/2
memeLen = meme.sequences$retrievedLen[1]
for(i in rownames(cloning.sequences)){
  cloning.sequences[i,"mutant"] = str_to_upper(paste(substr(cloning.sequences[i,"sequences"],1,tackOn),
                                        meme.sequences.sub[i,"mut.context"],
                                        substr(cloning.sequences[i,"sequences"],tackOn+memeLen+1,str_length(cloning.sequences[i,"sequences"])),
                                        sep=""))
}

wt.sequences = DNAStringSet(cloning.sequences$sequences)
names(wt.sequences) = rownames(cloning.sequences)
mutant.sequences = DNAStringSet(cloning.sequences$mutant)
names(mutant.sequences) = paste(rownames(cloning.sequences),"_mutant",sep="")
export(wt.sequences,"3utr_final_cloning.fa",format="fasta")
export(mutant.sequences,"3utr_final_mutants.fa",format="fasta")


                                 