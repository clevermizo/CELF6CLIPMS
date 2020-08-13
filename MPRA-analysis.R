# MPRA-processing.R
# These code takes the edited and annotated DGE and runs the stats and generates figures.
# -- prep omnibus stats
# -- prep table of log fold change 
# -- prep individual figures for each element
# -- redraws figures for paper
# FUNCTION CALLS:
# -- mpra_modeling_functions.R
# -- 

### Libraries ###
rm(list=ls())
source("~/Dropbox/alwaysLoad.R")
library(stringr)
library(edgeR)
library(lme4)
library(broom)
library(multcomp)
library(car)
library(r2glmm)
dsub=readRDS("~/tmp/CELF6MS/mpra/mpradge.filtered.rds")
#################### Analysis Object #################################################
danalysis=dsub$annotation[,c("logFC.inp","logFC.wt","p.inp","p.wt",
                             "chr","strand","pkloc","newpkloc",
                             "geneName","ensID","HD","rankcriterion","criteriapass")]
# New columns:
# Omnibus effects
danalysis[,c("p_sequence","p_condition","p_interaction","r2","icc","pve_sequence","pve_condition","pve_interaction")]=
  NA
# Log Fold Change for Graphs
danalysis$mpra.logFC=vector("list",length=nrow(danalysis))
  names(danalysis$mpra.logFC)=rownames(danalysis)
# Model Means
danalysis$mpra.means=vector("list",length=nrow(danalysis))
  names(danalysis$mpra.means)=rownames(danalysis)
# Barcode Collapsed Data for Individual Plots
danalysis$mpra.collapse=vector("list",length=nrow(danalysis))
   names(danalysis$mpra.collapse)=rownames(danalysis)

   
#####################Split input & trap & calculate TE##############################   
dsub$samples$condition=str_replace(as.character(dsub$samples$group),"-input","")
dsub$samples$condition=str_replace(as.character(dsub$samples$condition),"-trap","")
dsub$samples$fraction=str_split(as.character(dsub$samples$group),"-",simplify = TRUE)[,2]

dsub.input = list(samples=dsub$samples,logy=dsub$logy,peak_name=dsub$peak_name,allele=dsub$allele)
dsub.input$samples = dsub.input$samples[dsub.input$samples$fraction=="input",]
dsub.input$logy = dsub.input$logy[,as.character(dsub.input$samples$sample)]
colnames(dsub.input$logy)=dsub.input$samples$sample2

dsub.trap = list(samples=dsub$samples,logy=dsub$logy,peak_name=dsub$peak_name,allele=dsub$allele)
dsub.trap$samples = dsub.trap$samples[dsub.trap$samples$fraction=="trap",]
dsub.trap$logy = dsub.trap$logy[,as.character(dsub.trap$samples$sample)]
colnames(dsub.trap$logy)=dsub.trap$samples$sample2

dsub.te = dsub.trap
# make sure the columns are in the same order for both matrices before calculation:
tmp=dsub.input$logy[,colnames(dsub.te$logy)]
tmp=dsub.te$logy - tmp
dsub.te$logy = tmp
rm(tmp)
# Load LMM calculation function
# Assign groupLevels for output models & graphs
groupLevels=c("vec","c6","c3","c4","c5","c3c6","c4c6","c5c6")


################## Load Functions & Process Elements ##########################
source("~/tmp/CELF6MS/code/mpraanalysiscode/mpra_modeling_functions.R")
danalysis=list(input=danalysis,trap=danalysis,te=danalysis)
dsubs=list(input=dsub.input,trap=dsub.trap,te=dsub.te)
rm(dsub.input,dsub.trap,dsub.te)
# Main loop, process Input, TRAP, and TE
for(ds in names(dsubs)){
 for(i in rownames(danalysis[[ds]])){
  bcs=names(dsub$peak_name)[dsub$peak_name==i]
  print(paste0("Reshaping data for modeling ",i, " in ",ds))
  df.long=reshape.mpra(d=dsubs[[ds]],bcs=bcs,groupLevels=groupLevels,minBarcodes=3)
  if(is.data.frame(df.long)){
    print(paste0("...Calculating LMM"))
     mdl=lmer(logy~allele*condition+(1|element),df.long)
    
    print(paste0("...extracting ANOVA"))
     omni=extractOmnibusStatistics(mdl)
     
    print(paste0("...adding omni to danalysis"))
     danalysis[[ds]][i,names(omni)]=omni
    
    print(paste0("...extracting log FC reference - mutant"))
     logfc=extractModelLogFCEstimates(mdl,df.long)
     danalysis[[ds]]$mpra.logFC[[i]]=logfc
    
    print(paste0("...extracting descriptive statistics (means, ses, confints), for graphs"))
     desc.stats=extractDescriptiveStats(mdl,df.long)
     danalysis[[ds]]$mpra.means[[i]] = desc.stats
     
    print(paste0("...extracting the barcode collapsed sample data"))
     ag=collapseBarcodes(df.long)
     danalysis[[ds]]$mpra.collapse[[i]]=ag
  }else{
    print("....data does not have enough barcode representation per sequence element or missing levels,skipping")
  }
  print(paste0("Finished with ",i," in ",ds))
 }
}
saveRDS(danalysis,file="~/tmp/CELF6MS/mpra/mpradge.analysis.rds")