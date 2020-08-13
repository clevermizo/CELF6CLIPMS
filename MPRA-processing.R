# MPRA-processing.R
# These code takes the sequence table used to generate the library and:
# -- prepares edgeR-style DGE object --
# -- filters data based on counts, samples, barcodes
# -- adds harmonized table of annotation information & rankcriteria for the targets
# FUNCTION CALLS:
# -- mpraanalysiscode/setup_mpra_dge.R <- setting up object
# -- mpraanalysiscode/filter_mpra_dge.R <- filtering PTRE seq data
# -- mpraanalysiscode/alwaysLoad.R <- various functions for output files
rm(list=ls())
counts="~/tmp/CELF6MS/mpra/counts" # barcode counts (output of python code for counting barcodes
newannotation="~/tmp/CELF6MS/mpra/MPRA_3utr_targets.rds" # peaks included in the PTRE-Seq library with their gene & coord information (See Supplemental Table S2)

### Libraries ###
source("~/Dropbox/alwaysLoad.R")
library(stringr)
library(edgeR)
library(lme4)
library(broom)
library(multcomp)
library(car)

### Options ###
# at least 10 counts in 4 samples, with at least 3 barcodes represented per sample.
# Set up targets files and groups

mpra.options=list(minCounts=10,
				  minBarcodes=3,
				  minSamples=4)

mpra.targets=data.frame(files=list.files("~/tmp/CELF6MS/mpra/counts/",pattern="[0-9].tab"))
# read sample names (e.g., c6-3, vector-1)
mpra.targets$sample = str_replace(mpra.targets$files,".tab","")
mpra.targets$sample2 = str_replace(mpra.targets$sample,"-input|-trap","")
# read group names (e.g., c6-input, vector-trap)
mpra.targets$group = str_replace(mpra.targets$files,".tab","")
mpra.targets$group = str_replace(mpra.targets$group,"-[0-9]","")
# append the count folder to filename for access
mpra.targets$files=paste0("~/tmp/CELF6MS/mpra/counts/",mpra.targets$files)


### Read DGE ###
source("~/tmp/CELF6MS/code/mpraanalysiscode/setup_mpra_dge.R")
d  = setup_mpra_dge(mpra.targets,"~/tmp/CELF6MS/mpra/counts/MPRA_library_DNA.tab")
# Celf6 Project Cleanup
d$element_name=str_replace(d$element_name," ","") # remove funky spaces
d$element_name=str_replace(d$element_name,"c6targctl_ref","c6targ_ctl") # relabel funky category 
   names(d$element_name) = rownames(d$count)
d$peak_name = str_split(d$element_name,"_c6targ_",simplify=TRUE)
d$allele = d$peak_name[,2]
  names(d$allele) = rownames(d$counts)
d$peak_name = d$peak_name[,1]
  names(d$peak_name) = rownames(d$counts)
# d$peak_name is the row index into seqtable above
# d$allele_name is the allele 
d$samples$sample = factor(d$samples$sample) # alphabetical is fine
d$samples$group = factor(as.character(d$samples$group),
						 c("vec-input",
						   "c6-input",
						   "c3-input",
						   "c4-input",
						   "c5-input",
						   "c3c6-input",
						   "c4c6-input",
						   "c5c6-input",
						   "vec-trap",
						   "c6-trap",
						   "c3-trap",
						   "c4-trap",
						   "c5-trap",
						   "c3c6-trap",
						   "c4c6-trap",
						   "c5c6-trap"))
						 
saveRDS(d,"~/tmp/CELF6MS/mpra/mpradge.rds")

######### FILTER DGE ##########
source("~/tmp/CELF6MS/code/mpraanalysiscode/filter_mpra_dge.R")
goodbarcodes = filter_mpra_dge(d,mpra.options)
# goodbarcodes; at least minSamples with minCounts, at least elements with minCounts in DNA, 
# at least elements with at least minBarcodes per element. 
dsub=d[goodbarcodes,]
dsub$element_name = dsub$element_name[goodbarcodes]
dsub$dna_counts=dsub$dna_counts[goodbarcodes]
dsub$dna.cpm=dsub$dna.cpm[goodbarcodes]
dsub$dna.log=dsub$dna.log[goodbarcodes]
dsub$cpm=dsub$cpm[goodbarcodes,]
dsub$log=dsub$log[goodbarcodes,]
dsub$y = dsub$y[goodbarcodes,]
dsub$logy = dsub$logy[goodbarcodes,]
dsub$peak_name = dsub$peak_name[goodbarcodes]
dsub$allele = dsub$allele[goodbarcodes]

####### Add Annotation #######
dsub$annotation=readRDS(newannotation)
#newannotation has the "harmonization" data with the rank criteria from the new expression analysis
#now, reduce the dsub$annotation to only elements in the filtered dge
dsub$annotation = dsub$annotation[intersect(rownames(dsub$annotation),unique(dsub$peak_name)),]
# Remove peaks that are no longer sig in the re-analysis. These can just be cleaned
# from the annotation table as we will loop over this for the actual statistical treatment.
dsub$annotation = dsub$annotation[!is.na(dsub$annotation$rankcriterion),]
# Remove peaks that are greater than 25 nucleotides (a "bin") away from the old called peak
# These are not the same peak
dsub$annotation = dsub$annotation[dsub$annotation$newdist2pk<=25,]
pkloc=paste(dsub$annotation$chr,dsub$annotation$strand,dsub$annotation$pkloc,sep=":")
names(pkloc)=rownames(dsub$annotation)
dsub$annotation = dsub$annotation[names(pkloc)[!duplicated(pkloc)],]
saveRDS(dsub,"~/tmp/CELF6MS/mpra/mpradge.filtered.rds")

