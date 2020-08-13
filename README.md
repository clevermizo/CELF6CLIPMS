# CELF6CLIPMS
Code related to CELF6 CLIP Target identification &amp; analysis


Description:
Nearly all R code requires stringr & Biostrings packages. Also BSGenome is used to pull sequence data.

1. alwaysLoad.R <- a variety of preparation functions used throughout for "prettifying" tables, preparing statistics for reports, etc.
2. c6utrbccount.py <- Python code used to count barcodes in PTRE-Seq library sequencing FASTQ files
3. CallPeaks.bash <- BASH code used to call genome wide CLIP peaks and prepare for differential enrichment analysis
   CallPeaks.bash uses a gene-by-gene approach for background estimation & requires bedtools, samtools, & Piranha
4. makeChromosomeWindowsIntoFeaturesGTF.R <- used in CallPeaks.bash pipeline, converts windowed genome into GTF for stranded counting with subread featureCounts
5. makePeakFeatureGTFfromSigPeaks.R <- used in CallPeaks.bash pipeline, finds peaks meeting threshold criteria and prepares GTF for sample counting under peaks
6. cleanPiranhaWidths.R <- makes peakwidths called by Piranha uniform for counting under peaks
7. MEME-Enrichment-Processing.R <- calculates enriched CISBP-RNA database motif sequences and performs hierarchical clustering
8. create_PTRESEQ_mutantsequences.R <- mutates sequences for experimentation by analyzing position probability matrices for motifs (PPMs) ("PWM" is used everywhere in the code but strictly speaking these are PPMs)
9. MPRA-processing.R <- processing and filtering counted barcode reads in PTRE-Seq library experiments
10. MPRE-analysis.R <- statistical analysis of effects of condition & sequence mutation in PTRE-Seq library experiments

mpraanalysiscode has a number of functions
(requires lme4, multcomp, car, broom, edgeR R packages)

1. setup_mpra_dge <- use edgeR DGE object as wrapper for PTRE-Seq data (edgeR internal functions are not used for analysis, but edgeR provides a convenient wrapper for reading count data and filtering.
2. filter_mpra_dge <- filter by expression & # of surviving barcodes
3. mpra_modeling_functions.R <- A number of statistical analytic functions using linear mixed modeling (lme4) & ANOVA, as well as r2glmm implementation of the Nakagawa & Schielzeth estimates of R2 from LMMs/GLMMs. 

  reshape.mpra(): converts countdata to "long" format for repeated measures style analysis
  extractOmnibusStatistics(): computes ANOVA F tests for LMM, and estimates of R2, ICC, proportion of variance explained terms for fixed effects
  extractModelLogFCEstimates(): extract the log fold change estimates, SE, post-hoc pvals, and 95% confidence intervals for each condition between reference and mutated sequence
  extractAllTukeys(): extract all Tukey's posthoc pairwise comparisons with associated statistics
  extractDescriptiveStats(): extract means & ses & confidence intervals for each condition/sequence (the group means)
  collapseBarcodes(): calculates the barcode collapsed mean for each sample/condition/sequence in the data (the sample means)
  
