#!/bin/bash
#CallPeaks.bash
#### The following pipeline was used to call peaks across the genome, using by-gene based model of background counts ####
# Initialize #
wd=~/tmp/CELF6MS # project working directory
piranhapval="0.1" # cutoff Negative Binomial probability for peak maxima
peakwidth="100" # truncate peak bounds to this width, Piranha sometimes results in extremely wide boundaries that are several kb in size
bams=$wd/bam #location of bam files
countsdir=$wd/peakcalling/counts # Folder where counts across genomic windows are stored
peakfile=$wd/peakcalling/peaks/output.gtf # Final GTF file where rows are individual peaks. Used for Differential Expression, for counting reads across individual samples
chrData="mm10.genome" #tab-delimited file with chr sizes, first column is chr names, second is size in bp
# parameters for bedtools
windowSize="100" # window genome into windows of windowSize
stepSize="50"  # window overlap. Set to 0 for no overlap. 

# Pipeline steps:
# 1. Generate genomic windows
# 2. Generate merged BAM file
# 3. Generate a GTF for windows for stranded counting
# 4. Count reads strandedly in each window
# 5. Call peaks gene-by-gene
# Depends:
# bedtools, samtools, Piranha
# R code: 
#  makeChromosomeWindowsIntoFeatureGTF.R <- takes BED file from chromosome windowing & makes a GTF file that subread featureCounts can use
#  makePeakFeatureGTFfromSigPeaks.R <- takes BED file from Piranha's peak calling & makes a GTF file that subread featureCounts can use, subsetting to peaks with p less than piranhapval threshold above
#  cleanPiranhaWidths.R <- truncates peak features to a window symmetrical around the peak maximum, as determined with peakwidth above. Initially I didn't impose a restriction like this, but Piranha 
#  returned many peaks of unreasonable length that upon inspection were merges across wide regions, sometimes several kb. 

##################### 1. Generate Genomic Windows by Chromsome #############################################################
for i in $(cut -f1 $chrData)
do
	echo making windows for $i with width of $windowSize in steps of $stepSize nt
	grep -P "^$i\t" $chrData > tmp.chr
	bedtools makewindows -g tmp.chr -w $windowSize -s $stepSize > $i".windows"
done
rm tmp.chr

#################### 2. Generate merge of CLIP BAM files (Peaks are called on merged bam) ##################################
bamfiles=$(ls $bams/yfp-ip*.bam)
outputbam=fullmerge.bam
samtools merge $outputbam $bamfiles


################### 3. Convert BED files of chromosomal windows into stranded GTF for counting #############################
for i in $(cut -f1 $chrData)
do
	echo Making annotation GTF for bins for $i
	R CMD BATCH --no-save \
		--no-restore \
		"--args chrom=\"$i\"" \
		makeChromosomeWindowsIntoFeatureGTF.R
done
# Final GTFs have windows down the rows, with windows duplicated for both + & - strand. Reads are counted in stranded fashion. 

################### 4. Count reads in stranded fashion in each bin for each windows chromosome file###############################################

# featureCounts parameters:
# -t bin : counts by the feature column keyword "bin" in the GTFs, each is a window
# -f : count at the "feature" level (don't look for "genes" or metafeatures in the annotation column of the GTF)
# -O : assign reads to all overlapping bins. All windows overlap by 50% in this code. This double counts reads, but
#      the idea is to have  sliding window to look for peaks. Not concerned at the stage about double counting. 
# -s 2: reversely stranded (read 2 is sense, based on lib prep protocol)
# -p : count fragments/paired-end mode
# -B :only count reads with both mates properly aligned
# -a : the GTF with the bins (chromosome windows) from above
# -o : outputfile
# bamfile: here the merged bam above.

for i in $(cut -f1 $chrData)
do
	echo Beginning bin counting for $i
	gtf="$i.windows.gtf"
	outputCount="$countsdir/$i.count"
	echo Reading bins from $gtf
	echo Outputting counts to $outputCount
	featureCounts -t bin -f -O -s 2 -p -B\
		-a $gtf \
		-o $outputCount \
		$outputbam	
done
# Convert the output of featureCounts into a BED file
for i in $(cut -f1 $chrData)
do
	echo Creating $countpath/$i.bed for Piranha
	sed -e '1,2d' "$countpath/$i.count" | cut -f2-5,7 > "$countpath/$i.bed"
done



##################### 5. Call Peaks ####################################################################
# the input into this loop is "genes.gtf" which contains gene id annotation for mm10
# this works as follows:
# bedtools intersect finds the counts in the appropriately countfile (section 4) that overlap a gene
# 0s are truncated
# Piranha calls peaks along the bins overlapping the genes coordinates
# Output is subset based on p-value threshold(see beginning)
# Subsetted peaks are written to a final GTF, $peakfile (see beginning)
# This final GTF grows until all genes in the genes.gtf are processed
# The final GTF can then be used to count reads in individual sample & control bams for differential expression analysis by peak. 

# read parses the following gene-specific info needed later: chr, start, end, strand, comments
while read chr dataset feature start end null1 strand null2 comments 
do
 #Parse Gene Name for Readability
 genename=$(echo $comments | cut -f3 -d";" | sed "s/\ //g" | sed "s/gene_name//g")
 echo Processing peaks for $genename
 # Generate temporary BED & information file for the ith gene. This will be used with bedtools intersect to grab count data for that gene
 printf "%s\t%s\t%s\t%s\n" $chr $start $end $strand > $countsdir/tmplookup.bed # bed file representing the gene's coordinates
 echo $comments > $countsdir/tmpcomments.txt #gene id information for later, the GTF's 9th column
 lookupcounts=$countsdir/$chr.bed  # chr is the gene's chr 
 
 if [ -f $lookupcounts ] # there are counts on this chromosome
  then 
   # bedtools intersect with the -S flag reports that strandedness does not exist in the counts file. it does. Force strandedness with grep $strand in the counts file first. 
   grep -P "\t"$strand"\t" $lookupcounts > $countsdir"/tmpcounts.bed"
   
   # Intersect: -wa reports the features in the counts bed file since that's the 'a' file, 'b' file is tmplookup.bed: i.e., the gene's coords
   bedtools intersect -a $countsdir/tmpcounts.bed -b $countsdir/tmplookup.bed -wa | grep -v -P "\t0$" > $countsdir/tmpintersect.bed
  
   # tmpintersect.bed now are just the genomic window stranded counts for the windows/bins that overlap the gene of interest
   # grep -v -P "\t0$" removes 0 counts. If you don't remove 0s, Piranha calls a background way too low and too many peaks are called.
   # This is certainly an issue genome-wide, where most regions in CLIP don't have count data. But even across a single gene I've found this to be true. 
   
   if [ -s $countsdir/tmpintersect.bed ] # If there are counts at all run Piranha
   then
    echo Counts found, running Piranha
	# Piranha options: -o outputfile -s sort (just in case), -p 1.0 (return all ps, I will subset by threshold myself later), 
	# -c no correction, again we are using a threshold
	# -u 0 no peak merging. I find Piranha merges across very wide boundaries
	# -d ZeroTruncatedNegativeBinomial. Type of prob dist for p generation
    Piranha -o $countsdir/tmppeaks.bed -s -p 1.0 -c -u 0 -d ZeroTruncatedNegativeBinomial -v $countsdir/tmpintersect.bed
   else # there are no counts at all, make an empty file
    > $countsdir/tmppeaks.bed
	echo No counts, empty peaks file created.
   fi
  else # there are no counts for this chromosome, make an empty file
   > $countsdir/tmppeaks.bed
   echo No counts, empty peaks file created.
 fi
  
 # Now either an empty tmppeaks.bed file or a file with peaks has been created. If has data, then:
 if [ -s $countsdir/tmppeaks.bed ]
 then
  echo Using cutoff for Piranha based on $piranhapval
  R CMD BATCH --no-save \
	--no-restore \
	"--args gtfprefix=\"$countsdir/tmppeaks\" piranhacutoff=\"$piranhapval\" comments=\"$countsdir/tmpcomments.txt\"" \
	makePeakFeatureGTFfromSigPeaks.R
  # This creates $countsdir/tmppeaks.gtf. Each line is a peak, and the annotation column retains the gene name for access later.
  # tmppeaks.gtf shows peaks less than pvalue threshold (see beginning)
  newgtf="tmppeaks_p"$piranhapval".gtf"
  
  # Next, clean
  R CMD BATCH --no-save \
	--no-restore \
	"--args gtf=\"$countsdir/$newgtf\" counts=\"$countsdir/tmpintersect.bed\" peakwidth=\"$peakwidth\"" \
	cleanPiranhaWidths.R
  #outputfile is tmppeaks_p0.1.widthclean.gtf. If there is data here, add to growing output file:
  
   if [ -f $countsdir/"tmppeaks_p"$piranhapval".widthclean.gtf" ]
   then
    echo Adding $genename peaks to $peakfile
    cat $countsdir/"tmppeaks_p"$piranhapval".widthclean.gtf" >> $peakfile
   else
    echo 'No peaks at threshold ' $piranhapval 'for' $genename
   fi
 else #tmppeaks.bed is empty
  echo Doing nothing for $genename # there are no peaks, move on.
 fi
rm $countsdir/tmp* # remove temporary files. 
done < genes.gtf


























