setup_mpra_dge = function(targets,DNAcounts){
# Each target's count table is formatted as:
## barcode sequence sequenceName Count, the python counting code also outputs an unheadered column of rownames
## targets has minimally the column "files"
## DNA counts is formatted the same as targets
# read in counts
	require(edgeR)
	d = readDGE(targets,columns=c(1,4),
            header=TRUE,row.names=1,sep="\t")

	# Label count columns
	colnames(d$counts)=d$samples$sample

	# Read in 1 file's complete table as a normal data frame, read column 2 (barcodes) as rownames
	tmp = read.table(targets$files[1],
					 header=TRUE,sep="\t",row.names=2)
	dna = read.table(DNAcounts[1],
					 header=TRUE,sep="\t",row.names=2)
	
	
	d$element_name = tmp[rownames(d$counts),"seqName"] # Pull the element names
	 names(d$element_name) = rownames(d$counts)
	
	d$dna_counts = dna[rownames(d$counts),"bcCount"] # Pull dna counts
	 names(d$dna_counts) = rownames(d$counts)
	
	d$dna.cpm = as.vector(cpm(d$dna_counts))
	  names(d$dna.cpm)=names(d$dna_counts)
	
	d$dna.log = as.vector(cpm(d$dna_counts,log=TRUE))
	  names(d$dna.log)=names(d$dna_counts)
	
	d=calcNormFactors(d)
	d$cpm=cpm(d$counts,normalized.lib.sizes=TRUE)
	d$log=cpm(d$counts,normalized.lib.sizes=TRUE,log=TRUE)
	
	
	d$y = apply(d$cpm,2,function(x){x/d$dna.cpm})
	  rownames(d$y) = rownames(d$cpm)
	  badVals=is.infinite(d$y) | is.nan(d$y) # cleanout divide-by-zero errors
      d$y[badVals]=NA

	d$logy = apply(d$log,2,function(x){x-d$dna.log})
	  rownames(d$logy)=rownames(d$cpm) 
	  d$logy[badVals]=NA # cleanout divide-by-zero errors

	return(d)

}
