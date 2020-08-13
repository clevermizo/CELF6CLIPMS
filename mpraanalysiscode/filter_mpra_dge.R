filter_mpra_dge=function(d,mpra.options){
  # Remove Divide-By-Zeros
  goodbarcodes = names(d$dna_counts)[d$dna_counts >= mpra.options$minCounts]
  # goodbarcodes are barcodes where there are DNA counts more than minCounts
  
  # Now, keep rows with minSamples with minCounts
  tmp = d$counts
  tmp = tmp[goodbarcodes,]
  # tmp is now d$counts with only barcodes with DNA counts > minCounts
  
  tmp.test = rowSums(tmp >= mpra.options$minCounts)
  goodbarcodes = names(tmp.test)[tmp.test >= mpra.options$minSamples]
  
  # goodbarcodes now are barcodes with at least: minSamples with >= minCounts,
  # and DNA counts (no NAs)
  
  tmp = d$element_name[goodbarcodes]
  # tmp is now all the elements represented by goodbarcodes
  tmp.table = table(tmp)
  # tmp.table has the barcode count for each element
  goodelements = names(tmp.table)[tmp.table >= mpra.options$minBarcodes]
  tmp2 = tmp[tmp %in% goodelements]
  # QC: on 7/1 with our MPRA data, this last step excludes 11 elements
  # Removes 16 barcodes
  goodbarcodes = names(tmp2)
  return(goodbarcodes)
}