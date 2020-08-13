####################################################################################################
# reshape.mpra makes a long-format data frame for LMM computation
####################################################################################################
reshape.mpra = function(d,bcs,groupLevels,minBarcodes=3){
  table(d$alleles[bcs])
  df=as.data.frame(t(d$logy[bcs,])) # transpose logy, rows are samples, cols are barcodes
  colnames(df)=paste0("logy:",colnames(df))
  # Add Condition
  df$condition=d$samples$condition
  df$condition=factor(df$condition,levels=groupLevels,labels=groupLevels)
  # ( Add allele after converting to "longform")
  df.long=reshape(df, # the data
                  direction="long", # direction for reshaping
                  varying=colnames(df)[grepl("logy",colnames(df))], #repeated measure
                  sep=":", #repeated measure delimiter
                  ids=rownames(df), #sample names
                  idvar="sample", #name for sample variable
                  timevar="barcode") #name for repeated measure variable
  df.long$allele=d$allele[df.long$barcode] #add alleles
  # Additional Factors for modeling
  df.long$sample=factor(df.long$sample) #alphabetical ordering for sample (doesn't matter)
  df.long$allele=factor(df.long$allele,levels=c("ref","alt"))
  df.long$element=factor(paste0(df.long$sample,"X",df.long$allele))
  
  # Add "_" character to column names for contrasts matrix for effects factors to 
  # ensure they are readable in model output
  colnames(contrasts(df.long$condition))=paste0("_",colnames(contrasts(df.long$condition)))
  colnames(contrasts(df.long$allele))=paste0("_",colnames(contrasts(df.long$allele)))
  colnames(contrasts(df.long$sample))=paste0("_",colnames(contrasts(df.long$sample)))
  colnames(contrasts(df.long$element))=paste0("_",colnames(contrasts(df.long$element)))
  
  # Require all sampleXallele ("element") to have minBarcodes to proceed with modeling and require 
  # allele & condition to have the necessary levels
  if(sum(table(df.long$element) >= minBarcodes)==length(table(df.long$element))){
    if( length(unique(df.long$allele)) == length(levels(df.long$allele)) &
        length(unique(df.long$condition)) == length(levels(df.long$condition))){
    return(df.long)
      }else{
      print(paste0("There are missing factor levels, returning NA"))
      return(NA)
      }
  }else{
    print(paste0("Minimum barcodes per element (=",minBarcodes,") not met, returning NA"))
    return(NA)
  }
}

####################################################################################################
# extractOmnibusStatistics pulls ANOVA, and estimates R2, ICC, and partial variance explained terms
####################################################################################################
extractOmnibusStatistics = function(mdl){
  require(car)
  require(r2glmm)
  
  #Omnibus p values
  a = as.data.frame(Anova(mdl,test.statistic = "F"))
  a = a$`Pr(>F)`
  names(a)=c("p_sequence","p_condition","p_interaction")
  
  #Omnibus R2
  r2=r2beta(mdl,partial = FALSE,method="nsj")
  
  #Omnibus ICC
  icc=as.data.frame(VarCorr(mdl))$vcov[1]/sum(as.data.frame(VarCorr(mdl))$vcov)
  
  #pve
  # submodels
  mdl_nointeraction = update(mdl,.~.-allele:condition)
  mdl_nocondition = update(mdl,.~.-condition-allele:condition)
  mdl_noallele = update(mdl,.~.-allele-allele:condition)
  
  # submodels Nakagawa-Schielzeth R2
  r2_nointeraction = r2beta(mdl_nointeraction,partial=FALSE,method="nsj")
  r2_nocondition = r2beta(mdl_nocondition,partial=FALSE,method="nsj")
  r2_noallele = r2beta(mdl_noallele,partial=FALSE,method="nsj")
  
  pve=c(pve_sequence=r2_nointeraction$Rsq - r2_noallele$Rsq,
        pve_condition=r2_nointeraction$Rsq - r2_nocondition$Rsq,
        pve_interaction=r2$Rsq - r2_nointeraction$Rsq)
  
  omnistats=c(a,r2=r2$Rsq,icc=icc,pve)
  
  return(omnistats)
}

########################################################################################################
# extractModelLogFCEstimates pulls the estimates for log FC reference - mutant for each elementXcondition
##########################################################################################################

extractModelLogFCEstimates = function(mdl,df.long){
 require(multcomp)
 require(broom)
 
 # Make an exemplar dataframe with conditions & ref/alt
 exemplar=expand.grid(allele=levels(df.long$allele),condition=levels(df.long$condition))
  contrasts(exemplar$allele)=contrasts(df.long$allele)
  contrasts(exemplar$condition)=contrasts(df.long$condition)
 
 # Make the model matrix to transform the betas
  X=model.matrix(~allele*condition,exemplar)
  rownames(X)=paste(exemplar$condition,exemplar$allele,sep=":")
  
 # Computer the linear terms for the fold change REF minus ALT
  Xref=X[exemplar$allele=="ref",]
  Xalt=X[exemplar$allele=="alt",]
  
  linfct=Xref-Xalt
  rownames(linfct)=paste0(rownames(linfct),"-alt")

 # compute the values using multcomp's glht
  
  logfc=glht(mdl,linfct=linfct)
  pvals=tidy(summary(logfc))
  ci=tidy(confint(logfc))
  
  logfc.summary = data.frame(logfc=pvals$estimate,
                             se=pvals$std.error,
                             z=pvals$statistic,
                             p=pvals$p.value,
                             lo=ci$conf.low,
                             hi=ci$conf.high)
  rownames(logfc.summary)=pvals$lhs
  return(logfc.summary)
}

########################################################################################################
# extractAllTukeys pulls the estimates for all pairwise differences between means.
##########################################################################################################

extractAllTukeys = function(mdl,df.long){
 require(multcomp)
 require(broom)
 
 # Make an exemplar dataframe with conditions & ref/alt
 exemplar=expand.grid(allele=levels(df.long$allele),condition=levels(df.long$condition))
  contrasts(exemplar$allele)=contrasts(df.long$allele)
  contrasts(exemplar$condition)=contrasts(df.long$condition)
 
 # Make the model matrix to transform the betas
  X=model.matrix(~allele*condition,exemplar)
  rownames(X)=paste(exemplar$condition,exemplar$allele,sep=":")
  
 # Computer the linear terms for all pairwise comparisons between rows of X
  alltuk=combn(rownames(X),2)
  T = X[alltuk[2,],]-X[alltuk[1,],]
  rownames(T)=paste(alltuk[2,],alltuk[1,],sep="-")

 # compute the values using multcomp's glht
  
  logfc=glht(mdl,linfct=T)
  pvals=tidy(summary(logfc))
  ci=tidy(confint(logfc))
  
  logfc.summary = data.frame(logfc=pvals$estimate,
                             se=pvals$std.error,
                             z=pvals$statistic,
                             p=pvals$p.value,
                             lo=ci$conf.low,
                             hi=ci$conf.high)
  rownames(logfc.summary)=pvals$lhs
  return(logfc.summary)
}


########################################################################################################
# extractDescriptiveStats pulls the descriptive stat estimates for each elementXcondition
########################################################################################################

extractDescriptiveStats= function(mdl,df.long){
 require(multcomp)
 require(broom)
 # Make an exemplar dataframe with conditions & ref/alt
 exemplar=expand.grid(allele=levels(df.long$allele),condition=levels(df.long$condition))
  contrasts(exemplar$allele)=contrasts(df.long$allele)
  contrasts(exemplar$condition)=contrasts(df.long$condition)
 
 # Make the model matrix to transform the betas
  X=model.matrix(~allele*condition,exemplar)
  rownames(X)=paste(exemplar$condition,exemplar$allele,sep=":")
  
 mns=glht(mdl,linfct=X)
 pvals=tidy(summary(mns))
 ci=tidy(confint(mns))
  
 descstats=data.frame(mn=pvals$estimate,
                       se=pvals$std.error,
                       lo=ci$conf.low,
                       hi=ci$conf.high)
 rownames(descstats)=pvals$lhs 
 
 return(descstats)
}


########################################################################################################
# collapseBarcodes pulls the barcode averaged sample data
########################################################################################################

collapseBarcodes=function(df.long){
  
   ag = aggregate(logy~sample*condition*allele,df.long,mean)
   rownames(ag)=paste(ag$sample,ag$allele,sep="X")
   ag = ag[,-1]
   
   return(ag)
  
}













