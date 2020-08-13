options(stringsAsFactors = FALSE)
prettyPValue=Vectorize(function(pvalue){
  if(pvalue>=0.0001){
    pp=as.character(round(pvalue,digits=4))
    return(pp)
  }else{
    expnt=floor(log10(pvalue))
    coef=(10^log10(pvalue))/(10^expnt)
    coef=round(coef,digits = 2)
    pp=paste(coef,"E",expnt,sep="")
    return(pp)
  }
})
makePrettyCi=function(x,Ndec=4){
  if(is.vector(x) & length(x)==3){
    xrnd=round(x,digits = Ndec)
    xtext=as.character(xrnd[1])
    xtext2=str_c(xrnd[2:3],collapse =":")
    xtext3=paste(xtext," [",xtext2,"]",sep="")
    return(xtext3)
  }else{
    stop("Input must be triplet: c(estimate,lower,upper)")
  }
}
unwrapPrettyCi=function(xtext){
  if(length(grep("\\[",xtext))==1 & length(grep("\\]",xtext))==1){
    xrnd=str_split(xtext," ",simplify = TRUE)
    xrnd1=xrnd[1]
    xrnd2=as.vector(str_split(str_replace_all(xrnd[2],"[\\[\\]]",""),":",simplify=TRUE))
    return(as.numeric(c(xrnd1,xrnd2)))
  }else{
    stop("Input doesn't look like a pretty CI")
  }
}
computeConfIntOverlap=function(c1,c2){
  #rank order the confidence intervals
  #first is "a" second is "b"
  if(sum(is.na(c(c1,c2)))>0){
    overlap=NA
    return(overlap)
  }else{
    if(max(c2)>=max(c1)){
      a=c1
      b=c2
    }else if(max(c1)>=max(c2)){
      a=c2
      b=c1
    }
    # overlap them
    d1=max(a)-min(b)
    # if d1<=0 then there is no overlap
    d2=min(b)-min(a)
    # if d2<=0 then there is 100% overlap
    if(d1<=0&d2>0){
      overlap=0
      return(overlap)
    }else if (d1>0 & d2<=0){
      overlap=100
      return(overlap)
    }else if (d1>0 & d2>0){
      overlap=(d1/diff(a))*100
      return(overlap)
    }
  }
}
makePrettyPercent=Vectorize(function(pct,digits=2){
  if(pct==0){
    pct.out="0%"
  }else if(round(pct,digits)==0){
    pct.out=paste("<",1/(10^digits),"%",sep="")
  }else{
    pct.out=paste(round(pct,digits),"%",sep="")
  }
})
corSummaryTable=function(CorMatrix){
  Cnames=matrix(rep(colnames(CorMatrix),
                    length(colnames(CorMatrix))),
                nrow=nrow(CorMatrix),ncol=ncol(CorMatrix))
  Creport=data.frame(sample1=Cnames[upper.tri(Cnames)],
                     sample2=t(Cnames)[upper.tri(t(Cnames))],
                     R=CorMatrix[upper.tri(CorMatrix)])
  return(Creport)
}