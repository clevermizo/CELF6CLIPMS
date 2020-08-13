########## c6utrbccount.py ############
#### used for counting barcodes in PTRE-Seq samples ######
from Bio import SeqIO as seq
from Bio.Seq import Seq
import sys
import os
import subprocess
import numpy as np
import pandas as pd
import pdb
import re
import csv

r1=sys.argv[1]
outputFile=sys.argv[2]

print(r1)
print(outputFile)


bcseqlookup="bc_and_arraySeqs.tab"
bcnamelookup="bc_and_arrayNames.tab"
bcseq=pd.read_table(bcseqlookup,sep="\t",header=None)
bcseq.columns=["bc","sequence"]
bcnames=pd.read_table(bcnamelookup,sep="\t",header=None)
bcseq = bcseq.assign(seqName=bcnames[[1]])
c6pat=re.compile("c6targ")
c6match=list()

for i in range(len(bcseq.index)):
    got_a_match=c6pat.search(bcseq["seqName"][i])
    if got_a_match:
      c6match.append(True)
    else:
      c6match.append(False)

bcseq = bcseq[c6match]
Nbc = len(bcseq.index)

number_of_reads = int(os.popen("zcat " + r1 + ".fastq.gz | wc -l").read().split()[0])/4
os.system("gunzip -k "+r1+".fastq.gz")

adapter1="TCATGTACC"
adapter2="CAGGTGTACC"
adapter3="GTTCCTGTACC"
adapter4="AGCAGCTGTACC"

bcseq = bcseq.assign(bcCount=0)
for i in bcseq.index:
  bctmp = str.upper(str(Seq(bcseq.bc[i]).reverse_complement()))
  searchPattern = "^" + adapter1 + bctmp + "CTCGAG|^" + adapter2 + bctmp + "CTCGAG|^" + adapter3 + bctmp + "CTCGAG|^" + adapter4 + bctmp + "CTCGAG"
  lineCount=os.popen("grep -E '" + searchPattern + "' " + r1 + ".fastq | wc -l").read()
  lineCount=int(lineCount)
  bcseq=bcseq.set_value(i,'bcCount',lineCount)
  
# Write output data frame
bcseq.to_csv(outputFile,sep="\t",quoting=csv.QUOTE_NONE)
