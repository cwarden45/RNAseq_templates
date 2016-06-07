import sys
import re
import os

genome = "hg19"

bioconductorAnnotations = "TxDb_" + genome + "_exon_annotations.txt"
exonGtf = "TxDb_" + genome + "_gene.gtf"
source = "Bioconductor_" + genome

inHandle = open(bioconductorAnnotations)
lines = inHandle.readlines()

exonOut = open(exonGtf, 'w')

lineCount = 0

for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	lineCount +=1
	
	if lineCount > 1:
		lineInfo = line.split("\t")
		chr = lineInfo[0]
		chr = re.sub("\"","",chr)
		start = lineInfo[1]
		stop = lineInfo[2]
		symbol = lineInfo[4]
		symbol = re.sub("\"","",symbol)
		strand = lineInfo[5]
		strand = re.sub("\"","",strand)
		
		gtfLine = chr + "\t" + source + "\texon\t" + start + "\t" + stop + "\t0.0\t" + strand + "\t.\tgene_id \"" + symbol + "\"\n";
		exonOut.write(gtfLine)