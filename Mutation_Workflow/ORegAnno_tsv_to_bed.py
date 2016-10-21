import sys
import re
import os

build = "hg19"
inTSV = "ORegAnno_Combined_2015.12.22.tsv"
outBED = "hg19_ORegAnno.bed"

#build = "hg38"
#inTSV = "ORegAnno_Combined_2015.12.22.tsv"
#outBED = "hg38_ORegAnno.bed"

outHandle = open(outBED, "w")

inHandle = open(inTSV)
line = inHandle.readline()

while line:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	annID = lineInfo[0]
	species = lineInfo[1]
	type = lineInfo[3]
	gene = lineInfo[4]
	source=lineInfo[5]
	pubmed = lineInfo[11]
	annBuild = lineInfo[13]
	strand = lineInfo[14]
	if (strand == "1") or (strand == " 1"):
		strand = "+"
	elif strand == "-1":
		strand = "-"
	elif (strand != "Strand") and (strand != "NULL")and (strand != "N/A"):
		print "Need to provide strand mapping for |" + strand + "|"
		sys.exit()
	chr = lineInfo[15]
	start = lineInfo[16]
	stop = lineInfo[17]
	
	if build == annBuild:
		chr = re.sub("\s+","",chr)
		start = re.sub("\s+","",start)
		stop = re.sub("\s+","",stop)
		annID = re.sub("\s+","",annID)
		text = chr + "\t" + start + "\t" + stop + "\t" + annID + "\t0\t" + strand + "\n"
		outHandle.write(text)
	line = inHandle.readline()