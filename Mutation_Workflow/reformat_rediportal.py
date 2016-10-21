import sys
import re
import os

inputFile = "REDIportal_table1_full.txt"
vcfFile = "hg19_REDIportal.vcf"
bedFile = "hg19_REDIportal.bed"

inHandle = open(inputFile)
line = inHandle.readline()
			
outHandle = open(vcfFile, "w")
bedHandle = open(bedFile, "w")

while line:
	line = re.sub("\r","",line)
	line = re.sub("\n","",line)
	
	lineInfo = line.split("\t")
	
	chr = lineInfo[0]
	pos = lineInfo[1]
	ref = lineInfo[2]
	var = lineInfo[3]
	
	text = chr + "\t" + pos + "\t.\t" + ref + "\t" + var + "\t0\tPASS\t.\t.\n"
	outHandle.write(text)
	
	text = chr + "\t" + pos + "\t" + pos + "\n"
	bedHandle.write(text)
	
	line = inHandle.readline()