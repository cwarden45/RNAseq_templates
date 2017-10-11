import sys
import re
import os

bed = ""
annovar = ""

for arg in sys.argv:
	bedResult = re.search("^--bed=(.*)",arg)
	annovarResult = re.search("^--annovar=(.*)",arg)
	helpResult = re.search("^--help",arg)
	
	if bedResult:
		bed = bedResult.group(1)

	if annovarResult:
		annovar = annovarResult.group(1)
		
	if helpResult:
		print "Usage: python bed_to_ANNOVAR.py --bed=file.bed --annovar=file.avinput\n"
		print "--bed : .bed with with bedtools annotations\n"
		print "--annovar : annotations in ANNOVAR format, for summary\n"
		sys.exit()

outHandle = open(annovar,"w")
			
inHandle = open(bed)
line = inHandle.readline()
			
while line:
	lineInfo = line.split("\t")
	chr = lineInfo[0]
	start = lineInfo[1]
	stop = lineInfo[2]
	ref = lineInfo[3]
	var = lineInfo[4]
	type = lineInfo[5]
	somatic_pvalue = lineInfo[6]
	totalCov = lineInfo[7]
	oreg = lineInfo[11]
					
	text = chr + "\t" + start+ "\t" + stop + "\t" + ref + "\t" + var + "\t" + type + "\t" + somatic_pvalue+ "\t" + oreg + "\n"
	outHandle.write(text)
						
	line = inHandle.readline()
inHandle.close()
outHandle.close()