import sys
import re
import os

parameterFile = "parameters.txt"

readsFolder = ""
countsFile = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	
	if param == "Reads_Folder":
		readsFolder = value
		
	if param == "total_counts_file":
		countsFile = value


outHandle = open(countsFile, 'w')
text = "Sample\tTotal.Reads\n"
outHandle.write(text)

fastqcFolder = readsFolder + "/QC"
fileResults = os.listdir(readsFolder)

for file in fileResults:
	result = re.search("(.*)_S\d+_L\d{3}_R1_001.fastq.gz$",file)
	
	if result:
		sample = result.group(1)
		print sample
			
		#get total reads from FastQC
		fastqcPrefix = re.sub(".fastq.gz","",file)
		fastQCtext = fastqcFolder + "/" + fastqcPrefix + "_fastqc/fastqc_data.txt"
			
		inHandle = open(fastQCtext)
		line = inHandle.readline()
			
		lineCount = 0
			
		while line:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			lineCount += 1
				
			if lineCount == 7:
					
				totalResult = re.search("Total Sequences\t(\d+)",line)
				if totalResult:
					totalReads = totalResult.group(1)
				else:
					print "Problem parsing FastQC file!\n"
					sys.exit()
				
			line = inHandle.readline()
			
		inHandle.close()		

		text = sample + "\t" + totalReads + "\n";
		outHandle.write(text)