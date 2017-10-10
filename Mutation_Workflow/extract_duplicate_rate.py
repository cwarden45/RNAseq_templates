import sys
import re
import os
import subprocess

parameterFile = "parameters.txt"
statFile = "duplicate_statistics.txt"

readsFolder = ""
AlignmentFolder = ""

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

	if param == "Alignment_Folder":
		AlignmentFolder = value
		
if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()
	
if (AlignmentFolder == "") or (AlignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder'!"
	sys.exit()


statHandle = open(statFile,"w")
text = "SampleID\tSeqID\tuserID\tTotalReads\tPercent.Duplicate\n"
statHandle.write(text)
	
fastqcFolder = readsFolder + "/QC"
fileResults = os.listdir(readsFolder)

for file in fileResults:
	result = re.search("(.*)_S\d+_L\d{3}_R1_001.fastq.gz$",file)
	
	if result:
		sample = result.group(1)
		print sample
		
		r2 = re.search("^(\d+)_.*",sample)
		seqID = r2.group(1)
		
		shortID = re.sub(seqID + "_","",sample)
		
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
		
		#get duplicate rate from Picard MarkDuplicates
		sampleSubfolder = AlignmentFolder + "/" + sample
		duplicateMetrics = sampleSubfolder + "/MarkDuplicates_metrics.txt"

		inHandle = open(duplicateMetrics)
		line = inHandle.readline()
		
		lineCount = 0
		
		while line:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			lineCount += 1
			
			if lineCount == 8:
				lineInfo = line.split("\t")
				
				percentDuplicate = 100*float(lineInfo[8])
				percentDuplicate =  '{0:.2f}'.format(percentDuplicate) + "%"
				
			line = inHandle.readline()
		
		inHandle.close()
		
		text = sample + "\t" + seqID + "\t" + shortID + "\t" + totalReads + "\t" + percentDuplicate + "\n"
		statHandle.write(text)