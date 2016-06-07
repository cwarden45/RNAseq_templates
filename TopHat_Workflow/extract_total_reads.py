import sys
import re
import os

parameterFile = "parameters.txt"

alignmentFolder = ""
countsFile = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	
	if param == "Alignment_Folder_PC":
		alignmentFolder = value
		
	if param == "total_counts_file":
		countsFile = value


outHandle = open(countsFile, 'w')
text = "Sample\tTotal.Reads\n"
outHandle.write(text)

fileResults = os.listdir(alignmentFolder)

for subfolder in fileResults:
	result = re.search("^\d",subfolder)

	if result:		
		readInfoFile = alignmentFolder + subfolder + "/prep_reads.info"
		
		if os.path.isfile(readInfoFile):
			sample = re.sub("_\w{6}_L999$","",subfolder)
			print sample
			
			inHandle = open(readInfoFile)
			lines = inHandle.readlines()
			
			totalReads = ""
			
			for line in lines:
				line = re.sub("\n","",line)
				line = re.sub("\r","",line)

				result2 = re.search("^reads_in =(\d+)",line)
				if result2:
					totalReads = result2.group(1)			

					text = sample + "\t" + totalReads + "\n";
					outHandle.write(text)