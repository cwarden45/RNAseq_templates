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
		
	if param == "aligned_stats_file":
		countsFile = value
		
outHandle = open(countsFile, 'w')
text = "Sample\taligned.reads\n"
outHandle.write(text)

fileResults = os.listdir(alignmentFolder)

for subfolder in fileResults:
	result = re.search("^\d",subfolder)

	if result:		
		readInfoFile = alignmentFolder + "/" + subfolder + "/align_summary.txt"
		
		if os.path.isfile(readInfoFile):
			sample = subfolder
			print sample
			
			inHandle = open(readInfoFile)
			lines = inHandle.readlines()
			
			totalReads = ""
			
			for line in lines:
				line = re.sub("\n","",line)
				line = re.sub("\r","",line)

				#print line
				result2 = re.search("\s+Mapped\s+:\s+(\d+)",line)
				if result2:
					alignedReads = result2.group(1)
					print alignedReads

					text = sample + "\t" + alignedReads + "\n";
					outHandle.write(text)