import sys
import re
import os

alignmentFolder = "R:\\path\\to\\folder\\"
countsFile = "total_read_counts.txt"

outHandle = open(countsFile, 'w')
text = "Sample\tTotal.Reads\n"
outHandle.write(text)

fileResults = os.listdir(alignmentFolder)

for subfolder in fileResults:
	result = re.search("^\d",subfolder)

	if result:		
		readInfoFile = alignmentFolder + subfolder + "\\prep_reads.info"
		
		if os.path.isfile(readInfoFile):
			sample = subfolder
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