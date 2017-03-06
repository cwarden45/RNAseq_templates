import sys
import re
import os

parameterFile = "parameters.txt"

alignmentFolder = ""
libraryType = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	
	if param == "Alignment_Folder_MAC":
		alignmentFolder = value

	if param == "Read_Pairing":
		libraryType = value

if (libraryType == "") or (libraryType == "[required]"):
	print "Need to enter a value for 'Read_Pairing'!"
	sys.exit()
		
fileResults = os.listdir(alignmentFolder)

statsFile = "RSeQC_stats.txt"
statHandle = open(statsFile, 'w')
text = "Sample\tReverse.Strand.Bias\tmedTIN\n"
statHandle.write(text)

for file in fileResults:
	result = re.search("(.*).bam$",file)
	
	if result:
		sample = result.group(1)
		print sample
		subfolder = alignmentFolder + "/" + sample

		text = sample
		
		strandStat = subfolder + "/" + sample + "_infer_strand.txt"
		inHandle = open(strandStat)
		lines = inHandle.readlines()
					
		for line in lines:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			if libraryType == "SE":
				result2 = re.search("Fraction of reads explained by \"\+-,-\+\": (.*)",line)
				
				if result2:
					formatFrac =  '{0:.2g}'.format(float(result2.group(1)))
					text = text + "\t" +  formatFrac
			elif libraryType == "PE":
				result2 = re.search("Fraction of reads explained by \"1\+-,1-\+,2\+\+,2--\": (.*)",line)
				
				if result2:
					formatFrac =  '{0:.2g}'.format(float(result2.group(1)))
					text = text + "\t" +  formatFrac
			else:
				print "'Read_Pairing' must be 'SE' or 'PE'"
				sys.exit()
			
	
		tinStat = subfolder + "/" + sample + ".summary.txt"
		lineCount = 0
		inHandle = open(tinStat)
		lines = inHandle.readlines()
					
		for line in lines:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			lineCount += 1
			
			if lineCount == 2:
				lineInfo = line.split("\t")
				text = text + "\t" + '{0:.3g}'.format(float(lineInfo[2]))
				
		text = text + "\n"
		statHandle.write(text)
