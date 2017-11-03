import sys
import re
import os

alignmentFolder = "Oases_Assemblies"
statFile = "Oases_eXpress_stats.txt"

#alignmentFolder = "Oases_Assemblies_Trim"
#statFile = "Oases_eXpress_stats_Trim.txt"

def getPercentAligned(inputfile):
	inHandle = open(inputfile)
	line = inHandle.readline()
		
	totalReads = 0
	alignedReads = 0
		
	lineCount =0
	while line:
		lineCount += 1
			
		if lineCount == 1:
			countResult = re.search("^(\d+)",line)
			totalReads = int(countResult.group(1))

		if lineCount == 5:
			countResult = re.search("^(\d+)",line)
			alignedReads = int(countResult.group(1))
				
		line = inHandle.readline()
			
	inHandle.close()
	
	percentAligned = 100 * float(alignedReads)/float(totalReads)
	return('{0:.2f}'.format(percentAligned) + "%")

outHandle = open(statFile, "w")
text = "Sample\tLongest.Contig.KB\tMean.Contig.Length.KB\tPercent.Aligned\tNumber.of.Contigs\tTPM10.Contigs\n"
outHandle.write(text)

fileResults = os.listdir(alignmentFolder)

for folder in fileResults:
	sample = folder
	subfolder = os.path.join(alignmentFolder, folder)
	
	if os.path.isdir(subfolder):
		express_table = os.path.join(subfolder, "results.xprs")
		if os.path.isfile(express_table):
			print sample
			
			flagstat_table = os.path.join(subfolder, sample+"_flagstat.txt")
			percentAligned = getPercentAligned(flagstat_table)
			
			longestContig=0
			contigNum = 0
			contigSum = 0
			tpm10count = 0
			
			inHandle = open(express_table)
			line = inHandle.readline()
				
			totalReads = 0
			alignedReads = 0
				
			lineCount =0
			while line:
				lineCount += 1
					
				if lineCount > 1:
					line = re.sub("\r","",line)
					line = re.sub("\n","",line)
					
					lineInfo = line.split("\t")
					length = float(lineInfo[2]) / 1000
					solvable = lineInfo[13]
					tpm = float(lineInfo[14])
					if solvable == "T":
						contigNum += 1
						contigSum += length
						if tpm > 10:
							tpm10count += 1
							
						if length > longestContig:
							longestContig = length
						
				line = inHandle.readline()
					
			inHandle.close()
			
			meanLength = float(contigSum) / float(contigNum)
			text = sample + "\t" + str(longestContig) + "\t" + '{0:.2f}'.format(meanLength) + "kb\t" + percentAligned + "\t" + str(contigNum) + "\t" + str(tpm10count) + "\n"
			outHandle.write(text)