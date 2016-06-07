import sys
import re
import os

parameterFile = "parameters.txt"

alignmentFolder = ""
gtf_file = ""
lncRNAgtf = ""

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
		
	if param == "txGTF_MAC":
		gtf_file = value

	if param == "lncRNA_GTF_MAC":
		lncRNAgtf = value		

fileResults = os.listdir(alignmentFolder)

finishedSamples = ()

for file in fileResults:
	result = re.search(".bam$",file)
	fullPath = os.path.join(alignmentFolder, file)
	
	if result:
		sample = re.sub(".bam$","",file)
		if sample not in finishedSamples:
			print sample
		
			nameSortedBam = sample + ".name.sort.bam"
			nameSortedBam = os.path.join(alignmentFolder, nameSortedBam)
			command = "/opt/samtools-1.3/samtools sort -n -o " + nameSortedBam + " " + fullPath
			os.system(command)
		
			countsFile = sample + "_gene_counts.txt"
			command = "htseq-count -f bam -s no " + nameSortedBam + " " + gtf_file + " > " + countsFile
			os.system(command)
			
			countsFile = sample + "_lncRNA_counts.txt"
			command = "htseq-count -f bam -s no -i gene_name " + nameSortedBam + " " + lncRNAgtf + " > " + countsFile
			os.system(command)
	
			command = "rm " + nameSortedBam
			os.system(command)

		
