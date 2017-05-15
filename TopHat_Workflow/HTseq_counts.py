import sys
import re
import os

parameterFile = "parameters.txt"

alignmentFolder = ""
gtf_file = ""
lncRNAgtf = ""
strandType = ""

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

	if param == "strand":
		strandType = value

if (strandType != "yes") and (strandType != "no") and (strandType != "reverse"):
	print "Need to provide HT-Seq mapping for strand: " + strandType
	sys.exit()
		
fileResults = os.listdir(alignmentFolder)

finishedSamples = ()

for file in fileResults:
	result = re.search(".bam$",file)
	fullPath = os.path.join(alignmentFolder, file)
	
	if result:
		sample = re.sub(".bam$","",file)
		sortResult = re.search(".name.sort.bam",file)
		if (sample not in finishedSamples) and (not sortResult):
			print sample

			countsFile = sample + "_gene_counts.txt"
			command = "htseq-count -f bam -r pos -s " + strandType + " " + fullPath + " " + gtf_file + " > " + countsFile
			os.system(command)
			
			countsFile = sample + "_lncRNA_counts.txt"
			command = "htseq-count -f bam -r pos -s " + strandType + " -i gene_name " + fullPath + " " + lncRNAgtf + " > " + countsFile
			#os.system(command)
			
			## switch to code below to use name-sorted BAM ##
		
			#nameSortedBam = sample + ".name.sort.bam"
			#nameSortedBam = os.path.join(alignmentFolder, nameSortedBam)
			#command = "/opt/samtools-1.3/samtools sort -n -o " + nameSortedBam + " " + fullPath
			#os.system(command)
		
			#countsFile = sample + "_gene_counts.txt"
			#command = "htseq-count -f bam -s " + strandType + " " + nameSortedBam + " " + gtf_file + " > " + countsFile
			#os.system(command)
			
			#countsFile = sample + "_lncRNA_counts.txt"
			#command = "htseq-count -f bam -s " + strandType + " -i gene_name " + nameSortedBam + " " + lncRNAgtf + " > " + countsFile
			#os.system(command)
	
			#command = "rm " + nameSortedBam
			#os.system(command)

		
