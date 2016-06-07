import sys
import re
import os

root = "/path/to/output"
alignmentFolder = "/path/to/folder/"
gtf_file = "/path/to/file"
lncRNAgtf = "/path/to/file"

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
		
			countsFile = root + sample + "_gene_counts.txt"
			countsFile = os.path.join(root,countsFile)
			command = "htseq-count -f bam -s no " + nameSortedBam + " " + gtf_file + " > " + countsFile
			os.system(command)
			
			countsFile = root + sample + "_lncRNA_counts.txt"
			countsFile = os.path.join(root,countsFile)
			command = "htseq-count -f bam -s no -i gene_name " + nameSortedBam + " " + lncRNAgtf + " > " + countsFile
			os.system(command)
	
			command = "rm " + nameSortedBam
			os.system(command)

		
