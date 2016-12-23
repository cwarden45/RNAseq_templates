import sys
import re
import os

parameterFile = "parameters.txt"
finishedSamples = ()

alignmentFolder = ""
gtf_file = ""
lncRNAgtf = ""
strandType = ""
email = ""

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
		
	if param == "Cluster_Email":
		email = value

if (strandType != "yes") and (strandType != "no") and (strandType != "reverse"):
	print "Need to provide HT-Seq mapping for strand: " + strandType
	sys.exit()
		
fileResults = os.listdir(alignmentFolder)

submitAll = "master_htseq_queue.sh"
masterHandle = open(submitAll,"w")
text = "#!/bin/bash\n"
masterHandle.write(text)

jobCount = 0

for file in fileResults:
	result = re.search(".bam$",file)
	fullPath = os.path.join(alignmentFolder, file)
	
	if result:
		sample = re.sub(".bam$","",file)
		sortResult = re.search(".name.sort.bam",file)
		
		if not sortResult:
			jobCount += 1
		
		if (sample not in finishedSamples) and (not sortResult):
			print sample
	
			shellScript = "htseq_" + sample + ".sh"
			text = "qsub " + shellScript + "\n"
			masterHandle.write(text)
			
			outHandle = open(shellScript, "w")
			text = "#!/bin/bash\n"
			text = text + "#$ -M "+email+"\n"
			text = text + "#$ -m bea\n"
			text = text + "#$ -N cwHT"+str(jobCount)+"\n"
			text = text + "#$ -q all.q\n"
			text = text + "#$ -j yes\n"
			text = text + "#$ -o cwHT"+str(jobCount)+".log\n"
			text = text + "#$ -cwd\n"
			text = text + "#$ -V\n"
			outHandle.write(text)
	
			nameSortedBam = sample + ".name.sort.bam"
			sortPrefix = re.sub(".bam$","",nameSortedBam)
			text = "samtools sort -n " + fullPath + " " + sortPrefix + "\n"
			outHandle.write(text)
		
			countsFile = sample + "_gene_counts.txt"
			text = "htseq-count -f bam -s " + strandType + " " + nameSortedBam + " " + gtf_file + " > " + countsFile + "\n"
			outHandle.write(text)
			
			countsFile = sample + "_lncRNA_counts.txt"
			text = "htseq-count -f bam -s " + strandType + " -i gene_name " + nameSortedBam + " " + lncRNAgtf + " > " + countsFile + "\n"
			#os.system(command)
	
			text = "rm " + nameSortedBam + "\n"
			outHandle.write(text)

		
