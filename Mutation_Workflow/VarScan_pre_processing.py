import sys
import re
import os
import subprocess

parameterFile = "parameters.txt"
finishedSamples = ()

java_mem = ""
threads = ""
ref = ""
fa_ref = ""
alignmentFolder = ""
readsFolder = ""

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
		alignmentFolder = value
		
	if param == "STAR_Ref":
		ref = value
		
	if param == "FA_Ref":
		fa_ref = value
		
	if param == "Java_Mem":
		java_mem = value
		
	if param == "Threads":
		threads = value

if (java_mem == "") or (java_mem == "[required]"):
	print "Need to enter a value for 'Java_Mem'!"
	sys.exit()
		
if (threads == "") or (threads == "[required]"):
	print "Need to enter a value for 'Threads'!"
	sys.exit()

if (fa_ref == "") or (fa_ref == "[required]"):
	print "Need to enter a value for 'FA_Ref'!"
	sys.exit()
	
if (ref == "") or (ref == "[required]"):
	print "Need to enter a value for 'STAR_Ref'!"
	sys.exit()

if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()
	
if (alignmentFolder == "") or (alignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder'!"
	sys.exit()
	
fileResults = os.listdir(readsFolder)

for file in fileResults:
	resultGZ = re.search("(.*)_\w{6}_L\d{3}_R1_001.fastq.gz$",file)
	
	if resultGZ:
		sample = resultGZ.group(1)
		
		if (sample not in finishedSamples):
			print sample
			
			outputSubfolder = alignmentFolder +"/" + sample

			trimmedBam = alignmentFolder + "/" + sample + ".bam"	
			pileup = outputSubfolder + "/" + sample + ".pileup"
			
			command = "/opt/samtools-1.3.1/samtools mpileup -C50 -Q 20 -q 10 -f " + fa_ref + "  " + trimmedBam + " > " + pileup
			os.system(command)
