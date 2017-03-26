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
email = ""
STAR = ""

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

	if param == "Cluster_Email":
		email = value

	if param == "Path_to_STAR":
		STAR = value

if (STAR == "") or (STAR == "[required]"):
	print "Need to enter a value for 'Path_to_STAR'!"
	sys.exit()
		
if (email == "") or (email == "[required]"):
	print "Need to enter a value for 'Cluster_Email'!"
	sys.exit()
		
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

jobCount = 0
for file in fileResults:
	resultGZ = re.search("(.*)_\w{6}_L\d{3}_R1_001.fastq.gz$",file)
	
	if resultGZ:
		sample = resultGZ.group(1)
		jobCount += 1
		
		if (sample not in finishedSamples):
			print sample
			
			shellScript = sample + ".sh"
			outHandle = open(shellScript, "w")
			text = "#!/bin/bash\n"
			text = text + "#$ -M "+email+"\n"
			text = text + "#$ -m bea\n"
			text = text + "#$ -N CW"+str(jobCount)+"\n"
			text = text + "#$ -q all.q\n"
			text = text + "#$ -pe shared "+threads+"\n"
			text = text + "#$ -l vf=16G\n"
			text = text + "#$ -j yes\n"
			text = text + "#$ -o CW"+str(jobCount)+".log\n"
			text = text + "#$ -cwd\n"
			text = text + "#$ -V\n"
			outHandle.write(text)
				
			outputSubfolder = alignmentFolder +"/" + sample
			text = "mkdir " + outputSubfolder + "\n"
			outHandle.write(text)

			read1 = re.sub(".gz$","",readsFolder + "/" + file)
			text = "gunzip -c " + readsFolder + "/" + file + " > " + read1+ "\n"
			outHandle.write(text)		
			
			starPrefix = outputSubfolder + "/" + sample + "_"
			text = STAR + " --genomeDir " + ref+ " --readFilesIn " + read1 + " --runThreadN " +threads+ " --outFileNamePrefix " + starPrefix + " --twopassMode Basic --outSAMstrandField intronMotif\n"
			outHandle.write(text)
			
			#since code is for .gz files
			text = "rm " + read1 + "\n"
			outHandle.write(text)
