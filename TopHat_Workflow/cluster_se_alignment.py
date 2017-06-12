import sys
import re
import os
import subprocess

parameterFile = "parameters.txt"
finishedSamples = ()

threads = ""
ref = ""
alignmentFolder = ""
readsFolder = ""
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

	if param == "Reads_Folder_MAC":
		readsFolder = value
	
	if param == "Alignment_Folder_MAC":
		alignmentFolder = value
		
	if param == "Bowtie2_Ref":
		ref = value
		
	if param == "Threads":
		threads = value

	if param == "strand":
		strandType = value
		
	if param == "Cluster_Email":
		email = value
		
if (threads == "") or (threads == "[required]"):
	print "Need to enter a value for 'Threads'!"
	sys.exit()

if (strandType == "") or (strandType == "[required]"):
	print "Need to enter a value for 'strand'!"
	sys.exit()

if (email == "") or (email == "[required]"):
	print "Need to enter a value for 'Cluster_Email'!"
	sys.exit()
	
if (ref == "") or (ref == "[required]"):
	print "Need to enter a value for 'Bowtie2_Ref'!"
	sys.exit()

if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder_MAC'!"
	sys.exit()
	
if (alignmentFolder == "") or (alignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder_MAC'!"
	sys.exit()
	
fileResults = os.listdir(readsFolder)

submitAll = "master_tophat_queue.sh"
masterHandle = open(submitAll,"w")
text = "#!/bin/bash\n"
masterHandle.write(text)

jobCount = 0
for file in fileResults:
	result = re.search("(.*)_S\d+_L\d{3}_R1_001.fastq$",file)
	
	if result:
		sample = result.group(1)
		jobCount += 1
		
		if (sample not in finishedSamples):
			print sample
			shellScript = sample + ".sh"
			text = "qsub " + shellScript + "\n"
			masterHandle.write(text)

			outHandle = open(shellScript, "w")
			text = "#!/bin/bash\n"
			text = text + "#$ -M "+email+"\n"
			text = text + "#$ -m bea\n"
			text = text + "#$ -N CW"+str(jobCount)+"\n"
			text = text + "#$ -q all.q\n"
			text = text + "#$ -pe shared "+threads+"\n"
			text = text + "#$ -l vf=4G\n"
			text = text + "#$ -j yes\n"
			text = text + "#$ -o CW"+str(jobCount)+".log\n"
			text = text + "#$ -cwd\n"
			text = text + "#$ -V\n"
			outHandle.write(text)
				
			outputSubfolder = alignmentFolder +"/" + sample
			text = "mkdir " + outputSubfolder + "\n"
			outHandle.write(text)
										
			read1 = readsFolder + "/" + file

			tophatStrand = ""
			if strandType == "no":
				tophatStrand = "fr-unstranded"
			elif strandType == "reverse":
				tophatStrand = "fr-firststrand"
			elif strandType == "yes":
				tophatStrand = "fr-secondstrand"
			else:
				print "Need to provide TopHat mapping for strand: " + strandType
				sys.exit()
			
			text = "tophat -o " + outputSubfolder + " -p " + threads + " --no-coverage-search --library-type " + tophatStrand + " " + ref + " " + read1 + "\n" 
			outHandle.write(text)
									
			topHatBam = outputSubfolder + "/accepted_hits.bam"																			
			userBam = alignmentFolder + "/" + sample + ".bam"
			
			sortPrefix = re.sub(".bam$","",userBam)
			text = "samtools sort " + topHatBam + " " + sortPrefix + "\n"
			outHandle.write(text)
			
			text = "rm " + topHatBam + "\n"
			outHandle.write(text)

			text = "samtools index " + userBam + "\n"
			outHandle.write(text)

			text = "gzip " + read1
			outHandle.write(text)
				
			#run shell script to submit all samples.  Otherwise, you can run qsub manually.
			#command = "qsub " + shellScript
			#os.system(command)
			#subprocess.run(command, shell=True)
