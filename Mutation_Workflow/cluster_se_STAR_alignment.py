import sys
import re
import os
import subprocess

java = "/net/isi-dcnl/ifs/user_data/Seq/jre1.8.0_121/bin/java"
jar_path = "/net/isi-dcnl/ifs/user_data/Seq/"

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

submitAll = "master_STAR_queue.sh"
masterHandle = open(submitAll,"w")
text = "#!/bin/bash\n"
masterHandle.write(text)

jobCount = 0
for file in fileResults:
	resultGZ = re.search("(.*)_S\d+_L\d{3}_R1_001.fastq.gz$",file)
	
	if resultGZ:
		sample = resultGZ.group(1)
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
			text = text + "#$ -l vf="+java_mem+"\n"
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
			text = STAR + " --chimSegmentMin 12 --chimJunctionOverhangMin 12 --genomeDir " + ref+ " --readFilesIn " + read1 + " --runThreadN " +threads+ " --outFileNamePrefix " + starPrefix + " --twopassMode Basic --outSAMstrandField intronMotif\n"
			outHandle.write(text)

			tempDir = outputSubfolder + "/tmp"
			text = "mkdir " + tempDir + "\n"
			outHandle.write(text)	
			
			starSam = outputSubfolder + "/" + sample + "_Aligned.out.sam"																			
			rgBam = outputSubfolder + "/rg.bam"
			text = java + " -Xmx" + java_mem + " -Djava.io.tmpdir="+ tempDir + " -jar "+jar_path+"picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups I=" + starSam + " O=" + rgBam + " SO=coordinate RGID=1 RGLB=RNA-Seq RGPL=Illumina RGPU=COH RGSM=" + sample + "\n"
			outHandle.write(text)

			text = "rm " + starSam + "\n"
			outHandle.write(text)		
			
			duplicateMetrics = outputSubfolder + "/MarkDuplicates_metrics.txt"
			filteredBam = outputSubfolder + "/nodup.bam"
			text = java + " -Xmx" + java_mem + " -Djava.io.tmpdir="+ tempDir + " -jar "+jar_path+"picard-tools-2.5.0/picard.jar MarkDuplicates I=" + rgBam + " O=" + filteredBam + " M=" + duplicateMetrics+" REMOVE_DUPLICATES=true CREATE_INDEX=true\n"
			outHandle.write(text)

			text = "rm " + rgBam + "\n"
			outHandle.write(text)

			trimmedBam = alignmentFolder + "/" + sample + ".bam"	
			text = java + " -Xmx" + java_mem + " -jar "+jar_path+"GenomeAnalysisTK-3.6.jar -T SplitNCigarReads -R " + fa_ref + " -I " + filteredBam + " -o " + trimmedBam+" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS\n"
			outHandle.write(text)

			text = "rm " + filteredBam + "\n"
			outHandle.write(text)
			
			text = "rm " + re.sub(".bam$",".bai",filteredBam) + "\n"
			outHandle.write(text)

			text = "rm -R " + tempDir + "\n"
			outHandle.write(text)
			
			#since code is for .gz files
			text = "rm " + read1 + "\n"
			outHandle.write(text)
