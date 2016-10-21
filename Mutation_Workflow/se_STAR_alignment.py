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

jobCount = 0
for file in fileResults:
	resultGZ = re.search("(.*)_\w{6}_L\d{3}_R1_001.fastq.gz$",file)
	
	if resultGZ:
		sample = resultGZ.group(1)
		jobCount += 1
		
		if (sample not in finishedSamples):
			print sample
			
			read1 = re.sub(".gz$","",readsFolder + "/" + file)
			command = "gunzip -c " + readsFolder + "/" + file + " > " + read1
			os.system(command)
			
			outputSubfolder = alignmentFolder +"/" + sample
			command = "mkdir " + outputSubfolder + "\n"
			os.system(command)
			
			starPrefix = outputSubfolder + "/" + sample + "_"

			#got segmentation fault when I tried to let STAR create sorted .bam file
			command = "/opt/STAR/bin/Linux_x86_64/STAR --genomeDir " + ref+ " --readFilesIn " + read1 + " --runThreadN " +threads+ " --outFileNamePrefix " + starPrefix + " --twopassMode Basic --outSAMstrandField intronMotif"
			os.system(command)

			#since this code is for .gz files
			command = "rm " + read1
			os.system(command)
									
			starSam = outputSubfolder + "/" + sample + "_Aligned.out.sam"																			
			rgBam = outputSubfolder + "/rg.bam"
			command = "java -Xmx" + java_mem + " -jar /opt/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups I=" + starSam + " O=" + rgBam + " SO=coordinate RGID=1 RGLB=RNA-Seq RGPL=Illumina RGPU=COH RGSM=" + sample
			os.system(command)

			command = "rm " + starSam
			os.system(command)
			
			duplicateMetrics = outputSubfolder + "/MarkDuplicates_metrics.txt"
			filteredBam = outputSubfolder + "/nodup.bam"
			command = "java -Xmx" + java_mem + " -jar /opt/picard-tools-2.5.0/picard.jar MarkDuplicates I=" + rgBam + " O=" + filteredBam + " M=" + duplicateMetrics+" REMOVE_DUPLICATES=true CREATE_INDEX=true"
			os.system(command)

			command = "rm " + rgBam
			os.system(command)

			trimmedBam = alignmentFolder + "/" + sample + ".bam"	
			command = "java -Xmx" + java_mem + " -jar /opt/GenomeAnalysisTK-3.6.jar -T SplitNCigarReads -R " + fa_ref + " -I " + filteredBam + " -o " + trimmedBam+" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS"
			os.system(command)

			command = "rm " + filteredBam
			os.system(command)
			
			command = "rm " + re.sub(".bam$",".bai",filteredBam)
			os.system(command)