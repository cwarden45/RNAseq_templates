import sys
import re
import os

parameterFile = "parameters.txt"
finishedSamples = []

threads = ""
ref = ""
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

	if param == "Reads_Folder_MAC":
		readsFolder = value
	
	if param == "Alignment_Folder_MAC":
		alignmentFolder = value
		
	if param == "Bowtie2_Ref":
		ref = value
		
	if param == "Threads":
		threads = value
		
if (threads == "") or (threads == "[required]"):
	print "Need to enter a value for 'Threads'!"
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

for file in fileResults:
	result = re.search("(.*)_\w{6}_L\d{3}_R1_001.fastq$",file)
	
	if result:
		sample = result.group(1)
		if sample not in finishedSamples:
			print sample

			outputSubfolder = alignmentFolder +"/" + sample
			command = "mkdir " + outputSubfolder
			os.system(command)
									
			read1 = readsFolder + "/" + file
								
			print "\n\nAlign via TopHat\n\n"
			command = "tophat -o " + outputSubfolder + " -p " + threads + " --no-coverage-search " + ref + " " + read1 
			os.system(command)
								
			topHatBam = outputSubfolder + "/accepted_hits.bam";																			
			userBam = alignmentFolder + "/" + sample + ".bam";
								
			print "\n\nCreate Sorted BAM File\n\n"
			command = "/opt/samtools-1.3/samtools sort " + topHatBam + " -o " + userBam
			os.system(command)

			command = "rm " + topHatBam
			os.system(command)

			print "\n\nIndexing BAM File\n\n"
			command = "samtools index " + userBam
			os.system(command)
