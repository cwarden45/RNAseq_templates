import sys
import re
import os

parameterFile = "parameters.txt"

threads = ""
ref = ""
quantFolder = ""
readsFolder = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	
	if param == "Threads":
		threads = value
		
	if param == "Salmon_index":
		ref = value

	if param == "Raw_Code_MAC":
		quantFolder = value
		
	if param == "Reads_Folder_MAC":
		readsFolder = value	
		
finishedSamples = ()
fileResults = os.listdir(readsFolder)

for file in fileResults:
	result = re.search("(.*)_\w{6}_L\d{3}_R1_001.fastq$",file)
	result2 = re.search("(.*)_\w{6}_L\d{3}_R2_001.fastq$",file)
	
	if result:
		sample = result.group(1)
		if sample not in finishedSamples:
			print sample

			outputSubfolder = quantFolder +"/" + sample
			command = "mkdir " + outputSubfolder
			os.system(command)
									
			read1 = readsFolder + "/" + file
			read2 = re.sub("_R1_001.fastq$","_R2_001.fastq",read1)
			
			if os.path.isfile(read2):
				print "Reverse Read Found --> Running PE Quantification"
				command = "/opt/salmon/bin/salmon quant -l U -i " + ref +  " -p " + threads + " -o " + outputSubfolder + " -1 " + read1 + " -2 " + read2
				os.system(command)
			else:
				print "Running SE Quantification"
				command = "/opt/salmon/bin/salmon quant -l U -i " + ref +  " -p " + threads + " -o " + outputSubfolder + " -r " + read1
				os.system(command)

