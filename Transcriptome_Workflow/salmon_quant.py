import sys
import re
import os

threads = 2
ref = "/path/to/index_folder"
quantFolder = "salmon_quant"
readsFolder = "Reads"

finishedSamples = ()
fileResults = os.listdir(readsFolder)

for file in fileResults:
	result = re.search("(.*)_\w{6}_L\d{3}_R1_001.fastq$",file)
	
	if result:
		sample = result.group(1)
		if sample not in finishedSamples:
			print sample

			outputSubfolder = quantFolder +"/" + sample
			command = "mkdir " + outputSubfolder
			os.system(command)
									
			read1 = readsFolder + "/" + file
					
			command = "/opt/salmon/bin/salmon quant -l U -i " + ref +  " -p " + str(threads) + " -o " + outputSubfolder + " -r " + read1
			os.system(command)

