import sys
import re
import os

#don't really use unaligned reads, but can check for common causes of low alignment rate

readsFolder = "../../Reads"
fastq_screen = "/path/to/fastq_screen_v0.11.3/fastq_screen"
config_file = "/path/to/fastq_screen_v0.11.3/fastq_screen.conf"
finishedSamples = ("")

fileResults = os.listdir(readsFolder)

for file in fileResults:
	result = re.search("(.*)_S\d+_L\d{3}_R1_001.fastq.gz$",file)
	
	if result:
		sample = result.group(1)
		
		if (sample not in finishedSamples):
			print sample
			fastq = readsFolder + "/" + file
			
			command = fastq_screen + " --subset 1000 --threads 4 " + fastq + " --conf " + config_file
			os.system(command)