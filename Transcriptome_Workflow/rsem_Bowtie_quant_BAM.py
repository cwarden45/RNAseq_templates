import sys
import re
import os

threads = 4
quantFolder = "RSEM_Bowtie_FQ"
alignmentFolder = "Bowtie_alignment"

finishedSamples = ()
fileResults = os.listdir(alignmentFolder)

for file in fileResults:
	result = re.search("(.*).bam$",file)
	
	if result:
		sample = result.group(1)
		if sample not in finishedSamples:
			print sample

			outputSubfolder = quantFolder +"/" + sample
			command = "mkdir " + outputSubfolder
			os.system(command)
									
			read1 = alignmentFolder + "/" + file
			outputPrefix = outputSubfolder + "/" + sample
			#you can also add '—estimate-rspd –calc-pme –calc-ci'
			#however, this will increase the run-time and code doesn't currently use those columns
			command = "/opt/RSEM/rsem-calculate-expression --bam --no-bam-output -p " + str(threads) + "  " + read1 + "  " + rsemIndex + " " + outputPrefix
			os.system(command)
