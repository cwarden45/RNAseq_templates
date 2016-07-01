import sys
import re
import os


mapFile = "/path/to/RSEM_Bowtie_index/transcript_to_gene.txt"
refFa = "/path/to/RSEM_Bowtie_index/gencode.v24.transcriptID.fa"
rsemIndex = "/path/to/RSEM_Bowtie_index/Bowtie_gencode_v24"
command = "/opt/RSEM/rsem-prepare-reference --bowtie --bowtie-path /opt/bowtie --transcript-to-gene-map " + mapFile + " " + refFa + " " + rsemIndex
os.system(command)

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
