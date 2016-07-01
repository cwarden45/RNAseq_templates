import sys
import re
import os

threads = 2
ref = "/path/to/RSEM_Bowtie_index/Bowtie_gencode_v24"
alignmentFolder = "Bowtie_alignment"
readsFolder = "Reads"

finishedSamples = ()
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
								
			bowtieSam = outputSubfolder + "/alignment.sam";	
			command = "/opt/bowtie/bowtie" + " -p " + str(threads) + " -S " + " " + ref + " " + read1 + " > " + bowtieSam
			os.system(command)
						
			bowtieBam = outputSubfolder + "/alignment.bam";																			
			command = "/opt/samtools-1.3/samtools view -b " + bowtieSam + " > " + bowtieBam
			os.system(command)
								
			userBam = alignmentFolder + "/" + sample + ".bam";
			command = "/opt/samtools-1.3/samtools sort " + bowtieBam + " -o " + userBam
			os.system(command)
			
			command = "/opt/samtools-1.3/samtools index " + userBam
			os.system(command)
			
			command = "rm -R " + outputSubfolder
			os.system(command)
