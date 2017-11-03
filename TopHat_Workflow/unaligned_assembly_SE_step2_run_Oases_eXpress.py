import sys
import re
import os

unfinishedSamples = ("","")

readFolder = "Downsampled_Reads"	
alignmentFolder = "Oases_Assemblies"

#readFolder = "Downsampled_Reads/Trimmed_Reads"	
#alignmentFolder = "Oases_Assemblies_Trim"

fileResults = os.listdir(readFolder)

for file in fileResults:
	result = re.search("(.*)_R1.fastq$",file)
	#result = re.search("(.*)_R1_cutadapt.fastq$",file)
	fullPath = os.path.join(readFolder, file)
	
	if result:
		sample = result.group(1)
		if sample in unfinishedSamples:
			print sample
			subfolder = alignmentFolder + "/" + sample
			command = "mkdir " + subfolder
			os.system(command)
			
			#add -strand_specific
			command = "/opt/velvet_1.2.10/velveth " + subfolder + " 21 -strand_specific -short -fastq " + fullPath
			os.system(command)
			
			command = "/opt/velvet_1.2.10/velvetg " + subfolder+ " -cov_cutoff 10 -read_trkg yes"
			os.system(command)		

			command = "/opt/oases_0.2.8/oases " + subfolder
			os.system(command)

			assembly = subfolder + "/transcripts.fa"
			bowtieIndex = subfolder + "/" + sample
			command = "bowtie2-build " + assembly + " " + bowtieIndex
			os.system(command)
			
			alnSam = subfolder + "/" + sample + ".sam"
			command = "bowtie2 -x " + bowtieIndex + " -U " + fullPath + " -S " + alnSam
			os.system(command)
			
			alnBam = subfolder + "/" + sample + ".bam"
			command = "samtools sort " + alnSam + " -o " + alnBam
			os.system(command)

			command = "rm " + alnSam
			os.system(command)

			command = "rm " + subfolder + "/Graph2"
			os.system(command)
			
			command = "rm " + subfolder + "/LastGraph"
			os.system(command)			

			command = "rm " + subfolder + "/Sequences"
			os.system(command)

			command = "rm " + subfolder + "/Roadmaps"
			os.system(command)

			command = "rm " + subfolder + "/PreGraph"
			os.system(command)
			
			statFile = subfolder + "/" + sample + "_flagstat.txt"
			command = "samtools flagstat " + alnBam + " > " + statFile
			os.system(command)
			
			command = "/opt/express-1.5.1-linux_x86_64/express --no-bias-correct " + assembly + " " + alnBam + " -o " + subfolder
			os.system(command)