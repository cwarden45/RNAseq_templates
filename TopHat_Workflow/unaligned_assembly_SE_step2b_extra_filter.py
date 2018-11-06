import sys
import re
import os

finishedSamples = ()

readFolder = "Downsampled_Reads/Trimmed_Reads"	
alignmentFolder = "Oases_Assemblies_Trim"
filteredFolder = "BLAST_BWA_filtered"

threads = 4
bwaRef = "BLAST_BWA_filtered/BLAST_hits_1st_pass_2each_UNIQUE.fasta"

command = "bwa index " + bwaRef
os.system(command)

fileResults = os.listdir(readFolder)

for file in fileResults:
	result = re.search("(.*)_R1_cutadapt.fastq$",file)
	fullPath = os.path.join(readFolder, file)
	
	if result:
		sample = result.group(1)
		if sample not in finishedSamples:
			print sample
			subfolder = alignmentFolder + "/" + sample
			subfolderF = filteredFolder + "/" + sample
			command = "mkdir " + subfolderF
			os.system(command)
			
			assembly = subfolder + "/transcripts.fa"

			alnSam = subfolderF + "/aligned.sam"
			command = "bwa mem -t "+ str(threads) + " " + bwaRef + " " + assembly + " > " + alnSam + "\n"
			os.system(command)

			assemblyF = subfolderF + "/transcripts.fa"
			command = "samtools fasta -f 4 "+ alnSam + " > " + assemblyF + "\n"
			os.system(command)
			
			bowtieIndex = subfolderF + "/" + sample
			command = "bowtie2-build " + assemblyF + " " + bowtieIndex
			os.system(command)
			
			alnSam = subfolderF + "/" + sample + ".sam"
			command = "bowtie2 -x " + bowtieIndex + " -U " + fullPath + " -S " + alnSam
			os.system(command)
			
			alnBam = subfolderF + "/" + sample + ".bam"
			command = "samtools sort " + alnSam + " -o " + alnBam
			os.system(command)

			command = "rm " + alnSam
			os.system(command)
			
			statFile = subfolderF + "/" + sample + "_flagstat.txt"
			command = "samtools flagstat " + alnBam + " > " + statFile
			os.system(command)
			
			command = "/opt/express-1.5.1-linux_x86_64/express --no-bias-correct " + assemblyF + " " + alnBam + " -o " + subfolderF
			os.system(command)