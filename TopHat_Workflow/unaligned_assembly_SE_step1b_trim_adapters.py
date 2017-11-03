import os
import sys
import re
from Bio.Seq import Seq

readsFolder = "Downsampled_Reads"
finishedSamples = ()


cutadaptFolder = readsFolder + "/Trimmed_Reads"
command = "mkdir " + cutadaptFolder
os.system(command)

fileResults = os.listdir(readsFolder)

jobCount = 0
for file in fileResults:
	result = re.search("(.*)_R1.fastq$",file)
	
	if result:
		jobCount += 1
		sample = result.group(1)
		if sample not in finishedSamples:
			print sample
			
			read1 = readsFolder + "/" + file

			trim1 = cutadaptFolder + "/" + sample + "_R1_temp.fastq"

			Fadapter = "TACACTCTTTCCCTACACGACGCTCTTCCGATCT"
			Radapter = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"			
			seqObj = Seq(Radapter)
			revcomR =  str(seqObj.reverse_complement())
			seqObj = Seq(Fadapter)
			revcomF =  str(seqObj.reverse_complement())
			
			command = "cutadapt -a "+Radapter+" -g "+revcomF+" -m 40 --max-n 0 " +  read1 + " > " + trim1
			os.system(command)
			
			trim2 = cutadaptFolder + "/" + sample + "_R1_cutadapt.fastq"
			command = "cutadapt -a "+Fadapter+" -g "+revcomR+" -m 40 --max-n 0 " +  trim1 + " > " + trim2
			os.system(command)
			
			command = "rm  " +  trim1
			os.system(command)