import sys
import re
import os

parameterFile = "parameters.txt"

alignmentFolder = ""
bed_file = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	
	if param == "Alignment_Folder_MAC":
		alignmentFolder = value
		
	if param == "RSeQC_bed_MAC":
		bed_file = value	

fileResults = os.listdir(alignmentFolder)

finishedSamples = ()

for file in fileResults:
	result = re.search(".bam$",file)
	fullPath = os.path.join(alignmentFolder, file)
	
	if result:
		sample = re.sub(".bam$","",file)
		sortResult = re.search(".name.sort.bam",file)
		if (sample not in finishedSamples) and (not sortResult):
			print sample
			subfolder = alignmentFolder + "/" + sample
			
			print "Determining Strand for Housekeeping Genes"
			strandStat = subfolder + "/" + sample + "_infer_strand.txt"
			command = "infer_experiment.py -r " + bed_file + " -i " + fullPath + " > " + strandStat
			os.system(command)

			print "Calculating TIN scores"
			#skip since output file not created
			tinStat = subfolder + "/" + sample + "_tin_stat.txt"
			command = "tin.py -r " + bed_file + " -i " + fullPath + " > " + tinStat
			#os.system(command)

		
