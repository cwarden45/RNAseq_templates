import sys
import re
import os

parameterFile = "parameters.txt"

alignmentFolder = ""
bed_file = ""
sample_description = ""

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

	if param == "sample_description_file":
		sample_description = value	
		
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
			command = "tin.py -r " + bed_file + " -i " + fullPath
			os.system(command)
			
			tinOut1 = sample + ".summary.txt"
			command = "mv " + tinOut1 + " " + subfolder + "/" + tinOut1
			os.system(command)
			
			tinOut2 = sample + ".tin.xls"
			command = "mv " + tinOut2 + " " + subfolder + "/" + tinOut2
			os.system(command)

plotBams = ""
lineCount = 0
sampleIDindex = -1

inHandle = open(sample_description)
line = inHandle.readline()

while line:
	line = re.sub("\r","",line)
	line = re.sub("\n","",line)
	lineInfo = line.split("\t")
	
	lineCount +=1
	
	if lineCount ==1:
		for i in range(0, len(lineInfo)):
			if lineInfo[i] == "sampleID":
				sampleIDindex = i
				
		if sampleIDindex == -1:
			print "There must be a column in the sample description file called 'sampleID' to define ordered .bam files!"
			sys.exit()
	elif lineCount == 2:
		plotBams = alignmentFolder + "/" + lineInfo[sampleIDindex] + ".bam"
	else:
		plotBams = plotBams + "," + alignmentFolder + "/" + lineInfo[sampleIDindex] + ".bam"
	line = inHandle.readline()

command = "geneBody_coverage.py -r " + bed_file + " -i " + plotBams + " -o gene_coverage"
os.system(command)	
