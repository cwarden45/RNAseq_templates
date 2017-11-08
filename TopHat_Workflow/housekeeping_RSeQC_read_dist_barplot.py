import sys
import re
import os
import numpy as np
import matplotlib.pyplot as plt

figure_out_file = "housekeeping_gene_read_distribution.png"
x_shift=0.4
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
		
	if param == "RSeQC_bed_ngs":
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

			print "Calculate Read Distribution for Housekeeping Genes"
			readStat = subfolder + "/" + sample + "_read_distribution.txt"
			command = "read_distribution.py -r " + bed_file + " -i " + fullPath + " > " + readStat
			os.system(command)

plotBams = ""
lineCount = 0
sampleIDindex = -1
sampleLabelindex = -1

inHandle = open(sample_description)
line = inHandle.readline()

sampleLabels = []
percentCDS = []
percentUTR5 = []
percentUTR3 = []
percentIntron = []

while line:
	line = re.sub("\r","",line)
	line = re.sub("\n","",line)
	lineInfo = line.split("\t")
	
	lineCount +=1
	
	if lineCount ==1:
		for i in range(0, len(lineInfo)):
			if lineInfo[i] == "sampleID":
				sampleIDindex = i
			if lineInfo[i] == "userID":
				sampleLabelindex = i

		if sampleIDindex == -1:
			print "There must be a column in the sample description file called 'sampleID' to map outputfile!"
			sys.exit()
				
		if sampleLabelindex == -1:
			print "There must be a column in the sample description file called 'userID' for labeling plots!"
			sys.exit()
	else:
		sampleID = lineInfo[sampleIDindex]
		sampleLabel = lineInfo[sampleLabelindex]
		sampleLabels.append(sampleLabel)
		readStat = alignmentFolder + "/" + sampleID + "/" + sampleID + "_read_distribution.txt"

		sampleCDS = 0
		sampleUTR3 = 0
		sampleUTR5 = 0
		sampleIntron = 0
		
		statHandle = open(readStat)
		stat_lines = statHandle.readlines()
					
		for stat_line in stat_lines:
			stat_line = re.sub("\n","",stat_line)
			stat_line = re.sub("\r","",stat_line)
			if(re.search("^CDS_Exons",stat_line)):
				statResult = re.search("^CDS_Exons\s+\d+\s+(\d+)",stat_line)
				sampleCDS = int(statResult.group(1))
			elif(re.search("^5'UTR_Exons",stat_line)):
				statResult = re.search("^5'UTR_Exons\s+\d+\s+(\d+)",stat_line)
				sampleUTR5 = int(statResult.group(1))
			elif(re.search("^3'UTR_Exons",stat_line)):
				statResult = re.search("^3'UTR_Exons\s+\d+\s+(\d+)",stat_line)
				sampleUTR3 = int(statResult.group(1))
			elif(re.search("^Introns",stat_line)):
				statResult = re.search("^Introns\s+\d+\s+(\d+)",stat_line)
				sampleIntron = int(statResult.group(1))
		statHandle.close()
		
		tagSum = sampleCDS + sampleUTR3 + sampleUTR5 + sampleIntron
		percentCDS.append(float(sampleCDS) / float(tagSum))
		percentUTR5.append(float(sampleUTR5) / float(tagSum))
		percentUTR3.append(float(sampleUTR3) / float(tagSum))
		percentIntron.append(float(sampleIntron) / float(tagSum))
	line = inHandle.readline()

ind = np.arange(len(sampleLabels))
pI = plt.bar(ind, percentIntron, color="red")
prevSum = percentIntron
p5 = plt.bar(ind, percentUTR5, color="purple", bottom=prevSum)
prevSum = np.add(prevSum,percentUTR5)

p3 = plt.bar(ind, percentUTR3, color="green", bottom=prevSum)
prevSum = np.add(prevSum,percentUTR3)
pC = plt.bar(ind, percentCDS, color="blue", bottom=prevSum)

plt.ylabel('Relative Percent')
plt.xticks(np.add(ind,x_shift), sampleLabels, rotation=90)
plt.yticks(np.arange(0, 1.01, 0.1))

plt.legend((pC[0], p3[0], p5[0], pI[0]), ('CDS', '3`UTR','5`UTR','Intron'), bbox_to_anchor=(0., 1.02, 1., .102), ncol=4)

plt.savefig(figure_out_file)
