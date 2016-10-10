import sys
import re
import os

parameterFile = "parameters.txt"
finishedSamples = []

rawCountFolder = ""
gtfFile = ""
alignmentFolder = ""
javaMem = ""
pairing = ""
strand = ""
qortsJar = "/opt/QoRTs/QoRTs_1.1.8/QoRTs.jar"

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]

	if param == "strand":
		strand = value
	
	if param == "pairing":
		pairing = value
	
	if param == "QoRTs_Count_Folder":
		rawCountFolder = value

	if param == "Alignment_Folder":
		alignmentFolder = value
		
	if param == "GENCODE_GTF":
		gtfFile = value		

	if param == "Java_Mem":
		javaMem = value	

if (alignmentFolder == "") or (alignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder'!"
	sys.exit()
		
if (rawCountFolder == "") or (rawCountFolder == "[required]"):
	print "Need to enter a value for 'QoRTs_Count_Folder'!"
	sys.exit()
	
if (gtfFile == "") or (gtfFile == "[required]"):
	print "Need to enter a value for 'GENCODE_GTF'!"
	sys.exit()
	
if (javaMem == "") or (javaMem == "[required]"):
	print "Need to enter a value for 'Java_Mem'!"
	sys.exit()

if (pairing == "") or (pairing == "[required]"):
	print "Need to enter a value for 'pairing'!"
	sys.exit()
	
if (strand == "") or (strand == "[required]"):
	print "Need to enter a value for 'strand'!"
	sys.exit()
	
fileResults = os.listdir(alignmentFolder)

for file in fileResults:
	result = re.search(".bam$",file)
	fullPath = os.path.join(alignmentFolder, file)
	
	if result:
		sample = re.sub(".bam$","",file)
		
		if sample not in finishedSamples:
			print sample
			resultFolder = rawCountFolder + "/" + sample
			command = "mkdir " + resultFolder
			os.system(command)
			
			singleEndFlag = ""
			if pairing == "SE":
				singleEndFlag = "--singleEnded "
			
			strandedFlag = ""
			if strand == "reverse":
				strandedFlag = "--stranded "
			if strand == "yes":
				strandedFlag = "--stranded --stranded_fr_secondstrand "
				
			#omit --keepMultiMapped option listed in manual
			command = "java -Xmx"+javaMem+" -jar "+qortsJar+" QC "+singleEndFlag+strandedFlag+"--runFunctions writeKnownSplices,writeNovelSplices,writeSpliceExon " + fullPath + " " + gtfFile + " " + resultFolder
			os.system(command)
