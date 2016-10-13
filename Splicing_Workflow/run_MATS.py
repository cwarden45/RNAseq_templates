import sys
import re
import os

alignmentFolder = "/path/to/TopHat_Alignment"

matsFolder = "/path/to/rMATS.3.2.5/"
matsGtf = matsFolder + "gtf/Homo_sapiens.knownGene.hg19.gtf"

#see 'MATS_extract_gene_acc_from_GTF' template for creation of UniProt and RefSeq ID mappings

#at time, this created mapping has 'From' (UniProt) and 'To' (gene symbol) columns
uniprotConversionFile = matsFolder + "gtf/hg19_UniProt_to_Gene_Symbol.txt"

#use : https://biodbnet-abcc.ncifcrf.gov/db/db2db.php
#	although header is different same fuction used to parse both UniProt and RefSeq files
refseqConversionFile = matsFolder + "gtf/hg19_RefSeq_to_Gene_Symbol.txt"

#can be 'paired' or 'single'
libType = ""
readLength = ""

#can be 'P' (for paired) or 'U' (for unpaired)
analysisType = ""

#can be 'fr-unstranded' (for unstranded), 'fr-firststrand' ('reverse' for HT-Seq, typical for stranded Illumina data) or 'fr-secondstrand' ('yes' for HT-Seq)
strand = ""

resultFolder = ""
#can add replicates with comma-separated .bam files (line 40)
CONTROLbam = os.path.join(alignmentFolder, ".bam")
TRTbam = os.path.join(alignmentFolder, ".bam")

fdr_cutoff = 0.05
inclusion_cutoff =0.1

### for 1-vs-1 comparison, modify code above this line ##
#using `-c 0.1` to require a 10% splicing difference, but you'll still need to filter output table by IncLevelDifference and FDR
#enable novel splice site detection with `novelSS 1`
#make sure you use python 2.7

command = "python " + matsFolder + "RNASeq-MATS.py -b1 " + TRTbam + " -b2 " + CONTROLbam+ " -gtf " + matsGtf + " -o " + resultFolder + " -c 0.1 -novelSS 1 -t " +libType+" -len " + readLength + " -analysis " + analysisType + " -libType " + strand
#os.system(command)

matsOutputFolder = resultFolder + "/MATS_output"
reformatFolder = resultFolder + "/MATS_output_reformat"
if not os.path.exists(reformatFolder):
	command = "mkdir " + reformatFolder
	os.system(command)

if True:
	#run annotation functions - requires prepaparation of UniProt ID mapping and RefSeq access mapping
	#so, make sure you have created / downloaded `uniprotConversionFile` and `refseqConversionFile`
	def makeGeneHash(inputfile):
		dict = {}

		inHandle = open(inputfile)
		line = inHandle.readline()

		lineCount = 0

		while line:
			lineCount += 1

			
			if(lineCount > 1):
				line = line.replace("\n","")
				line = line.replace("\r","")
		
				lineInfo = line.split("\t")
		
				matsID = lineInfo[0]
				geneName = lineInfo[1]
				
				dict[matsID] = geneName
			
			line = inHandle.readline()

		return dict
		
	geneDict = makeGeneHash(uniprotConversionFile)
	geneDict2 = makeGeneHash(refseqConversionFile)
	
	files = os.listdir(matsOutputFolder)

	for file in files:
		print file
		
		upCount = 0
		downCount = 0
		
		fdrIndex = 19
		inclusionIndex = 22
		
		MXEresult = re.search("^MXE",file)
		if MXEresult:
			fdrIndex = 21
			inclusionIndex = 24
		
		inputfile = matsOutputFolder + "/" + file
		outputfile = reformatFolder + "/" + file
		outputfile = re.sub(".txt$",".csv",outputfile)
		
		inHandle = open(inputfile)
		line = inHandle.readline()

		outHandle = open(outputfile, "w")

		lineCount = 0
		
		while line:
			line = line.replace("\n","")
			line = line.replace("\r","")
		
			lineCount += 1
			
			if lineCount == 1:
				text = re.sub("\t",",",line) + ",Status\n"
				outHandle.write(text)
			else:	
				lineInfo = line.split("\t")
			
				gene = lineInfo[2]
				gene = gene.replace("\"","")
				#gene = gene.replace("-\d+$","")
				
				if gene in geneDict:
					lineInfo[2] = geneDict[gene]
				elif gene in geneDict2:
					lineInfo[2] = geneDict2[gene]
				
				status = "No Change"
				
				if (float(lineInfo[inclusionIndex]) >= inclusion_cutoff) and (float(lineInfo[fdrIndex]) < fdr_cutoff):
					status = "Increased Inclusion"
					upCount += 1
				if (float(lineInfo[inclusionIndex]) <= -inclusion_cutoff) and (float(lineInfo[fdrIndex]) < fdr_cutoff):
					status = "Decreased Inclusion"
					downCount += 1
					
				text = ",".join(lineInfo) + "," + status + "\n"
				outHandle.write(text)
		

			line = inHandle.readline()
		print "Increased Inclusion : " + str(upCount)
		print "Decreased Inclusion : " + str(downCount)