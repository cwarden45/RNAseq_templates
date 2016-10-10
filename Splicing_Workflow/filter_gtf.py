import sys
import re
import os

gencodeGTF = "human/hg19/gencode.v19.annotation.gtf"
reducedGTF = "human/hg19/gencode.v19.annotation.FILTERED.gtf"
geneTable = "human/RSEM_Bowtie_Index/gene_info.txt"

#kept gene hash
geneHash = {}

inHandle = open(geneTable)
line = inHandle.readline()

lineCount = 0

while line:
	line = re.sub("\n","",line)
	line = re.sub("\"","",line)
	line = re.sub("\r","",line)
	lineInfo = line.split("\t")
	
	lineCount += 1
	
	if lineCount > 1:
		transcript = lineInfo[0]
		geneHash[transcript]=1
	line = inHandle.readline()

inHandle.close()
print "Select exons from " + str(len(geneHash.keys())) + " transcripts"

#filter GTF

outHandle = open(reducedGTF,'w')

inHandle = open(gencodeGTF)
line = inHandle.readline()

while line:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	lineInfo = line.split("\t")
	
	commentResult = re.search("^#",line)
	
	if commentResult:
		text = line + "\n"
		outHandle.write(text)
	else:
		#could really use same filters as .fa file to create filtered .gtf, but I want to make sure the two tables will match
		#all samples in filtered table have transcript IDs that start with ENST (not ENSTG or ENSTGR)
		transcriptResult = re.search("transcript_id \"(ENS[T|TR]\d+.\d+)\";",lineInfo[8])
		if transcriptResult:
			transcript = transcriptResult.group(1)
			if transcript in geneHash:
				text = line + "\n"
				outHandle.write(text)

	line = inHandle.readline()
	