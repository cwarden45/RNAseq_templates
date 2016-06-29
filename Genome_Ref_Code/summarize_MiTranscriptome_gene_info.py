import sys
import re
import os

gtf = "hg19/hg19_cancer_specific.gtf"
txTable = "hg19/hg19_gene_info.txt"

txOut = open(txTable, 'w')
text = "symbol\ttranscript.id\tchr\tstart\tend\tstrand\ttx.max.length\n"
txOut.write(text)

inHandle = open(gtf)
lines = inHandle.readlines()

lineCount = 0

geneHash = {}
txChrHash = {}
txStartHash = {}
txStopHash = {}
txStrandHash = {}
txLengthHash = {}

for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	lineCount +=1
	
	ignoreResult = re.search("^#",line)
	
	if not ignoreResult:
		lineInfo = line.split("\t")
		chr = lineInfo[0]
		annType = lineInfo[2]
		start = int(lineInfo[3])
		stop = int(lineInfo[4])
		symbol = lineInfo[4]
		strand = lineInfo[6]
		
		if annType == "exon":
			annText = lineInfo[8]
			annInfo = annText.split("; ")
			txID = "NA"
			geneSymbol = "NA"
			for value in annInfo:
				annResult1 = re.search("gene_id (\".*\")", value)
				annResult2 = re.search("transcript_id (\".*\")", value)
				if annResult1:
					geneSymbol = annResult1.group(1)
				if annResult2:
					txID = annResult2.group(1)
					
			if (txID == "NA") or (geneSymbol == "NA"):
				print "Issue parsing " + line
				print "txID : " + txID
				print "symbol : " + geneSymbol
			
			geneHash[txID] = geneSymbol
			txChrHash[txID] = chr
			txStrandHash[txID] = strand
			
			if txID in txLengthHash:
				if start < txStartHash[txID]:
					txStartHash[txID] = start
				
				if stop > txStopHash[txID]:
					txStopHash[txID] = stop
					
				txLengthHash[txID] = txLengthHash[txID] + stop - start
			else:
				txStartHash[txID] = start
				txStopHash[txID] = stop
				txLengthHash[txID] = stop - start			

for txID in geneHash:
	geneSymbol = geneHash[txID]
	
	txChr = txChrHash[txID]
	txStart = txStartHash[txID]
	txStop = txStopHash[txID]
	txStrand = txStrandHash[txID]
	txLength = txLengthHash[txID]
		
	text = geneSymbol + "\t" + txID + "\t" + txChr + "\t" + str(txStart) + "\t" + str(txStop) + "\t" + txStrand + "\t" + str(txLength) + "\n"
	txOut.write(text)