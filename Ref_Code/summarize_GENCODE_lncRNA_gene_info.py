import sys
import re
import os

gtf = "hg19/gencode.v24lift37.long_noncoding_RNAs.gtf"
geneTable = "hg19/hg19_gene_info.txt"
txTable = "hg19/hg19_transcript_info.txt"

geneOut = open(geneTable, 'w')
text = "symbol\tgene.id\tchr\tstart\tend\tstrand\ttx.max.length\n"
geneOut.write(text)

txOut = open(txTable, 'w')
text = "transcript.id\tsymbol\tchr\tstart\tend\tstrand\tlength\n"
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
			geneID = "NA"
			txID = "NA"
			geneSymbol = "NA"
			for value in annInfo:
				annResult1 = re.search("gene_id (\".*\")", value)
				annResult2 = re.search("transcript_id (\".*\")", value)
				annResult3 = re.search("gene_name (\".*\")", value)
				if annResult1:
					geneID = annResult1.group(1)
				if annResult2:
					txID = annResult2.group(1)
				if annResult3:
					geneSymbol = annResult3.group(1)
					
			if (geneID == "NA") or (txID == "NA") or (geneSymbol == "NA"):
				print "Issue parsing " + line
				print "geneID : " + geneID
				print "txID : " + txID
				print "symbol : " + geneSymbol
			
			geneHash[txID] = geneSymbol + "\t" + geneID
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

geneIDhash = {}
geneChrHash = {}
geneStartHash = {}
geneStopHash = {}
geneStrandHash = {}
geneLengthHash = {}

for txID in geneHash:
	geneText = geneHash[txID]
	geneInfo = geneText.split("\t")
	geneSymbol = geneInfo[0]
	geneID = geneInfo[1]
	geneIDhash[geneSymbol] = geneID
	
	txChr = txChrHash[txID]
	txStart = txStartHash[txID]
	txStop = txStopHash[txID]
	txStrand = txStrandHash[txID]
	txLength = txLengthHash[txID]
	
	geneChrHash[geneSymbol] = txChr
	geneStrandHash[geneSymbol] = txStrand
	
	if geneSymbol in geneStartHash:
		if txStart < geneStartHash[geneSymbol]:
			geneStartHash[geneSymbol] = txStart
	else:
		geneStartHash[geneSymbol] = txStart

	if geneSymbol in geneStopHash:
		if txStop > geneStopHash[geneSymbol]:
			geneStopHash[geneSymbol] = txStop
	else:
		geneStopHash[geneSymbol] = txStop

	if geneSymbol in geneLengthHash:
		if txLength > geneLengthHash[geneSymbol]:
			geneLengthHash[geneSymbol] = txLength
	else:
		geneLengthHash[geneSymbol] = txLength
		
	text = txID + "\t" + geneSymbol + "\t" + txChr + "\t" + str(txStart) + "\t" + str(txStop) + "\t" + txStrand + "\t" + str(txLength) + "\n"
	txOut.write(text)
	
for geneSymbol in geneIDhash:
	geneID = geneIDhash[geneSymbol]
	
	geneChr = geneChrHash[geneSymbol]
	geneStart = geneStartHash[geneSymbol]
	geneStop = geneStopHash[geneSymbol]
	geneStrand = geneStrandHash[geneSymbol]
	geneLength = geneLengthHash[geneSymbol]
		
	text = geneSymbol + "\t" + geneID + "\t" + geneChr + "\t" + str(geneStart) + "\t" + str(geneStop) + "\t" + geneStrand + "\t" + str(geneLength) + "\n"
	geneOut.write(text)