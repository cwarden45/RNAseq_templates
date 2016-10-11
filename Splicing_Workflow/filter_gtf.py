import sys
import re
import os

gencodeGTF = "human/hg19/gencode.v19.annotation.gtf"
reducedGTF = "human/hg19/gencode.v19.annotation.FILTERED.gtf"
geneTable = "human/hg19/gencode.v19.gene_info.txt"

fullTableHandle = open(geneTable,'w')
text = "Ensembl.TranscriptID\tEnsembl.GeneID\tTranscript.Name\tGene.Name\tTranscript.Type\n"
fullTableHandle.write(text)

outHandle = open(reducedGTF,'w')

inHandle = open(gencodeGTF)
line = inHandle.readline()

transcriptHash = {}
geneHash = {}

while line:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	lineInfo = line.split("\t")
	
	commentResult = re.search("^#",line)
	
	if commentResult:
		text = line + "\n"
		outHandle.write(text)
	else:
		annotationType = lineInfo[2]
		
		if (annotationType == "gene") or (annotationType == "exon")or (annotationType == "transcript")or (annotationType == "CDS"):
			annText = lineInfo[8]
			annInfo = annText.split(";")
			
			ensGene = ""
			ensTrans = ""
			geneName = ""
			transName = ""
			transType = ""
			
			for ann in annInfo:
				geneIDResult = re.search("gene_id (.*)",ann)
				transIDResult = re.search("transcript_id (.*)",ann)
				geneNameResult = re.search("gene_name (.*)",ann)
				transNameResult = re.search("transcript_name (.*)",ann)
				transTypeResult = re.search("transcript_type (.*)",ann)
			
				if geneIDResult:
					ensGene = geneIDResult.group(1)
				elif transIDResult:
					ensTrans = transIDResult.group(1)
				elif geneNameResult:
					geneName = geneNameResult.group(1)
				elif transNameResult:
					transName = transNameResult.group(1)				
				elif transTypeResult :
					transType = transTypeResult.group(1)
					
			if (ensGene == "") or (ensTrans == "") or (geneName == "")or (transName == "")or (transType == ""):
				print "Revise code for parsing: " + annText
				sys.exit()
			
			#inclusive search
			pcResult = re.search("protein_coding",transType)
			lncResult = re.search("lncRNA",transType)
			miRNAResult = re.search("miRNA",transType)
			scarnaResult = re.search("scaRNA",transType)
			snrnaResult = re.search("snRNA",transType)
			vaultrnaResult = re.search("vaultRNA",transType)
			ribozymeResult = re.search("ribozyme",transType)
			snornaResult = re.search("snoRNA",transType)
			lincrnaResult = re.search("lincRNA",transType)
			
			if pcResult or lncResult or miRNAResult or scarnaResult or vaultrnaResult or ribozymeResult or snornaResult or lincrnaResult:
				text = line + "\n"
				outHandle.write(text)
				
				geneHash[geneName] = 1
				
				if ensTrans not in transcriptHash:
					transcriptHash[ensTrans]=1
					
					text = ensTrans + "\t"+ensGene+"\t"+transName+"\t"+geneName+"\t"+transType+"\n"
					fullTableHandle.write(text)
		elif (annotationType != "start_codon") and (annotationType != "stop_codon")and (annotationType != "UTR")and (annotationType != "Selenocysteine"):
			print "Keep or drop annotation type? : " + annotationType
			sys.exit()

	line = inHandle.readline()
	
print str(len(transcriptHash.keys())) + " transcripts considered in final transcript reference"
print str(len(geneHash.keys())) + " passing genes (includes lincRNAs)"
	