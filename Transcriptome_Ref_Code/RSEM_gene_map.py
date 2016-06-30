import sys
import re
import os

#gencodeFa = "mouse/gencode.vM9.transcripts.fa"
#rsemFa = "mouse/RSEM_Bowtie_Index/gencode.vM9.transcriptID.fa"
#geneTable = "mouse/RSEM_Bowtie_Index/gene_info.txt"
#rsemMap = "mouse/RSEM_Bowtie_Index/transcript_to_gene.txt"

gencodeFa = "human/gencode.v24.transcripts.fa"
rsemFa = "human/RSEM_Bowtie_Index/gencode.v24.transcriptID.fa"
geneTable = "human/RSEM_Bowtie_Index/gene_info.txt"
rsemMap = "human/RSEM_Bowtie_Index/transcript_to_gene.txt"

newFaHandle = open(rsemFa,'w')
fullTableHandle = open(geneTable,'w')
text = "Ensembl.TranscriptID\tEnsembl.GeneID\tHAVANA.TranscriptID\tHAVANA.GeneID\tTranscript.Name\tGene.Name\tNCBI.GeneID\tGene.Type\n"
fullTableHandle.write(text)
rsemHandle = open(rsemMap,'w')

inHandle = open(gencodeFa)

typeHash = {}
geneHash = {}
transcriptCount = 0

line = inHandle.readline()

writeFlag = 1

while line:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	lineInfo = line.split("\t")
	
	headerResult = re.search("^>(.*)",line)
	
	if headerResult:
		header =  headerResult.group(1)
		headerInfo = header.split("|")
		transcriptID = headerInfo[0]
		geneSymbol = headerInfo[5]
		geneType = headerInfo[7]
		
		#inclusive search
		pcResult = re.search("protein_coding",geneType)
		lncResult = re.search("lncRNA",geneType)
		miRNAResult = re.search("miRNA",geneType)
		scarnaResult = re.search("scaRNA",geneType)
		snrnaResult = re.search("snRNA",geneType)
		vaultrnaResult = re.search("vaultRNA",geneType)
		ribozymeResult = re.search("ribozyme",geneType)
		snornaResult = re.search("snoRNA",geneType)
		lincrnaResult = re.search("lincRNA",geneType)
	
		#if pcResult or lncResult or miRNAResult:	
		if pcResult or lncResult or miRNAResult or scarnaResult or snrnaResult or vaultrnaResult or ribozymeResult or snornaResult or lincrnaResult:
			writeFlag = 1
			transcriptCount += 1
			typeHash[geneType] = 1
			geneHash[geneSymbol] = 1

			text = ">" + transcriptID + "\n"
			newFaHandle.write(text)

			text = geneSymbol + "\t" + transcriptID +  "\n"
			rsemHandle.write(text)
		
			text = "\"" + "\"\t\"".join(headerInfo[0:8]) + "\"\n"
			fullTableHandle.write(text)
		else:
			writeFlag = 0
		
		#exclusive search
		pseudoResult = re.search("pseudogene",geneType)
		rRNAResult = re.search("rRNA",geneType)
		tecResult = re.search("TEC",geneType)
		nmdResult = re.search("nonsense_mediated_decay",geneType)
		nmsResult = re.search("non_stop_decay",geneType)
		intronResult = re.search("retained_intron",geneType)
		trnaResult = re.search("retained_intron",geneType)
		asResult = re.search("antisense",geneType)
		trnaResult = re.search("tRNA",geneType)
		otherNcResult = re.search("^non_coding",geneType)
		srnaResult = re.search("sRNA",geneType)
		miscrnaResult = re.search("misc_RNA",geneType)
	
		#if pseudoResult or rRNAResult or tecResult or nmdResult or nmsResult or intronResult or asResult or trnaResult or otherNcResult or srnaResult or miscrnaResult:
		#	writeFlag = 0
		#else:
		#	writeFlag = 1
		#	transcriptCount += 1
		#	typeHash[geneType] = 1
		#	geneHash[geneSymbol] = 1

		#	text = ">" + transcriptID + "\n"
		#	newFaHandle.write(text)

		#	text = geneSymbol + "\t" + transcriptID +  "\n"
		#	rsemHandle.write(text)
		
		#	text = "\"" + "\"\t\"".join(headerInfo[0:8]) + "\"\n"
		#	fullTableHandle.write(text)
	elif writeFlag == 1:
		text = line + "\n"
		newFaHandle.write(text)
	line = inHandle.readline()
	
print " , ".join(typeHash.keys()) + " considered in final transcript reference"
print str(transcriptCount) + " passing transcripts"
print str(len(geneHash.keys())) + " passing genes (includes lincRNAs)"
