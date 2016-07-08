import sys
import re
import os

headerPrefix = "hg19_ct_MiTranscriptome_5659_"

tableBrowserFa = "mitranscriptome.fa"
rsemFa = "RSEM_Bowtie_Index/mitranscriptome.fa"
geneTable = "RSEM_Bowtie_Index/gene_info.txt"
rsemMap = "RSEM_Bowtie_Index/transcript_to_gene.txt"

newFaHandle = open(rsemFa,'w')
fullTableHandle = open(geneTable,'w')
text = "Ensembl.TranscriptID\tTranscript.Name\tGene.Name\tNCBI.GeneID\tGene.Type\n"
fullTableHandle.write(text)
rsemHandle = open(rsemMap,'w')

inHandle = open(tableBrowserFa)

geneHash = {}
transcriptCount = 0

line = inHandle.readline()

writeFlag = 1

while line:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	lineInfo = line.split("\t")
	
	headerResult = re.search("^>" + headerPrefix + "(\S+)",line)
	
	if headerResult:
		transcriptCount += 1
		transcriptName =  headerResult.group(1)
		
		multiGeneResult = re.search("(.*)\.\d+$",transcriptName)
		if multiGeneResult:
			geneName =  multiGeneResult.group(1)
		else:
			geneName = transcriptName
		geneHash[geneName]=1	
		
		text = ">" + transcriptName + "\n"
		newFaHandle.write(text)

		text = geneName + "\t" + transcriptName +  "\n"
		rsemHandle.write(text)
		
		text = "\"" + transcriptName + "\"\t\"" + transcriptName + "\"\t\"" + geneName + "\"\tNA\tlncRNA\n"
		fullTableHandle.write(text)
	else:
		text = line + "\n"
		newFaHandle.write(text)		
	line = inHandle.readline()
	
print str(transcriptCount) + " MiTranscriptome transcripts"
print str(len(geneHash.keys())) + " MiTranscriptome genes"