import sys
import re
import os

#gtfFile = "Mus_musculus.knownGene.mm9.gtf"
#geneList = "mm9_UniProt_IDs.txt"

gtfFile = "Homo_sapiens.knownGene.hg19.gtf"
geneList = "hg19_UniProt_IDs.txt"

#get mapping from from UniProt website: http://www.uniprot.org/uploadlists/
#to match annotation script, you'll want to map to 'gene name'

inHandle = open(gtfFile)
line = inHandle.readline()

geneDict = {}

while line:
	line = line.replace("\n","")
	line = line.replace("\r","")
	
	lineInfo = line.split("\t")
	
	result = re.search("; gene_name \"(.+)\";$", lineInfo[8])
	
	if result:
		gene = result.group(1)
		geneDict[gene] = 1
		
	line = inHandle.readline()
	
outHandle = open(geneList, "w")

for gene in geneDict:
	text = gene + "\n"
	outHandle.write(text)
