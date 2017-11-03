import sys
import re
import os
from Bio import SeqIO

folder = "Oases_Assemblies"
num_top_transcripts = 2
inputfile = folder+ "_top_"+str(num_top_transcripts)+"_transcript_IDs_for_BLAST.txt"
outputfile = folder+ "_top_"+str(num_top_transcripts)+"_transcript_seq_for_BLAST.fa"

#folder = "Oases_Assemblies_Trim"
#num_top_transcripts = 2
#inputfile = folder+ "_top_"+str(num_top_transcripts)+"_transcript_IDs_for_BLAST_Trim.txt"
#outputfile = folder+ "_top_"+str(num_top_transcripts)+"_transcript_seq_for_BLAST_Trim.fa"

outHandle = open(outputfile,"w")

sampleHash = {}

inHandle = open(inputfile)
line = inHandle.readline()

lineCount = 0

while line:
	line = re.sub("\r","",line)
	line = re.sub("\n","",line)
	
	lineCount +=1
	
	if lineCount > 1:
		lineInfo = line.split("\t")
		sample = lineInfo[0]
		transcript = lineInfo[1]
		tpm = lineInfo[2]
		
		eXpress_file = folder + "/" + sample + "/transcripts.fa"
		
		if eXpress_file in sampleHash:
			tempHash = sampleHash[eXpress_file]
			tempHash[transcript] = tpm
			sampleHash[eXpress_file] = tempHash
		else:
			tempHash = {}
			tempHash[transcript] = tpm
			sampleHash[eXpress_file] = tempHash		
	
	line = inHandle.readline()
	
for inputfile in sampleHash:
	sample = inputfile
	sample = re.sub(folder + "/","",sample)
	sample = re.sub("/transcripts.fa","",sample)
	print sample
	tempHash = sampleHash[inputfile]
	
	fasta_sequences = SeqIO.parse(inputfile,'fasta')
	for fasta in fasta_sequences:
		transName = fasta.id
		transSequence = str(fasta.seq)
		
		if transName in tempHash:
			text = ">" + sample + "_" + transName + "\n"
			text = text + transSequence + "\n"
			outHandle.write(text)