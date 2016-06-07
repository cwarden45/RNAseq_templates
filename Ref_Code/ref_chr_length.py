import sys
import re
import os

def chrLength(inputfile, outputfile):
	chrHash = {}

	name = ""
	length = 0
	
	inHandle = open(inputfile)
	line = inHandle.readline()

	while line:	
		line = re.sub("\n","",line)
		line = re.sub("\r","",line)
		
		result = re.search("^>",line)
		
		if result:
			if name != "":
				chrHash[name] = length
			name = re.sub("^>","",line)
			length = 0
		else:
			length += len(line)
		line = inHandle.readline()
		
	chrHash[name] = length

	outHandle = open(outputfile, "w")
	text = "Chr\tLength\n"
	outHandle.write(text)
	
	for chr in chrHash:
		length = chrHash[chr]
		text = chr + "\t" + str(length) + "\n"
		outHandle.write(text)

faRef = "/path/to/ref.fa"
lengthFile = "ref_chr_length.txt"
chrLength(faRef, lengthFile)
		
