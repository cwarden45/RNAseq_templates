import sys
import re
import os

alignmentFolder = ""
genome = ""
igvtools = "/opt/igvtools_2.3.91/igvtools"

fileResults = os.listdir(alignmentFolder)

for file in fileResults:
	bamResult = re.search("(.*).bam$",file)
	if bamResult:
		print file
		fullPath = alignmentFolder + "/" + file

		tdf = fullPath + ".tdf"
		command = igvtools + " count -z 7 " + fullPath + " " + tdf + " " + genome
		os.system(command)