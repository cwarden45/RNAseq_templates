import sys
import re
import os
import subprocess

#input and output is for hard filtering applied in this code (RepeatMasker and RNA-editing filters derived from `vcfExtOut`)
vcfExtIn = ".GATK.HC.best.practices.flagged.vcf"
vcfExtOut = ".GATK.HC.best.practices.filtered.vcf"
summary_file = "GATK_best_practices_variant_counts.txt"

parameterFile = "parameters.txt"
finishedSamples = ()

alignmentFolder = ""
repeatBED = ""
rnaEditBED = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	
	if param == "Alignment_Folder":
		alignmentFolder = value

	if param == "RepeatMasker_BED":
		repeatBED = value

	if param == "RNA_Editing_BED":
		rnaEditBED = value		
	
if (alignmentFolder == "") or (alignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder'!"
	sys.exit()
	
if (repeatBED == "") or (repeatBED == "[required]"):
	print "Need to enter a value for 'RepeatMasker_BED'!"
	sys.exit()
	
if (rnaEditBED == "") or (rnaEditBED == "[required]"):
	print "Need to enter a value for 'RNA_Editing_BED'!"
	sys.exit()
	
fileResults = os.listdir(alignmentFolder)

statHandle = inHandle = open(summary_file,"w")
text = "Sample\tSNP.count\tIns.count\tDel.count\tnoRepeat.SNP.count\tnoRepeat.Ins.count\tnoRepeat.Del.count\tnoRNAedit.SNP.count\tnoRNAedit.Ins.count\tnoRNAedit.Del.count\n"
statHandle.write(text)

for file in fileResults:
	result = re.search("(.*).bam$",file)
	
	if result:
		sample = result.group(1)
		
		if (sample not in finishedSamples):
			print sample
			filteredBam = alignmentFolder + "/" + file
			
			outputSubfolder = alignmentFolder +"/" + sample
			
			vcfFull = outputSubfolder + "/" + sample + vcfExtIn																			
			vcfFlagFilter = outputSubfolder + "/" + sample + vcfExtOut
			
			snpCounts = 0
			insertionCounts = 0
			deletionCounts = 0
			
			outHandle = open(vcfFlagFilter,"w")
			
			inHandle = open(vcfFull)
			line = inHandle.readline()
			
			while line:
				
				commentResult = re.search("^#",line)
				
				if commentResult:
					outHandle.write(line)
				else:
					lineInfo = line.split("\t")
					chr = lineInfo[0]
					variantStatus = lineInfo[6]
					
					if (variantStatus == "PASS"):

						refSeq = lineInfo[3]
						varSeq = lineInfo[4]
						
						if (len(refSeq) == 1) and (len(varSeq) == 1):
							snpCounts += 1
						else:
							multVarResult = re.search(",",varSeq)
							
							if multVarResult:
								refVars = varSeq.split(",")
								for testVar in refVars:
									if testVar != refSeq:
										if (len(refSeq) == 1) and (len(testVar) == 1):
											snpCounts += 1
										elif len(testVar) > len(refSeq):
											insertionCounts += 1
										elif len(testVar) < len(refSeq):
											deletionCounts += 1
										else:
											print "Modify code to count ref: " + refSeq + ", var: " + varSeq										
							elif len(varSeq) > len(refSeq):
								insertionCounts += 1
							elif len(varSeq) < len(refSeq):
								deletionCounts += 1
							else:
								print "Modify code to count ref: " + refSeq + ", var: " + varSeq
							
						outHandle.write(line)	
					elif (variantStatus != "SnpCluster") and (variantStatus != "QC")and (variantStatus != "QC;SnpCluster"):
						print "Modify code to Keep or Skip variant status: " + variantStatus
						sys.exit()
				
				line = inHandle.readline()
			inHandle.close()
			outHandle.close()
			
			text = sample + "\t" + str(snpCounts) + "\t" + str(insertionCounts) + "\t" + str(deletionCounts)

			#RepeatMasker filer
			noRepeatSnps = 0
			noRepeatIns = 0
			noRepeatDel = 0
			
			noRepeatVCF = re.sub(".vcf$",".sansRepeatMasker.vcf",vcfFlagFilter)
			
			command = "/opt/bedtools2/bin/bedtools intersect -header -v -a " + vcfFlagFilter + " -b " + repeatBED + " > " + noRepeatVCF
			os.system(command)

			inHandle = open(noRepeatVCF)
			line = inHandle.readline()
			
			while line:
				
				commentResult = re.search("^#",line)
				
				if not commentResult:
					lineInfo = line.split("\t")
					chr = lineInfo[0]
					variantStatus = lineInfo[6]
					
					if (variantStatus == "PASS"):

						refSeq = lineInfo[3]
						varSeq = lineInfo[4]
						
						if (len(refSeq) == 1) and (len(varSeq) == 1):
							noRepeatSnps += 1
						else:
							multVarResult = re.search(",",varSeq)
							
							if multVarResult:
								refVars = varSeq.split(",")
								for testVar in refVars:
									if testVar != refSeq:
										if (len(refSeq) == 1) and (len(testVar) == 1):
											noRepeatSnps += 1
										elif len(testVar) > len(refSeq):
											noRepeatIns += 1
										elif len(testVar) < len(refSeq):
											noRepeatDel += 1
										else:
											print "Modify code to count ref: " + refSeq + ", var: " + varSeq										
							elif len(varSeq) > len(refSeq):
								noRepeatIns += 1
							elif len(varSeq) < len(refSeq):
								noRepeatDel += 1
							else:
								print "Modify code to count ref: " + refSeq + ", var: " + varSeq	
					elif (variantStatus != "SnpCluster") and (variantStatus != "QC")and (variantStatus != "QC;SnpCluster"):
						print "Modify code to Keep or Skip variant status: " + variantStatus
						sys.exit()
				
				line = inHandle.readline()
			inHandle.close()
			
			text = text + "\t" + str(noRepeatSnps) + "\t" + str(noRepeatIns) + "\t" + str(noRepeatDel)
			
			#RNA-editing filter
			noRNAeditSnps = 0
			noRNAeditIns = 0
			noRNAeditDel = 0
			
			noRNAeditVCF = re.sub(".vcf$",".sansRNAedit.vcf",vcfFlagFilter)

			command = "/opt/bedtools2/bin/bedtools intersect -header -v -a " + noRepeatVCF + " -b " + rnaEditBED + " > " + noRNAeditVCF
			os.system(command)

			inHandle = open(noRNAeditVCF)
			line = inHandle.readline()
			
			while line:
				
				commentResult = re.search("^#",line)
				
				if not commentResult:
					lineInfo = line.split("\t")
					chr = lineInfo[0]
					variantStatus = lineInfo[6]
					
					if (variantStatus == "PASS"):

						refSeq = lineInfo[3]
						varSeq = lineInfo[4]
						
						if (len(refSeq) == 1) and (len(varSeq) == 1):
							noRNAeditSnps += 1
						else:
							multVarResult = re.search(",",varSeq)
							
							if multVarResult:
								refVars = varSeq.split(",")
								for testVar in refVars:
									if testVar != refSeq:
										if (len(refSeq) == 1) and (len(testVar) == 1):
											noRNAeditSnps += 1
										elif len(testVar) > len(refSeq):
											noRNAeditIns += 1
										elif len(testVar) < len(refSeq):
											noRNAeditDel += 1
										else:
											print "Modify code to count ref: " + refSeq + ", var: " + varSeq										
							elif len(varSeq) > len(refSeq):
								noRNAeditIns += 1
							elif len(varSeq) < len(refSeq):
								noRNAeditDel += 1
							else:
								print "Modify code to count ref: " + refSeq + ", var: " + varSeq	
					elif (variantStatus != "SnpCluster") and (variantStatus != "QC")and (variantStatus != "QC;SnpCluster"):
						print "Modify code to Keep or Skip variant status: " + variantStatus
						sys.exit()
				
				line = inHandle.readline()
			inHandle.close()
			
			text = text + "\t" + str(noRNAeditSnps) + "\t" + str(noRNAeditIns) + "\t" + str(noRNAeditDel) + "\n"
			statHandle.write(text)
