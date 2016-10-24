import sys
import re
import os
import subprocess

#input and output is for hard filtering applied in this code (RepeatMasker and RNA-editing filters derived from `vcfOut`)
vcfIn = "joint_variant_calls.GATK.HC.best.practices.flagged.vcf"
vcfOut = "joint_variant_calls.GATK.HC.best.practices.filtered.vcf"
summary_file = "GATK_joint_variant_counts.txt"

parameterFile = "parameters.txt"

repeatBED = ""
rnaEditBED = ""

sampleTextHash = {}

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]

	if param == "RepeatMasker_BED":
		repeatBED = value

	if param == "RNA_Editing_BED":
		rnaEditBED = value		
	
if (repeatBED == "") or (repeatBED == "[required]"):
	print "Need to enter a value for 'RepeatMasker_BED'!"
	sys.exit()
	
if (rnaEditBED == "") or (rnaEditBED == "[required]"):
	print "Need to enter a value for 'RNA_Editing_BED'!"
	sys.exit()

statHandle = inHandle = open(summary_file,"w")
text = "Sample\tSNP.count\tIns.count\tDel.count\tnoRepeat.SNP.count\tnoRepeat.Ins.count\tnoRepeat.Del.count\tnoRNAedit.SNP.count\tnoRNAedit.Ins.count\tnoRNAedit.Del.count\n"
statHandle.write(text)

sample = "Joint"

snpCounts = 0
insertionCounts = 0
deletionCounts = 0
			
outHandle = open(vcfOut,"w")
			
inHandle = open(vcfIn)
line = inHandle.readline()
	
sampleIndexHash ={}
	
while line:
				
	commentResult = re.search("^#",line)
				
	if commentResult:
		outHandle.write(line)
		
		sampleIDresult = re.search("^#CHR",line)
		if sampleIDresult:
			lineInfo = line.split("\t")
			for i in xrange(9,len(lineInfo)):
				lineInfo[i] = re.sub("\r","",lineInfo[i])
				lineInfo[i] = re.sub("\n","",lineInfo[i])
				sampleIndexHash[lineInfo[i]]=i
	else:
		lineInfo = line.split("\t")
		chr = lineInfo[0]
		variantStatus = lineInfo[6]
					
		if (variantStatus == "PASS"):
			refSeq = lineInfo[3]
			varSeq = lineInfo[4]
						
			if (len(refSeq) == 1) and (len(varSeq) == 1):
				snpCounts += 1
				
				for sampleID in sampleIndexHash:
					lineInfo[sampleIndexHash[sampleID]] = re.sub("\r","",lineInfo[sampleIndexHash[sampleID]])
					lineInfo[sampleIndexHash[sampleID]] = re.sub("\n","",lineInfo[sampleIndexHash[sampleID]])
					genotype = lineInfo[sampleIndexHash[sampleID]][0:3]
					if genotype != "0/0":
						if sampleID in sampleTextHash:
							countInfo = sampleTextHash[sampleID].split("\t")
							countInfo[0] = str(int(countInfo[0]) + 1)
							sampleTextHash[sampleID] = "\t".join(countInfo)
						else:
							sampleTextHash[sampleID] = "1\t0\t0\t0\t0\t0\t0\t0\t0"
			else:
				multVarResult = re.search(",",varSeq)
				
				print "Need to add code for providing joint mixed variant counts"
			
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
								#print line
								print "Modify code to count ref: " + refSeq + ", var: " + varSeq										
				elif len(varSeq) > len(refSeq):
					insertionCounts += 1
					for sampleID in sampleIndexHash:
						lineInfo[sampleIndexHash[sampleID]] = re.sub("\r","",lineInfo[sampleIndexHash[sampleID]])
						lineInfo[sampleIndexHash[sampleID]] = re.sub("\n","",lineInfo[sampleIndexHash[sampleID]])
						genotype = lineInfo[sampleIndexHash[sampleID]][0:3]
						if genotype != "0/0":
							if sampleID in sampleTextHash:
								countInfo = sampleTextHash[sampleID].split("\t")
								countInfo[1] = str(int(countInfo[1]) + 1)
								sampleTextHash[sampleID] = "\t".join(countInfo)
							else:
								sampleTextHash[sampleID] = "0\t1\t0\t0\t0\t0\t0\t0\t0"
				elif len(varSeq) < len(refSeq):
					deletionCounts += 1
					for sampleID in sampleIndexHash:
						lineInfo[sampleIndexHash[sampleID]] = re.sub("\r","",lineInfo[sampleIndexHash[sampleID]])
						lineInfo[sampleIndexHash[sampleID]] = re.sub("\n","",lineInfo[sampleIndexHash[sampleID]])
						genotype = lineInfo[sampleIndexHash[sampleID]][0:3]
						if genotype != "0/0":
							if sampleID in sampleTextHash:
								countInfo = sampleTextHash[sampleID].split("\t")
								countInfo[2] = str(int(countInfo[2]) + 1)
								sampleTextHash[sampleID] = "\t".join(countInfo)
							else:
								sampleTextHash[sampleID] = "0\t0\t1\t0\t0\t0\t0\t0\t0"
				else:
					#print line
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
			
noRepeatVCF = re.sub(".vcf$",".sansRepeatMasker.vcf",vcfOut)
			
command = "/opt/bedtools2/bin/bedtools intersect -header -v -a " + vcfOut + " -b " + repeatBED + " > " + noRepeatVCF
#os.system(command)

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
				for sampleID in sampleIndexHash:
					lineInfo[sampleIndexHash[sampleID]] = re.sub("\r","",lineInfo[sampleIndexHash[sampleID]])
					lineInfo[sampleIndexHash[sampleID]] = re.sub("\n","",lineInfo[sampleIndexHash[sampleID]])
					genotype = lineInfo[sampleIndexHash[sampleID]][0:3]
					if genotype != "0/0":
						countInfo = sampleTextHash[sampleID].split("\t")
						countInfo[3] = str(int(countInfo[3]) + 1)
						sampleTextHash[sampleID] = "\t".join(countInfo)
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
					for sampleID in sampleIndexHash:
							lineInfo[sampleIndexHash[sampleID]] = re.sub("\r","",lineInfo[sampleIndexHash[sampleID]])
							lineInfo[sampleIndexHash[sampleID]] = re.sub("\n","",lineInfo[sampleIndexHash[sampleID]])
							genotype = lineInfo[sampleIndexHash[sampleID]][0:3]
							if genotype != "0/0":
								countInfo = sampleTextHash[sampleID].split("\t")
								countInfo[4] = str(int(countInfo[4]) + 1)
								sampleTextHash[sampleID] = "\t".join(countInfo)
				elif len(varSeq) < len(refSeq):
					noRepeatDel += 1
					for sampleID in sampleIndexHash:
							lineInfo[sampleIndexHash[sampleID]] = re.sub("\r","",lineInfo[sampleIndexHash[sampleID]])
							lineInfo[sampleIndexHash[sampleID]] = re.sub("\n","",lineInfo[sampleIndexHash[sampleID]])
							genotype = lineInfo[sampleIndexHash[sampleID]][0:3]
							if genotype != "0/0":
								countInfo = sampleTextHash[sampleID].split("\t")
								countInfo[5] = str(int(countInfo[5]) + 1)
								sampleTextHash[sampleID] = "\t".join(countInfo)
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
			
noRNAeditVCF = re.sub(".vcf$",".sansRNAedit.vcf",vcfOut)

command = "/opt/bedtools2/bin/bedtools intersect -header -v -a " + noRepeatVCF + " -b " + rnaEditBED + " > " + noRNAeditVCF
#os.system(command)

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
				for sampleID in sampleIndexHash:
					lineInfo[sampleIndexHash[sampleID]] = re.sub("\r","",lineInfo[sampleIndexHash[sampleID]])
					lineInfo[sampleIndexHash[sampleID]] = re.sub("\n","",lineInfo[sampleIndexHash[sampleID]])
					genotype = lineInfo[sampleIndexHash[sampleID]][0:3]
					if genotype != "0/0":
						countInfo = sampleTextHash[sampleID].split("\t")
						countInfo[6] = str(int(countInfo[6]) + 1)
						sampleTextHash[sampleID] = "\t".join(countInfo)
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
					for sampleID in sampleIndexHash:
						lineInfo[sampleIndexHash[sampleID]] = re.sub("\r","",lineInfo[sampleIndexHash[sampleID]])
						lineInfo[sampleIndexHash[sampleID]] = re.sub("\n","",lineInfo[sampleIndexHash[sampleID]])
						genotype = lineInfo[sampleIndexHash[sampleID]][0:3]
						if genotype != "0/0":
							countInfo = sampleTextHash[sampleID].split("\t")
							countInfo[7] = str(int(countInfo[7]) + 1)
							sampleTextHash[sampleID] = "\t".join(countInfo)
				elif len(varSeq) < len(refSeq):
					noRNAeditDel += 1
					for sampleID in sampleIndexHash:
						lineInfo[sampleIndexHash[sampleID]] = re.sub("\r","",lineInfo[sampleIndexHash[sampleID]])
						lineInfo[sampleIndexHash[sampleID]] = re.sub("\n","",lineInfo[sampleIndexHash[sampleID]])
						genotype = lineInfo[sampleIndexHash[sampleID]][0:3]
						if genotype != "0/0":
							countInfo = sampleTextHash[sampleID].split("\t")
							countInfo[8] = str(int(countInfo[8]) + 1)
							sampleTextHash[sampleID] = "\t".join(countInfo)
				else:
					print "Modify code to count ref: " + refSeq + ", var: " + varSeq	
		elif (variantStatus != "SnpCluster") and (variantStatus != "QC")and (variantStatus != "QC;SnpCluster"):
			print "Modify code to Keep or Skip variant status: " + variantStatus
			sys.exit()
				
	line = inHandle.readline()
inHandle.close()
			
text = text + "\t" + str(noRNAeditSnps) + "\t" + str(noRNAeditIns) + "\t" + str(noRNAeditDel) + "\n"
statHandle.write(text)

for sampleID in sampleTextHash:
	text =sampleID + "\t" + sampleTextHash[sampleID] + "\n"
	statHandle.write(text)