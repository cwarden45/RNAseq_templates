import sys
import re
import os
import subprocess

parameterFile = "parameters.txt"
relatedSamples = ()

java_mem = ""
fa_ref = ""
alignmentFolder = ""
strand = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]

	if param == "strand":
		strand = value
	
	if param == "Alignment_Folder":
		alignmentFolder = value
		
	if param == "FA_Ref":
		fa_ref = value
		
	if param == "Java_Mem":
		java_mem = value
		
if (strand == "") or (strand == "[required]"):
	print "Need to enter a value for 'strand'!"
	sys.exit()
		
if (java_mem == "") or (java_mem == "[required]"):
	print "Need to enter a value for 'Java_Mem'!"
	sys.exit()

if (fa_ref == "") or (fa_ref == "[required]"):
	print "Need to enter a value for 'FA_Ref'!"
	sys.exit()
	
	
if (alignmentFolder == "") or (alignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder'!"
	sys.exit()
	
fileResults = os.listdir(alignmentFolder)

combinedFiles = ""

for file in fileResults:
	result = re.search("(.*).bam$",file)
	
	if result:
		sample = result.group(1)
		
		if (sample in relatedSamples):
			print sample
			filteredBam = alignmentFolder + "/" + file
			
			outputSubfolder = alignmentFolder +"/" + sample

			fullVCF = outputSubfolder + "/" + sample + ".GATK.HC.full.AS.g.vcf"
			combinedFiles = combinedFiles + " --variant " + fullVCF
			command = "java -Xmx" + java_mem + " -jar /opt/GenomeAnalysisTK-3.6.jar -T HaplotypeCaller -R " + fa_ref + " -I " + filteredBam + " -o " + fullVCF + " -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 --emitRefConfidence GVCF -G Standard -G AS_Standard"
			os.system(command)			

combinedVCF = "joint_variant_calls.GATK.HC.best.practices.all.vcf"																			
command = "java -Xmx" + java_mem + " -jar /opt/GenomeAnalysisTK-3.6.jar -T GenotypeGVCFs -R " + fa_ref + combinedFiles + " -o " + combinedVCF
os.system(command)
			
strandFilter = " -filterName FS -filter \"FS > 30.0\""
if (strand == "reverse") or (strand == "yes"):
	strandFilter = ""
			
#QD = quality score / depth
#so, QD > 2.0 is not a very strict filter
bpFilteredVCF = "joint_variant_calls.GATK.HC.best.practices.flagged.vcf"																			
command = "java -Xmx" + java_mem + " -jar /opt/GenomeAnalysisTK-3.6.jar -T VariantFiltration -R " + fa_ref + " -V " + combinedVCF + " -o " + bpFilteredVCF + " -window 35 -cluster 3 -filterName QC -filter \"QD < 2.0\"" + strandFilter
os.system(command)