import sys
import re
import os
import subprocess

parameterFile = "parameters.txt"
finishedSamples = ()

java_mem = ""
fa_ref = ""
alignmentFolder = ""
strand = ""
max_alt_allele = ""

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

	if param == "Max_Alt_Alleles":
		max_alt_allele = value
		
if (strand == "") or (strand == "[required]"):
	print "Need to enter a value for 'strand'!"
	sys.exit()
		
if (java_mem == "") or (java_mem == "[required]"):
	print "Need to enter a value for 'Java_Mem'!"
	sys.exit()

if (fa_ref == "") or (fa_ref == "[required]"):
	print "Need to enter a value for 'FA_Ref'!"
	sys.exit()

if (max_alt_allele == "") or (max_alt_allele == "[required]"):
	print "Need to enter a value for 'Max_Alt_Alleles'!"
	sys.exit()	
	
if (alignmentFolder == "") or (alignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder'!"
	sys.exit()
	
fileResults = os.listdir(alignmentFolder)

for file in fileResults:
	result = re.search("(.*).bam$",file)
	
	if result:
		sample = result.group(1)
		
		if (sample not in finishedSamples):
			print sample
			filteredBam = alignmentFolder + "/" + file
			
			outputSubfolder = alignmentFolder +"/" + sample
								
			strandFilter = " -filterName FS -filter \"FS > 30.0\""
			if (strand == "reverse") or (strand == "yes"):
				strandFilter = ""
								
			fullVCF = outputSubfolder + "/" + sample + ".GATK.HC.full.vcf"																			
			#command = "java -Xmx" + java_mem + " -jar /opt/GenomeAnalysisTK-3.6.jar -T HaplotypeCaller -R " + fa_ref + " -I " + filteredBam + " -o " + fullVCF + " -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 --max_alternate_alleles " + max_alt_allele
			command = "java -Xmx" + java_mem + " -jar /opt/GATK-3.7/GenomeAnalysisTK.jar -T HaplotypeCaller -R " + fa_ref + " -I " + filteredBam + " -o " + fullVCF + " -dontUseSoftClippedBases -stand_call_conf 20.0 --max_alternate_alleles " + max_alt_allele
			os.system(command)
			
			#QD = quality score / depth
			#so, QD > 2.0 is not a very strict filter
			bpFilteredVCF = outputSubfolder + "/" + sample + ".GATK.HC.best.practices.flagged.vcf"																			
			#command = "java -Xmx" + java_mem + " -jar /opt/GenomeAnalysisTK-3.6.jar -T VariantFiltration -R " + fa_ref + " -V " + fullVCF + " -o " + bpFilteredVCF + " -window 35 -cluster 3 -filterName QC -filter \"QD < 2.0\"" + strandFilter
			command = "java -Xmx" + java_mem + " -jar /opt/GATK-3.7/GenomeAnalysisTK.jar -T VariantFiltration -R " + fa_ref + " -V " + fullVCF + " -o " + bpFilteredVCF + " -window 35 -cluster 3 -filterName QC -filter \"QD < 2.0\"" + strandFilter
			os.system(command)
