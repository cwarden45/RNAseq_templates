import sys
import re
import os
import subprocess

parameterFile = "parameters.txt"
finishedSamples = ()

java = "/net/isi-dcnl/ifs/user_data/Seq/jre1.8.0_121/bin/java"
jar_path = "/net/isi-dcnl/ifs/user_data/Seq/"

java_mem = ""
fa_ref = ""
alignmentFolder = ""
strand = ""
email = ""

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

	if param == "Cluster_Email":
		email = value
		
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

if (email == "") or (email == "[required]"):
	print "Need to enter a value for 'Cluster_Email'!"
	sys.exit()
	
fileResults = os.listdir(alignmentFolder)

jobCount = 0
for file in fileResults:
	result = re.search("(.*).bam$",file)
	
	if result:
		sample = result.group(1)
		
		if (sample not in finishedSamples):
			jobCount += 1
			print sample
			shellScript = "separate_GATK_" + sample + ".sh"
			outHandle = open(shellScript, "w")
			text = "#!/bin/bash\n"
			text = text + "#$ -M "+email+"\n"
			text = text + "#$ -m bea\n"
			text = text + "#$ -N sGATK"+str(jobCount)+"\n"
			text = text + "#$ -q short.q\n"
			text = text + "#$ -l vf="+re.sub("g","G",java_mem)+"\n"
			text = text + "#$ -j yes\n"
			text = text + "#$ -o sGATK"+str(jobCount)+".log\n"
			text = text + "#$ -cwd\n"
			text = text + "#$ -V\n"
			outHandle.write(text)
			
			filteredBam = alignmentFolder + "/" + file
			
			outputSubfolder = alignmentFolder +"/" + sample
								
			strandFilter = " -filterName FS -filter \"FS > 30.0\""
			if (strand == "reverse") or (strand == "yes"):
				strandFilter = ""
								
			fullVCF = outputSubfolder + "/" + sample + ".GATK.HC.full.vcf"																			
			text = java + " -Xmx" + java_mem + " -jar "+jar_path+"GenomeAnalysisTK-3.7.jar -T HaplotypeCaller -R " + fa_ref + " -I " + filteredBam + " -o " + fullVCF + " -dontUseSoftClippedBases -stand_call_conf 20.0 --max_alternate_alleles " + max_alt_allele + "\n"
			outHandle.write(text)
			
			#QD = quality score / depth
			#so, QD > 2.0 is not a very strict filter
			bpFilteredVCF = outputSubfolder + "/" + sample + ".GATK.HC.best.practices.flagged.vcf"																			
			text = java + " -Xmx" + java_mem + " -jar "+jar_path+"GenomeAnalysisTK-3.7.jar -T VariantFiltration -R " + fa_ref + " -V " + fullVCF + " -o " + bpFilteredVCF + " -window 35 -cluster 3 -filterName QC -filter \"QD < 2.0\"" + strandFilter + "\n"
			outHandle.write(text)
