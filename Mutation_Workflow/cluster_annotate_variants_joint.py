import sys
import re
import os

jointVCF = "joint_variant_calls.GATK.HC.best.practices.filtered.sansRNAedit.vcf"
resultFolder = "../../Result/Annotated_Variants/Joint_GATK_Variant_Calls/"

finishedSamples = ()
parameterFile = "parameters.txt"
annovarPath = ""
build = ""
threads = ""
email = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	

	if param == "ANNOVAR_Path":
		annovarPath = value	

	if param == "Threads":
		threads = value
		
	if param == "genome":
		build = value

	if param == "Cluster_Email":
		email = value
		
if (build == "") or (build == "[required]"):
	print "Need to enter a value for 'genome'!"
	sys.exit()
	
if (threads == "") or (threads == "[required]"):
	print "Need to enter a value for 'Threads'!"
	sys.exit()
		
if (annovarPath == "") or (annovarPath == "[required]"):
	print "Need to enter a value for 'ANNOVAR_Path'!"
	sys.exit()

if (email == "") or (email == "[required]"):
	print "Need to enter a value for 'Cluster_Email'!"
	sys.exit()
	
submitAll = "master_ANNOVAR_queue.sh"
masterHandle = open(submitAll,"w")
text = "#!/bin/bash\n"
masterHandle.write(text)	

outPrefix = resultFolder + "joint"
command = annovarPath + "convert2annovar.pl -format vcf4 " + jointVCF + " -allsample -outfile " + outPrefix
os.system(command)

fileResults = os.listdir(resultFolder)

jobCount = 0
for file in fileResults:
	result = re.search("joint.(.*).avinput$",file)
	
	if result:
		sample = result.group(1)
		
		if (sample not in finishedSamples):
			jobCount += 1
			shellScript = "ANNOVAR_" + sample + ".sh"
			text = "qsub " + shellScript + "\n"
			masterHandle.write(text)

			outHandle = open(shellScript, "w")
			text = "#!/bin/bash\n"
			text = text + "#$ -M "+email+"\n"
			text = text + "#$ -m bea\n"
			text = text + "#$ -N ann"+str(jobCount)+"\n"
			text = text + "#$ -q single.q\n"
			text = text + "#$ -l vf=4G\n"
			text = text + "#$ -j yes\n"
			text = text + "#$ -o ann"+str(jobCount)+".log\n"
			text = text + "#$ -cwd\n"
			text = text + "#$ -V\n"
			outHandle.write(text)

			resultSubfolder = resultFolder + "/" + sample
			text = "mkdir " + resultSubfolder + "\n"
			outHandle.write(text)
			
			annovarVar = resultFolder + "joint." + sample + ".avinput"

			annotationPrefix = resultSubfolder + "/" + sample
			if build == "hg19":
				text = annovarPath + "table_annovar.pl --otherinfo " + annovarVar +" " + annovarPath + "humandb/ -csvout -buildver "+build+" -out " + annotationPrefix +" -protocol refGene,clinvar_20160302,cosmic70,nci60,kaviar_20150923,gnomad_genome,hrcr1,dbnsfp30a,avsnp147,cadd13gt20,gwava -operation g,f,f,f,f,f,f,f,f,f,f -nastring NA --thread " + threads + "\n"
				outHandle.write(text)
			else:
				text = annovarPath + "table_annovar.pl --otherinfo " + annovarVar +" " + annovarPath + "humandb/ -csvout -buildver "+build+" -out " + annotationPrefix +" -protocol refGene,clinvar_20160302,cosmic70,nci60,kaviar_20150923,gnomad_genome,hrcr1,dbnsfp30a,avsnp147 -operation g,f,f,f,f,f,f,f,f -nastring NA --thread " + threads + "\n"
				outHandle.write(text)

			bedAnn = build + "_GWAScatalog.bed"
			annotationPrefix = resultSubfolder + "/" + sample + "_annovar_GWAS_Catalog"
			text = annovarPath + "annotate_variation.pl " + annovarVar +" " + annovarPath + "humandb/ -buildver "+build+" -out " + annotationPrefix +" -bedfile " + bedAnn + " -dbtype bed -regionanno -colsWanted 4\n"
			outHandle.write(text)

			bedAnn = build + "_RepeatMasker.bed"
			annotationPrefix = resultSubfolder + "/" + sample + "_annovar_RepeatMasker"
			text = annovarPath + "annotate_variation.pl " + annovarVar +" " + annovarPath + "humandb/ -buildver "+build+" -out " + annotationPrefix +" -bedfile " + bedAnn + " -dbtype bed -regionanno -colsWanted 4\n"
			outHandle.write(text)
			
			#bed annotation mostly works OK, but bedtools was more accurate for ORegAnno hit (probably because it isn't sorted)
			bedAnn = annovarPath + "/humandb/" + build + "_ORegAnno.bed"
			bedOut = resultSubfolder + "/" + sample + "_bedtools_ORegAnno.bed"
			text = "bedtools intersect -wa -wb -a " + annovarVar +" -b " + bedAnn + " > " + bedOut + "\n"
			outHandle.write(text)
			
			#minor formatting modification to work with variant summary script
			#make sure bed_to_ANNOVAR.py is in folder
			oregannoANNOVAR= resultSubfolder + "/" + sample + "_bedtools_ORegAnno.avinput\n"
			text = "python bed_to_ANNOVAR.py --bed="+bedOut+" --annovar="+oregannoANNOVAR+"\n"
			outHandle.write(text)