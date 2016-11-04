import sys
import re
import os

varscanPrefixes = ()
resultFolder = "Annotated_Variants/VarScan2"
countFile = "varscan_somatic_variant_count.txt"

parameterFile = "parameters.txt"
finishedSamples = ()

annovarPath = ""
build = ""
threads = ""

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

if (build == "") or (build == "[required]"):
	print "Need to enter a value for 'genome'!"
	sys.exit()
	
if (threads == "") or (threads == "[required]"):
	print "Need to enter a value for 'Threads'!"
	sys.exit()
	
if (annovarPath == "") or (annovarPath == "[required]"):
	print "Need to enter a value for 'ANNOVAR_Path'!"
	sys.exit()

statHandle = open(countFile,"w")
text = "pairID\tSNPcount\tDel.Count\tIns.Count\n"
statHandle.write(text)
	
for varscanPrefix in varscanPrefixes:
	print varscanPrefix
	
	resultSubfolder = resultFolder + "/" + varscanPrefix
	command = "mkdir " + resultSubfolder
	os.system(command)
	
	varScanSnp = varscanPrefix + ".snp"
	varScanIndel = varscanPrefix + ".indel"
	snpCount = 0
	insCount = 0
	delCount = 0
			
	annovarVar = resultSubfolder + "/" + varscanPrefix + ".avinput"
	outHandle = open(annovarVar,"w")
	#chr	start	stop	ref	var	status	qual(p-value)	cov

	#convert SNP file
	inHandle = open(varScanSnp)
	line = inHandle.readline()
	
	lineCount = 0
	
	while line:
		lineCount += 1
		
		if lineCount > 1:
			lineInfo = line.split("\t")
			chr = lineInfo[0]
			pos = lineInfo[1]
			ref = lineInfo[2]
			var = lineInfo[3]
			status = lineInfo[12]
			somatic_pvalue = lineInfo[14]
			tumorRefCov = int(lineInfo[8])
			tumorVarCov = int(lineInfo[9])
			totalCov = tumorRefCov + tumorVarCov
			
			if (len(ref) == 1 or len(var) == 1) and (status != "Germline"):
				if status == "Somatic":
					snpCount += 1
				#skip variants with multiple variant alleles --> separate frequencies aren't provided anyways
				text = chr + "\t" + pos+ "\t" + pos + "\t" + ref + "\t" + var + "\t" + status + "\tp=" + somatic_pvalue + "\t" + str(totalCov) + "\n"
				outHandle.write(text)
				
		line = inHandle.readline()
	inHandle.close()

	#convert indel file
	inHandle = open(varScanIndel)
	line = inHandle.readline()
	
	lineCount = 0
	
	while line:
		lineCount += 1
		
		if lineCount > 1:
			lineInfo = line.split("\t")
			chr = lineInfo[0]
			pos = int(lineInfo[1])
			ref = lineInfo[2]
			var = lineInfo[3]
			status = lineInfo[12]
			somatic_pvalue = lineInfo[14]
			tumorRefCov = int(lineInfo[8])
			tumorVarCov = int(lineInfo[9])
			totalCov = tumorRefCov + tumorVarCov
			
			insFlag = re.search("^\+([ACGT]+)",var)
			delFlag = re.search("^-([ACGT]+)",var)
			
			if insFlag and (status != "Germline"):
				if status == "Somatic":
					insCount += 1
					
				var = insFlag.group(1)
				text = chr + "\t" + str(pos)+ "\t" + str(pos+len(var)-1) + "\t-\t" + var + "\t" + status + "\tp=" + somatic_pvalue + "\t" + str(totalCov) + "\n"
				outHandle.write(text)

			if delFlag and (status != "Germline"):
				if status == "Somatic":
					delCount += 1
				var = delFlag.group(1)
				text = chr + "\t" + str(pos+1)+ "\t" + str(pos+1+len(var)-1) + "\t" + var + "\t-\t" + status + "\tp=" + somatic_pvalue + "\t" + str(totalCov) + "\n"
				outHandle.write(text)
				
		line = inHandle.readline()
	inHandle.close()
	
	outHandle.close()
	
	text = varscanPrefix + "\t"+str(snpCount)+"\t"+str(delCount)+"\t"+str(insCount)+"\n"
	statHandle.write(text)

	annotationPrefix = resultSubfolder + "/" + varscanPrefix
	if build == "hg19":
		command = annovarPath + "table_annovar.pl --otherinfo " + annovarVar +" " + annovarPath + "humandb/ -csvout -buildver "+build+" -out " + annotationPrefix +" -protocol refGene,clinvar_20160302,cosmic70,nci60,kaviar_20150923,hrcr1,dbnsfp30a,exac03,avsnp147,cadd13gt20,gwava -operation g,f,f,f,f,f,f,f,f,f,f -nastring NA --thread " + threads
		os.system(command)
	else:
		command = annovarPath + "table_annovar.pl --otherinfo " + annovarVar +" " + annovarPath + "humandb/ -csvout -buildver "+build+" -out " + annotationPrefix +" -protocol refGene,clinvar_20160302,cosmic70,nci60,kaviar_20150923,hrcr1,dbnsfp30a,exac03,avsnp147 -operation g,f,f,f,f,f,f,f,f -nastring NA --thread " + threads
		os.system(command)

	bedAnn = build + "_GWAScatalog.bed"
	annotationPrefix = resultSubfolder + "/" + varscanPrefix + "_annovar_GWAS_Catalog"
	command = annovarPath + "annotate_variation.pl " + annovarVar +" " + annovarPath + "humandb/ -buildver "+build+" -out " + annotationPrefix +" -bedfile " + bedAnn + " -dbtype bed -regionanno -colsWanted 4"
	os.system(command)

	#bed annotation mostly works OK, but bedtools was more accurate for ORegAnno hit (probably because it isn't sorted)
	bedAnn = annovarPath + "/humandb/" + build + "_ORegAnno.bed"
	bedOut = resultSubfolder + "/" + varscanPrefix + "_bedtools_ORegAnno.bed"
	command = "/opt/bedtools2/bin/bedtools intersect -wa -wb -a " + annovarVar +" -b " + bedAnn + " > " + bedOut
	os.system(command)
	
	#minor formatting modification to work with variant summary script
	oregannoANNOVAR= resultSubfolder + "/" + varscanPrefix + "_bedtools_ORegAnno.avinput"
	outHandle = open(oregannoANNOVAR,"w")
	
	inHandle = open(bedOut)
	line = inHandle.readline()
	
	while line:
		lineInfo = line.split("\t")
		chr = lineInfo[0]
		start = lineInfo[1]
		stop = lineInfo[2]
		ref = lineInfo[3]
		var = lineInfo[4]
		type = lineInfo[5]
		somatic_pvalue = lineInfo[6]
		totalCov = lineInfo[7]
		oreg = lineInfo[11]
			
		text = chr + "\t" + start+ "\t" + stop + "\t" + ref + "\t" + var + "\t" + type + "\t" + somatic_pvalue+ "\t" + oreg + "\n"
		outHandle.write(text)
				
		line = inHandle.readline()
	inHandle.close()
	outHandle.close()	