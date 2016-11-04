import sys
import re
import os

normal = "path/to/tumor.pileup"
tumor = "path/to/normal.pileup"
varscanPrefix = ""

varscan_jar = "/path/to/VarScan.v2.4.2.jar"
java_mem = "16g"
email = ""
strand = ""

strandFilter = "1"
if strand != "no":
	strandFilter = "0"

shellScript = varscanPrefix + ".sh"
outHandle = open(shellScript, "w")
text = "#!/bin/bash\n"
text = text + "#$ -M "+email+"\n"
text = text + "#$ -m bea\n"
text = text + "#$ -N VSsomatic\n"
text = text + "#$ -q all.q\n"
text = text + "#$ -l vf=16G\n"
text = text + "#$ -j yes\n"
text = text + "#$ -o VSsomatic.log\n"
text = text + "#$ -cwd\n"
text = text + "#$ -V\n"
outHandle.write(text)
	
#you can create .vcf with --output-vcf, but then you have to parse out somatic variants and combine SNP and indel files
text = "java -Xmx" + java_mem+ " -jar "+varscan_jar+" somatic " + normal + " " + tumor + " " + varscanPrefix + " --min-var-freq 0.3 --min-avg-qual 20 --p-value 0.01 --somatic-p-value 0.01 --min-coverage-normal 10 --min-coverage-tumor 10 --strand-filter " + strandFilter
outHandle.write(text)
