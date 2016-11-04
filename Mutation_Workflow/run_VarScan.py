import sys
import re
import os

normal = "path/to/tumor.pileup"
tumor = "path/to/normal.pileup"
varscanPrefix = ""
strand = ""
java_mem = "10g"

strandFilter = "1"
if strand != "no":
	strandFilter = "0"

#you can create .vcf with --output-vcf, but then you have to parse out somatic variants and combine SNP and indel files
command = "java -Xmx" + java_mem+ " -jar /opt/varscan/VarScan.v2.4.2.jar somatic " + normal + " " + tumor + " " + varscanPrefix + " --min-var-freq 0.3 --min-avg-qual 20 --p-value 0.01 --somatic-p-value 0.01 --min-coverage-normal 10 --min-coverage-tumor 10 --strand-filter " + strandFilter
os.system(command)
