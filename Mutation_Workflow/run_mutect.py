import sys
import re
import os

java_mem = "10g"

normal = "/path/to/normal.bam"
tumor = "/path/to/tumor.bam"
vcf = "sample.mutect2.vcf"

fa_ref = "/path/to/hg19.fa"

command = "java -Xmx" + java_mem + " -jar /opt/GenomeAnalysisTK-3.6.jar -T MuTect2 -R " + fa_ref + " -I:normal " + normal + " -I:tumor " + tumor + " -o " + vcf + " -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0  "
os.system(command)
