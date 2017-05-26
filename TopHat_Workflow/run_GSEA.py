import sys
import os
import re

expression_table = "../Result/GSEA/log2_fpkm.gct"
category_table = "../Result/GSEA/log2_fpkm.cls"

outputFolder = "../Result/GSEA/"

comparisons = ["Trt_versus_Cont"]

signature_folder = "/path/to/MSigDB"
signature_files = ["c5.bp.v6.0.symbols.gmt","c5.mf.v6.0.symbols.gmt","c6.all.v6.0.symbols.gmt","h.all.v6.0.symbols.gmt"]
#signature_files = ["c2.all.v6.0.symbols.gmt"]

mac_GSEA_jar = "/path/to/gsea2-2.2.2.jar"
javaMem = "8g"

for comparison in comparisons:
	for signature in signature_files:
		print "Testing " + comparison + " for " + signature
		
		gene_set = signature_folder + "/" + signature
		
		cls_shortID = os.path.basename(category_table)
		#print cls_shortID

		short_gs_ID = re.sub(".symbols.gmt","",signature)
		
		command = "java -Xmx"+javaMem+" -cp "+mac_GSEA_jar + " xtools.gsea.Gsea"
		command = command + " -res " + expression_table
		command = command + " -cls " + category_table + "#" + comparison
		command = command + " -gmx " + gene_set
		command = command + " -rpt_label " + short_gs_ID
		command = command + " -out " + outputFolder + "/" + comparison
		command = command + " -permute gene_set -collapse false -rnd_seed 0"
		
		os.system(command)

