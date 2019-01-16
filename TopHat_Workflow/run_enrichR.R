#1st Tested with Data from Liling Hong / Xu Dong from Western University
#I think it is best to first use the Enrichr web-interface, but this may be helfpul in some scenarios

compID = ""
deg.file = ""
up.status = " Up"
down.status = " Down"

output.folder = "../../Results/Round1/enrichR_test"

library(enrichR)

deg.table = read.table(deg.file, head=T, sep="\t")
bg.genes = as.character(deg.table$symbol)
up.genes = bg.genes[deg.table$status == up.status]
down.genes = bg.genes[deg.table$status == down.status]

#dbs = listEnrichrDbs()
dbs = c("KEGG_2016","Reactome_2016","BioCarta_2016",
		"ChEA_2016", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
		"MSigDB_Oncogenic_Signatures",
		"GO_Biological_Process_2018","GO_Cellular_Component_2018","GO_Molecular_Function_2018")

#up-regulated genes
enriched = enrichr(up.genes, dbs)

for (db in dbs){
	print(db)
	result.table = enriched[db]
	
	result.file = paste(output.folder,"/",compID,"_",db,"_UP.txt",sep="")
	write.table(result.table, result.file, quote=F, sep="\t", row.names=F)
}#end for (db in dbs)

#down-regulated genes
enriched = enrichr(down.genes, dbs)

for (db in dbs){
	print(db)
	result.table = enriched[db]
	
	result.file = paste(output.folder,"/",compID,"_",db,"_DOWN.txt",sep="")
	write.table(result.table, result.file, quote=F, sep="\t", row.names=F)
}#end for (db in dbs)
