compID = ""
degPrefix = ""
outputFolder = "../Results/Round1/GO/web_input_files"
up.status = " Up"
down.status = " Down"
genome = ""

up.file = paste(outputFolder,"/",compID,"_UP_symbol.txt",sep="")
down.file = paste(outputFolder,"/",compID,"_DOWN_symbol.txt",sep="")
bg.file = paste(outputFolder,"/",compID,"_background_symbol.txt",sep="")
deg.table = read.table(paste(degPrefix,".txt",sep=""),head=T, sep="\t")

bg.genes = as.character(deg.table$symbol)
write.table(data.frame(bg.genes),bg.file, col.names=F, row.names=F, sep="\t", quote=F)

up.genes = bg.genes[deg.table$status == up.status]
write.table(data.frame(up.genes),up.file, col.names=F, row.names=F, sep="\t", quote=F)

down.genes = bg.genes[deg.table$status == down.status]
write.table(data.frame(down.genes),down.file, col.names=F, row.names=F, sep="\t", quote=F)

#DAVID requires geneID (such as ENTREZ_GENE_ID), if defining background list
#...My preference would be enrichr or GATHER (which use gene symbols, but without custom set of background genes)
up.file = paste(outputFolder,"/",compID,"_UP_geneID.txt",sep="")
down.file = paste(outputFolder,"/",compID,"_DOWN_geneID.txt",sep="")
bg.file = paste(outputFolder,"/",compID,"_background_geneID.txt",sep="")

if((genome == "hg19")|(genome == "hg38")){
	library(org.Hs.eg.db)
	orgdb = org.Hs.eg.db
}else if((genome == "mm9")|(genome == "mm10")){
	library(org.Mm.eg.db)
	orgdb = org.Mm.eg.db
}else{
	stop(paste("Define gene mapping for genome: ",genome,ep=""))
}

gene.symbols = keys(orgdb, keytype="SYMBOL")
genecols = c("SYMBOL", "ENTREZID","GENENAME")
genetable = select(orgdb, keys=gene.symbols, columns=genecols, keytype="SYMBOL")

bg.ID = genetable$ENTREZID[match(bg.genes, as.character(genetable$SYMBOL))]
bg.ID = bg.ID[!is.na(bg.ID)]
write.table(data.frame(bg.ID),bg.file, col.names=F, row.names=F, sep="\t", quote=F)

up.ID = genetable$ENTREZID[match(up.genes, as.character(genetable$SYMBOL))]
up.ID = up.ID[!is.na(up.ID)]
write.table(data.frame(up.ID),up.file, col.names=F, row.names=F, sep="\t", quote=F)

down.ID = genetable$ENTREZID[match(down.genes, as.character(genetable$SYMBOL))]
down.ID = down.ID[!is.na(down.ID)]
write.table(data.frame(down.ID),down.file, col.names=F, row.names=F, sep="\t", quote=F)