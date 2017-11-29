genome = "mm9"
annotation.table = paste("TxDb_",genome,"_exon_annotations.txt",sep="")

if(genome == "hg38"){
	library(TxDb.Hsapiens.UCSC.hg38.knownGene)
	txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
	
	library(org.Hs.eg.db)
	orgdb = org.Hs.eg.db
} else if(genome == "hg19"){
	library(TxDb.Hsapiens.UCSC.hg19.knownGene)
	txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

	library(org.Hs.eg.db)
	orgdb = org.Hs.eg.db
} else if(genome == "hg18"){
	library(TxDb.Hsapiens.UCSC.hg18.knownGene)
	txdb = TxDb.Hsapiens.UCSC.hg18.knownGene

	library(org.Hs.eg.db)
	orgdb = org.Hs.eg.db
} else if (genome == "mm8"){
	library(TxDb.Mmusculus.UCSC.mm8.knownGene)
	txdb = TxDb.Mmusculus.UCSC.mm8.knownGene

	library(org.Mm.eg.db)
	orgdb = org.Mm.eg.db
} else if (genome == "mm9"){
	library(TxDb.Mmusculus.UCSC.mm9.knownGene)
	txdb = TxDb.Mmusculus.UCSC.mm9.knownGene

	library(org.Mm.eg.db)
	orgdb = org.Mm.eg.db
} else if (genome == "mm10"){
	library(TxDb.Mmusculus.UCSC.mm10.knownGene)
	txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

	library(org.Mm.eg.db)
	orgdb = org.Mm.eg.db
} else {
	stop("Need to add annotations for reference!")
}
 
gene.keys = keys(txdb, keytype="GENEID")
txcols = c("EXONCHROM", "EXONSTART","EXONEND","GENEID","EXONID","EXONSTRAND")
txtable = select(txdb, keys=gene.keys, columns=txcols, keytype="GENEID")
print(dim(txtable))
exon.length = txtable$EXONEND - txtable$EXONSTART + 1
txtable = data.frame(txtable, LENGTH=exon.length)

gene.symbols = keys(orgdb, keytype="SYMBOL")
genecols = c("SYMBOL", "ENTREZID","GENENAME")
genetable = select(orgdb, keys=gene.symbols, columns=genecols, keytype="SYMBOL")
print(dim(genetable))

temp.annotation = merge(txtable, genetable, by.x = "GENEID", by.y = "ENTREZID")
print(dim(temp.annotation))
temp.annotation = unique(temp.annotation)
print(dim(temp.annotation))
combined.annotation = data.frame(chr = temp.annotation$EXONCHROM,
								start = temp.annotation$EXONSTART,
								end = temp.annotation$EXONEND,
								width = temp.annotation$LENGTH,
								symbol = temp.annotation$SYMBOL,
								strand = as.character(temp.annotation$EXONSTRAND),
								description = temp.annotation$GENENAME)
ID=paste(combined.annotation$chr,combined.annotation$symbol,combined.annotation$strand,temp.annotation$EXONID, sep=".")
#print(dim(combined.annotation))
#combined.annotation=combined.annotation[match(unique(ID),ID),]
#ID=ID[match(unique(ID),ID)]
#print(dim(combined.annotation))
rownames(combined.annotation) = ID

write.table(combined.annotation, annotation.table, sep="\t", row.names=F)
