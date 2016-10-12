genes = c()

param.table = read.table("parameters.txt", header=T, sep="\t")
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
sample.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
threads = as.character(param.table$Value[param.table$Parameter == "Threads"])

library(SGSeq)
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
	stop(paste("Need to add annotation database for", genome))
}

gene.symbols = keys(orgdb, keytype="SYMBOL")
genecols = c("SYMBOL", "ENTREZID","GENENAME")
genetable = select(orgdb, keys=gene.symbols, columns=genecols, keytype="SYMBOL")
print(dim(genetable))

gene.keys = keys(txdb, keytype="GENEID")
txcols = c("EXONCHROM", "EXONSTART","EXONEND","GENEID","EXONID","EXONSTRAND")
txtable = select(txdb, keys=gene.keys, columns=txcols, keytype="GENEID")
print(dim(txtable))
gene.length = txtable$EXONEND - txtable$EXONSTART + 1
txtable = data.frame(txtable, LENGTH=gene.length)

geneIDs = genetable$ENTREZID[match(genes, genetable$SYMBOL)]
rm(genetable)
gene.chr = txtable$EXONCHROM[match(geneIDs,txtable$GENEID, nomatch=0)]
gene.strand = txtable$EXONSTRAND[match(geneIDs,txtable$GENEID, nomatch=0)]
gene.start = tapply(txtable$EXONSTART, txtable$GENEID, min)
gene.start = gene.start[match(geneIDs, names(gene.start))]
gene.stop = tapply(txtable$EXONEND, txtable$GENEID, max)
gene.stop = gene.stop[match(geneIDs, names(gene.stop))]
rm(txtable)
gr = GRanges(Rle(gene.chr),
              IRanges(start=gene.start, end=gene.stop),
              Rle(strand(gene.strand)))

si = read.table(sample.file, sep="", header=T)
si$file_bam = as.character(si$file_bam)
si$sample_name = as.character(si$userID)
print(si)
txf_ucsc = convertToTxFeatures(txdb)

sgfc_pred = analyzeFeatures(si, which=gr, cores=as.numeric(threads))
sgfc_pred = annotate(sgfc_pred, txf_ucsc)

SGS.folder = paste(user.folder,"/SGSeq_plots",sep="")
dir.create(SGS.folder)

for (i in 1:length(genes)){
	gene = genes[i]
	geneID = as.character(geneIDs[i])
	splicing.heatmap = paste(SGS.folder,"/",gene,"_junction_heatmap.png",sep="")
	png(splicing.heatmap)
	gene.stats = plotFeatures(sgfc_pred, geneName=geneID, color_novel="red", toscale = "none", include="junctions")
	dev.off()

	splice.graph = paste(SGS.folder,"/",gene,"_splice_graph.png",sep="")
	png(splice.graph)
	par(mfrow = c(nrow(si)+1, 1), mar = c(1, 3, 1, 1))
	plotSpliceGraph(rowRanges(sgfc_pred), geneName= geneID, toscale = "none", color_novel = "red")
	for (j in 1:nrow(si)) {
  	plotCoverage(sgfc_pred[, j], geneName= geneID, toscale = "none")
	}
	dev.off()	
}

