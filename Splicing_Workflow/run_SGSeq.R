max.genes = 500

count.na = function(subset.arr){
	return(length(subset.arr[is.na(subset.arr)]))
}#end count.na

param.table = read.table("parameters.txt", header=T, sep="\t")
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
sample.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
fdr.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "fdr_cutoff"]))
comp.name=as.character(param.table$Value[param.table$Parameter == "comp_name"])
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
threads = as.character(param.table$Value[param.table$Parameter == "Threads"])

gene.result = read.table(paste(user.folder,"/DSG/sig_gene_stats_",comp.name,".txt",sep=""), head=T, sep="\t")
sig.symbol = gene.result$geneName[gene.result$geneWisePadj < fdr.cutoff]

if(length(sig.symbol) > max.genes){
	print(paste("There are ",length(sig.symbol)," genes with FDR < ",fdr.cutoff,sep=""))
	stop("Please impose stricter JunctionSeq FDR cutoff (recommended) or change `max.genes` value in template")
}

library(SGSeq)
if(genome == "hg19"){
	library(TxDb.Hsapiens.UCSC.hg19.knownGene)
	txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

	library(org.Hs.eg.db)
	orgdb = org.Hs.eg.db
	
	library(BSgenome.Hsapiens.UCSC.hg19)
	varEffDb = Hsapiens
	seqlevelsStyle(varEffDb) = "NCBI"
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

geneIDs = genetable$ENTREZID[match(sig.symbol, genetable$SYMBOL)]
sig.symbol=sig.symbol[!is.na(geneIDs)]
geneIDs=geneIDs[!is.na(geneIDs)]
geneID.mapping = genetable[match(geneIDs,genetable$ENTREZID),]
rm(genetable)
gene.chr = txtable$EXONCHROM[match(geneIDs,txtable$GENEID, nomatch=0)]

gene.strand = txtable$EXONSTRAND[match(geneIDs,txtable$GENEID, nomatch=0)]

gene.start = tapply(txtable$EXONSTART, txtable$GENEID, min)
gene.start = gene.start[match(geneIDs, names(gene.start))]

gene.stop = tapply(txtable$EXONEND, txtable$GENEID, max)
gene.stop = gene.stop[match(geneIDs, names(gene.stop))]

#if your gene list is long enough, you'll probably need to filter non-canonical chromosome aligments
unaligned.genes = sig.symbol[grep("_",gene.chr)]
if(length(unaligned.genes) > 0){
	gene.strand=gene.strand[-grep("_",gene.chr)]
	gene.start=gene.start[-grep("_",gene.chr)]
	gene.stop=gene.stop[-grep("_",gene.chr)]
	gene.chr=gene.chr[-grep("_",gene.chr)]
}#end if(length(unaligned.genes) > 0

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

#novel exons / junctions will have empty annotation values

SGS.folder = paste(user.folder,"/SGSeq",sep="")
dir.create(SGS.folder)

feature.file = paste(SGS.folder,"/features_",comp.name,".txt",sep="")
feature.table = data.frame(rowRanges(sgfc_pred))
write.table(feature.table, feature.file, quote=F, sep="\t", row.names=F)

feature.table = read.table(feature.file, head=T, sep="\t")
feature.table$txName[feature.table$txName == "character(0)"]=NA

novel.status = tapply(feature.table$txName, feature.table$geneName, count.na)
novel.geneID = names(novel.status)
novel.status = as.character(novel.status)
novel.status[novel.status == "1"]="some novel"
novel.status[novel.status == "0"]="all known"
novel.symbol = geneID.mapping$SYMBOL[match(novel.geneID,geneID.mapping$ENTREZID)]

feature.junction.table = feature.table[feature.table$type == "J",]
novel.junction.status = tapply(feature.junction.table$txName, feature.junction.table$geneName, count.na)
novel.junction.status = novel.junction.status[match(novel.geneID, names(novel.junction.status))]
novel.junction.status = as.character(novel.junction.status)
novel.junction.status[novel.junction.status == "1"]="some novel"
novel.junction.status[novel.junction.status == "0"]="all known"

feature.exon.table = feature.table[feature.table$type == "E",]
novel.exon.status = tapply(feature.exon.table$txName, feature.exon.table$geneName, count.na)
novel.exon.status = novel.exon.status[match(novel.geneID, names(novel.exon.status))]
novel.exon.status = as.character(novel.exon.status)
novel.exon.status[novel.exon.status == "1"]="some novel"
novel.exon.status[novel.exon.status == "0"]="all known"

feature.donor.table = feature.table[feature.table$type == "D",]
novel.donor.status = tapply(feature.donor.table$txName, feature.donor.table$geneName, count.na)
novel.donor.status = novel.donor.status[match(novel.geneID, names(novel.donor.status))]
novel.donor.status = as.character(novel.donor.status)
novel.donor.status[novel.donor.status == "1"]="some novel"
novel.donor.status[novel.donor.status == "0"]="all known"

feature.acceptor.table = feature.table[feature.table$type == "A",]
novel.acceptor.status = tapply(feature.acceptor.table$txName, feature.acceptor.table$geneName, count.na)
novel.acceptor.status = novel.acceptor.status[match(novel.geneID, names(novel.acceptor.status))]
novel.acceptor.status = as.character(novel.acceptor.status)
novel.acceptor.status[novel.acceptor.status == "1"]="some novel"
novel.acceptor.status[novel.acceptor.status == "0"]="all known"

novel.gene.table = data.frame(novel.symbol, novel.geneID, novel.overall.status = novel.status,
								novel.junction.status, novel.exon.status, novel.donor.status, novel.acceptor.status)
gene.file = paste(SGS.folder,"/genes_novel_feature_status_",comp.name,".txt",sep="")
write.table(novel.gene.table, gene.file, quote=F, sep="\t", row.names=F)

rm(geneID.mapping)

data.image = paste(SGS.folder,"/",comp.name,".Rdata",sep="")
save.image(data.image)

###uncomment for variant effect predictions###
#vep = predictVariantEffects(sgv_pred, txdb, varEffDb)

###uncomment for feature counts###
#sgvc = getSGVariantCounts(sgv_pred, sample_info = si)
#x = counts(sgvc)
#vid = variantID(sgvc)
#eid = eventID(sgvc)

###uncomment for inclusion percentage###
#variantFreq(sgvc_pred)

###plot commands###
plot.folder = paste(SGS.folder,"/plots",sep="")
dir.create(SGS.folder)

for (i in 1:length(sig.symbol)){
	gene = as.character(sig.symbol[i])
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
