count.na.values = function(arr)
{
	return(length(arr[is.na(arr)]))
}#end def count.values

param.table = read.table("parameters.txt", header=T, sep="\t")

comp.name=as.character(param.table$Value[param.table$Parameter == "comp_name"])
count.folder=as.character(param.table$Value[param.table$Parameter == "QoRTs_Merged_Folder"])
deg.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "deg_groups"]), split=","))
trt.group = as.character(param.table$Value[param.table$Parameter == "treatment_group"])
goseq.flag = as.character(param.table$Value[param.table$Parameter == "run_goseq"])
threads = as.character(param.table$Value[param.table$Parameter == "Threads"])
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
pvalue.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "pvalue_cutoff"]))
fdr.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "fdr_cutoff"]))
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
gene.mapping.file = as.character(param.table$Value[param.table$Parameter == "GENCODE_Gene_Info"])

library(JunctionSeq)

setwd(output.folder)

gff.file = paste(count.folder,"/withNovel.forJunctionSeq.gff.gz",sep="")

sample.description.table = read.table(sample.description.file,head=T, sep="\t")
sampleIDs = as.character(sample.description.table$sample.ID)
sample.labels = as.character(sample.description.table$userID)

design = sample.description.table[,deg.groups]
if (length(deg.groups) == 1){
	sampleIDs = sampleIDs[!is.na(design)]
	sample.labels = sample.labels[!is.na(design)]
	design = data.frame(condition = as.factor(design[!is.na(design),]))
} else {
	deg.grp.na.counts = apply(design, 1, count.na.values)
	sampleIDs = sampleIDs[deg.grp.na.counts == 0]
	sample.labels = sample.labels[deg.grp.na.counts == 0]
	design = design[deg.grp.na.counts == 0,]
	for (i in 1:ncol(design)){
		design[,i] = as.factor(design[,i])
	}
}
new.col.names = colnames(design)
new.col.names[1]= "condition"
colnames(design) = new.col.names

gene.info.table = read.table(gene.mapping.file, head=T, sep="\t")
geneID.to.symbol.file = data.frame(ENS_GENEID = gene.info.table$Ensembl.GeneID, geneSymbol=gene.info.table$Gene.Name)
write.table(geneID.to.symbol.file, "ensid.2.symbol.txt", quote=F, sep="\t", row.names=F)

countFiles = paste(count.folder,"/",sampleIDs,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz",sep="")

jscs = readJunctionSeqCounts(countfiles = countFiles, samplenames = sample.labels, design = design,
								flat.gff.file = gff.file, gene.names = geneID.to.symbol.file)
jscs = estimateJunctionSeqSizeFactors(jscs)
jscs = estimateJunctionSeqDispersions(jscs, nCores = as.numeric(threads))
jscs = fitJunctionSeqDispersionFunction(jscs)
jscs = testForDiffUsage(jscs, nCores = as.numeric(threads))
jscs = estimateEffectSizes(jscs, nCores = as.numeric(threads))

raw.output.folder = paste(count.folder,"/",comp.name,"/",sep="")
dir.create(raw.output.folder)
plot.folder = paste(user.folder,"/DSG/",comp.name,"/",sep="")
dir.create(plot.folder)
								
writeCompleteResults(jscs, outfile.prefix=raw.output.folder, save.jscs = TRUE,
						FDR.threshold = fdr.cutoff, gzip.output = FALSE)
					
buildAllPlots(jscs=jscs, outfile.prefix = plot.folder,
				use.plotting.device = "png",FDR.threshold = 0.01)