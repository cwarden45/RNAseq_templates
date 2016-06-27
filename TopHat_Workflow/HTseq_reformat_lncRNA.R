normalizeTotalExpression <- function (geneExpr, totalReads) {
	return(geneExpr / totalReads)
}#end def normalizeTotalExpression

count.defined.values <- function(arr, expr.cutoff)
{
	return(length(arr[arr > expr.cutoff]))
}#end def count.values

trimmed.counts <- function(counts, min.percent, max.percent)
{
	total.counts = sum(counts)
	counts = counts[order(counts)]
	min.index = min.percent * length(counts)
	max.index = max.percent * length(counts)
	counts = counts[min.index:max.index]
	trimmed.counts = sum(counts)
	trimmed.percent = round(100 * trimmed.counts/total.counts, digits=1)
	return(trimmed.percent)
}#end def count.values

aligned.trimmed.counts <- function(counts, min.percent, max.percent)
{
	counts = counts[order(counts)]
	min.index = min.percent * length(counts)
	max.index = max.percent * length(counts)
	counts = counts[min.index:max.index]
	return(sum(counts))
}#end def aligned.trimmed.counts

param.table = read.table("parameters_lncRNA.txt", header=T, sep="\t")
lncRNA.source=as.character(param.table$Value[param.table$Parameter == "lncRNA_ann_name"])
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
min.expression = as.numeric(as.character(param.table$Value[param.table$Parameter == "rpkm_expression_cutoff"]))
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
htseq.anno.folder = as.character(param.table$Value[param.table$Parameter == "HTseq_input_folder"])
sample.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
total.reads.file = as.character(param.table$Value[param.table$Parameter == "total_counts_file"])
exonic.stat.file = as.character(param.table$Value[param.table$Parameter == "aligned_stats_file"])
counts.file = as.character(param.table$Value[param.table$Parameter == "counts_file"])
rpkm.file = as.character(param.table$Value[param.table$Parameter == "rpkm_file"])
full.annotation.file = paste(htseq.anno.folder,"\\",lncRNA.source,"\\",genome,"\\",genome,"_gene_info.txt",sep="")

setwd(output.folder)

sample.table = read.table(sample.file, header=T, sep="\t")
sampleID = as.character(sample.table$sampleID)
sample.label = as.character(sample.table$userID)
count.files = as.character(sample.table$HTseq.file)

count.files = as.character(sample.table$HTseq.file)
print(count.files)

total.reads.table = read.table(total.reads.file, header=T, sep="\t")
totalID = as.character(total.reads.table$Sample)
total.reads.table = total.reads.table[match(sampleID, totalID),]
print(total.reads.table$Total.Reads)

temp.file = count.files[[1]]
temp.table = read.table(temp.file, sep="\t", header=F)
genes = as.character(temp.table[[1]])

count.mat = matrix(nrow=nrow(temp.table), ncol=length(sampleID))
colnames(count.mat) = sample.label

for (i in 1:length(sampleID)){
	temp.file = count.files[[i]]
	temp.table = read.table(temp.file, sep="\t", header=F)
	if(!identical(as.character(temp.table[[1]]), genes)){
		stop("Need to revise code! Genes not in same order.")
	}

	count.mat[,i] = temp.table[[2]]
}#end for (i in 1:length(sampleIDs))

total.aligned.reads = apply(count.mat, 2, sum)

irrelevant.counts = c("__no_feature", "__ambiguous", "__too_low_aQual","__not_aligned","__alignment_not_unique")
extra.stats = count.mat[match(irrelevant.counts, genes),]
count.mat = count.mat[-match(irrelevant.counts, genes),]
genes = genes[-match(irrelevant.counts, genes)]
rownames(count.mat) = genes

gene.info  = read.table(full.annotation.file, header=T, sep="\t")
print(dim(gene.info))
gene.info = gene.info[match(genes, gene.info$symbol),]
print(dim(gene.info))

gene.length.kb = gene.info$tx.max.length / 1000

annotated.rpkm = data.frame(gene.info, count.mat)
write.table(annotated.rpkm, file = counts.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder,"lncRNA",counts.file,sep="/")
write.table(annotated.rpkm, file=result.file, row.names=F, quote=F, sep="\t")

exonic.stat.table = read.table(exonic.stat.file, header=T, sep="\t")
exonicID = as.character(exonic.stat.table$Sample)
exonic.stat.table = exonic.stat.table[match(sampleID, exonicID),]
total.reads = as.numeric(total.reads.table$Total.Reads)
aligned.reads = as.numeric(exonic.stat.table$aligned.reads)
percent.aligned.reads = round(100 * aligned.reads / total.reads, digits=1)

intergenic.reads = extra.stats[irrelevant.counts == "__no_feature", ]
exonic.reads = apply(count.mat, 2, sum)
unique.reads = exonic.reads + intergenic.reads
percent.exonic.reads = round(100 * exonic.reads / unique.reads, digits=1)

total.million.aligned.reads = aligned.reads / 1000000
print(total.million.aligned.reads)

rpk = matrix(ncol=ncol(count.mat), nrow=nrow(count.mat))
for (i in 1:ncol(count.mat)){
	counts = as.numeric(count.mat[,i])
	temp.rpk = counts / as.numeric(gene.length.kb)
	rpk[,i] = temp.rpk 
}
RPKM = log2(t(apply(rpk, 1, normalizeTotalExpression, totalReads = total.million.aligned.reads)) + min.expression)
colnames(RPKM) = sample.label

trimmed.percent = apply(count.mat, 2, trimmed.counts, min.percent=0.3, max.percent=0.95)

expressed.gene.counts = apply(RPKM, 2, count.defined.values, expr.cutoff = log2(min.expression))
percent.expressed.genes = round( 100 * expressed.gene.counts / nrow(RPKM), digits=1)
coverage.table = data.frame(Sample = sample.label, total.reads = total.reads,
							aligned.reads=aligned.reads, percent.aligned.reads =percent.aligned.reads,
							htseq.nofeature.reads =intergenic.reads, percent.unique.exonic.reads = percent.exonic.reads,
							Expressed.Genes = expressed.gene.counts, Percent.Expressed.Genes = percent.expressed.genes,
							trimmed.percent=trimmed.percent)
write.table(coverage.table, file="gene_coverage_stats_lncRNA.txt", quote=F, row.names=F, sep="\t")

annotated.rpkm = data.frame(gene.info, RPKM)
write.table(annotated.rpkm, file = rpkm.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder, "lncRNA",rpkm.file, sep="/")
write.table(annotated.rpkm, file=result.file, row.names=F, quote=F, sep="\t")