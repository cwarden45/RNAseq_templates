normalizeTotalExpression = function (geneExpr, totalReads) {
	return(geneExpr / totalReads)
}#end def normalizeTotalExpression

count.defined.values = function(arr, expr.cutoff)
{
	sig.figures = 1
	if (expr.cutoff > 0)
		sig.figures = 0
	expr.cutoff = round(expr.cutoff, digits=sig.figures)
	arr = round(arr, digits=sig.figures)
	return(length(arr[arr > expr.cutoff]))
}#end def count.values

trimmed.counts = function(counts, min.percent, max.percent)
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

reduced.exon.length = function(text.arr){
	#print(text.arr)
	mat = matrix(unlist(strsplit(text.arr,split="\t")),ncol=4, byrow=TRUE)
	#print(mat)
	exon.chr = mat[,1]
	non.canonical = exon.chr[grep("_",exon.chr)]
	if ((length(non.canonical) > 0)&(length(non.canonical) != nrow(mat))){
		if(nrow(mat) - length(non.canonical) == 1){
			mat = mat[-grep("_",exon.chr),]
			#print(mat)
			return(as.numeric(mat[3]) - as.numeric(mat[2]))
		}else{
			mat = mat[-grep("_",exon.chr),]
			#print(mat)
			exon.chr = mat[,1]
		}
	}#end if ((length(non.canonical) > 0)&(length(non.canonical) != nrow(mat)))
	exon.start = as.numeric(mat[,2])
	exon.stop = as.numeric(mat[,3])
	exon.strand = mat[,4]
	
	#gr = reduce(GRanges(Rle(exon.chr),
    	#		IRanges(start=exon.start, end=exon.stop),
    	#		Rle(strand(exon.strand))))
	#reduced.start = start(gr)
	#reduced.stop = end(gr)
	
	ir = reduce(IRanges(start=exon.start, end=exon.stop))
	reduced.start = start(ir)
	reduced.stop = end(ir)
	merged.exon.length = reduced.stop - reduced.start
	return(sum(merged.exon.length))
}#end def reduced.exon.length

param.table = read.table("parameters.txt", header=T, sep="\t")
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
min.expression = as.numeric(as.character(param.table$Value[param.table$Parameter == "fpkm_expression_cutoff"]))
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
htseq.anno.folder = as.character(param.table$Value[param.table$Parameter == "HTseq_input_folder"])
sample.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
total.reads.file = as.character(param.table$Value[param.table$Parameter == "total_counts_file"])
exonic.stat.file = as.character(param.table$Value[param.table$Parameter == "aligned_stats_file"])
counts.file = as.character(param.table$Value[param.table$Parameter == "counts_file"])
rpkm.file = as.character(param.table$Value[param.table$Parameter == "fpkm_file"])
aligned.type = as.character(param.table$Value[param.table$Parameter == "FPKM_norm"])
full.annotation.file = paste(htseq.anno.folder,"\\TxDb_",genome,"_exon_annotations.txt",sep="")

library(GenomicRanges)
setwd(output.folder)

sample.table = read.table(sample.file, header=T, sep="\t")
sampleID = as.character(sample.table$sampleID)
sample.label = as.character(sample.table$userID)
dash.flag = grep("-",sample.label)
if(length(dash.flag) > 0){
	print(paste(paste(sample.label[dash.flag],collapse=",")," samples labels have dashes in their labels",sep=""))
}
num.flag = grep("^[0-9]",sample.label)
if(length(num.flag) > 0){
	print(paste(paste(sample.label[num.flag],collapse=",")," samples labels start with numbers",sep=""))
}

if((length(dash.flag) > 0)|(length(num.flag) > 0)){
	stop("Please make sure sample labels do not have dashes and do not start with numbers")
}

count.files = as.character(sample.table$HTseq.file)

total.reads.table = read.table(total.reads.file, header=T, sep="\t")
totalID = as.character(total.reads.table$Sample)
total.reads.table = total.reads.table[match(sampleID, totalID),]
print(total.reads.table$Total.Reads)

temp.file = count.files[[1]]
temp.table = read.table(temp.file, sep="\t", header=F)
genes = as.character(temp.table[[1]])

count.mat = matrix(nrow=nrow(temp.table), ncol=length(sampleID))
colnames(count.mat) = sample.label

matched.genes = c()

for (i in 1:length(sampleID)){
	temp.file = count.files[[i]]
	temp.table = read.table(temp.file, sep="\t", header=F)
	temp.genes = as.character(temp.table[[1]])
	
	if(length(matched.genes) != 0){
		matched.genes = temp.genes[match(matched.genes, temp.genes, nomatch=0)]
	} else if(!identical(temp.genes, genes)){
		print("Genes not in same order - most likely, different quantification method was used for different samples")
		userAns = readline(prompt="Do you wish to proceed with a subset of matched gene symbols? (y/n): ")
		userAns = tolower(substr(userAns, 1, 1))
		if (userAns != "y"){
			stop("Please re-run htseq-count with all of your samples")
		} else {
			matched.genes = temp.genes[match(genes, temp.genes, nomatch=0)]
		}#end else
	} else {
		count.mat[,i] = temp.table[[2]]
	}#end else
}#end for (i in 1:length(sampleID))

if (length(matched.genes) != 0){
	count.mat = matrix(nrow=length(matched.genes), ncol=length(sampleID))
	colnames(count.mat) = sample.label
	
	for (i in 1:length(sampleID)){
	temp.file = count.files[[i]]
	temp.table = read.table(temp.file, sep="\t", header=F)
	temp.genes = as.character(temp.table[[1]])
	temp.counts = temp.table[[2]]
	
	count.mat[,i] = temp.counts[match(matched.genes, temp.genes, nomatch=0)]
	}#end for (i in 1:length(sampleID))
	
	genes = matched.genes
}#end if (length(matched.genes) != 0)

irrelevant.counts = c("__no_feature", "__ambiguous", "__too_low_aQual","__not_aligned","__alignment_not_unique")
extra.stats = count.mat[match(irrelevant.counts, genes),]
count.mat = count.mat[-match(irrelevant.counts, genes),]
genes = genes[-match(irrelevant.counts, genes)]
rownames(count.mat) = genes

exon.info = read.table(full.annotation.file, header=T, sep="\t")
print(dim(exon.info))
exon.length = exon.info$width

gene.symbol = as.character(levels(as.factor(as.character(exon.info$symbol))))
gene.chr = exon.info$chr[match(gene.symbol,exon.info$symbol)]
gene.strand= exon.info$strand[match(gene.symbol,exon.info$symbol)]
gene.description = exon.info$description[match(gene.symbol,exon.info$symbol)]
gene.start= tapply(as.numeric(exon.info$start), as.character(exon.info$symbol), min)
gene.end = tapply(as.numeric(exon.info$end), as.character(exon.info$symbol), max)
pre.gr = paste(exon.info$chr,exon.info$start,exon.info$end,exon.info$strand,sep="\t")
gene.length.kb = tapply(pre.gr, as.character(exon.info$symbol), reduced.exon.length) / 1000

gene.info = data.frame(symbol = gene.symbol, chr = gene.chr, start = gene.start, stop = gene.end, length.kb = gene.length.kb,
						strand = gene.strand, description=gene.description)
common.genes = gene.info$symbol[match(rownames(count.mat), gene.info$symbol,nomatch=0)]
gene.info = gene.info[match(common.genes, gene.info$symbol, nomatch=0),]
count.mat = count.mat[match(common.genes, rownames(count.mat), nomatch=0),]
gene.length.kb = gene.info$length.kb

annotated.rpkm = data.frame(gene.info, count.mat)
write.table(annotated.rpkm, file = counts.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder,counts.file,sep="/")
write.table(annotated.rpkm, file=result.file, row.names=F, quote=F, sep="\t")

exonic.stat.table = read.table(exonic.stat.file, header=T, sep="\t")
exonicID = as.character(exonic.stat.table$Sample)
exonic.stat.table = exonic.stat.table[match(sampleID, exonicID),]
total.reads = as.numeric(total.reads.table$Total.Reads)

if(aligned.type == "aligned"){
	aligned.reads = as.numeric(exonic.stat.table$aligned.reads)
} else if(aligned.type =="quantified"){
	aligned.reads=apply(count.mat, 2, sum)
}else {
	stop("Print RPKM_norm must be either 'aligned' or 'quantified'")
}#end else

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
RPKM = round(log2(t(apply(rpk, 1, normalizeTotalExpression, totalReads = total.million.aligned.reads)) + min.expression), digits=2)
colnames(RPKM) = sample.label

trimmed.percent = apply(count.mat, 2, trimmed.counts, min.percent=0.3, max.percent=0.95)

expressed.gene.counts = apply(RPKM, 2, count.defined.values, expr.cutoff = log2(min.expression))
percent.expressed.genes = round( 100 * expressed.gene.counts / nrow(RPKM), digits=1)
coverage.table = data.frame(Sample = sample.label, total.reads = total.reads,
							aligned.reads=aligned.reads, percent.aligned.reads =paste(percent.aligned.reads,"%",sep=""),
							htseq.nofeature.reads =intergenic.reads, percent.unique.exonic.reads = paste(percent.exonic.reads,"%",sep=""),
							Expressed.Genes = expressed.gene.counts, Percent.Expressed.Genes = paste(percent.expressed.genes,"%",sep=""),
							trimmed.percent=paste(trimmed.percent,"%",sep=""))
write.table(coverage.table, file="gene_coverage_stats.txt", quote=F, row.names=F, sep="\t")

	
#tables have different file formats for downstream R analysis versus other applications that involve parsing text files
#	--> don't set the Result folder to the working directory, or you may skip genes during DEG analysis
annotated.rpkm = data.frame(gene.info, RPKM)
write.table(annotated.rpkm, file = rpkm.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder, rpkm.file, sep="/")
write.table(annotated.rpkm, file=result.file, row.names=F, quote=F, sep="\t")
