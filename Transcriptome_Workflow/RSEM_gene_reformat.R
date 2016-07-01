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
}#end def trimmed.counts

count.defined.values <- function(arr, expr.cutoff)
{
	return(length(arr[arr > expr.cutoff]))
}#end def count.values

gene.summary <- function(gene.arr, gene.symbols)
{
	result = tapply(as.numeric(gene.arr), as.character(gene.symbols), sum)
}#end def gene.summary

param.table = read.table("parameters.txt", header=T, sep="\t")
comp.name=as.character(param.table$Value[param.table$Parameter == "comp_name"])
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
min.expression = as.numeric(as.character(param.table$Value[param.table$Parameter == "tpm_expression_cutoff"]))
gene.annotation.file = as.character(param.table$Value[param.table$Parameter == "gene_annotation_file"])
sample.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
total.reads.file = as.character(param.table$Value[param.table$Parameter == "total_counts_file"])
exonic.stat.file = as.character(param.table$Value[param.table$Parameter == "aligned_stats_file"])
gene.counts.file = as.character(param.table$Value[param.table$Parameter == "transcript_counts_file"])
gene.tpm.file = as.character(param.table$Value[param.table$Parameter == "transcript_tpm_file"])
gene.counts.file = as.character(param.table$Value[param.table$Parameter == "gene_counts_file"])
gene.tpm.file = as.character(param.table$Value[param.table$Parameter == "gene_tpm_file"])

setwd(output.folder)

sample.table = read.table(sample.file, header=T, sep="\t")
sampleID = as.character(sample.table$sampleID)
sample.label = as.character(sample.table$userID)
dash.flag = grep("-",sample.label)
if(length(dash.flag) > 0){
	print(paste(paste(sample.label[dash.flag],collapse=",")," samples labels have dashes in their labels",sep=""))
}
num.flag = grep("^[0-9]",sample.label)
if(length(dash.flag) > 0){
	print(paste(paste(sample.label[num.flag],collapse=",")," samples labels start with numbers",sep=""))
}

if((length(dash.flag) > 0)|(length(num.flag) > 0)){
	stop("Please make sure sample labels do not have dashes and do not start with numbers")
}

count.files = as.character(sample.table$quant.file)

temp.file = count.files[[1]]
temp.table = read.table(temp.file, sep="\t", header=T)
genes = as.character(temp.table$gene_id)

gene.counts.mat = matrix(nrow=nrow(temp.table), ncol=length(sampleID))
colnames(gene.counts.mat) = sample.label

gene.tpm.mat = matrix(nrow=nrow(temp.table), ncol=length(sampleID))
colnames(gene.tpm.mat) = sample.label

matched.transcripts = c()

for (i in 1:length(sampleID)){
	print(sampleID[i])
	temp.file = count.files[[i]]
	temp.table = read.table(temp.file, sep="\t", header=T)
	temp.transcripts = as.character(temp.table$Name)
	
	if(length(matched.transcripts) != 0){
		matched.transcripts = temp.transcripts[match(matched.transcripts, temp.transcripts, nomatch=0)]
	} else if(!identical(temp.transcripts, transcripts)){
		print("transcripts not in same order - most likely, different quantification method was used for different samples")
		userAns = readline(prompt="Do you wish to proceed with a subset of matched gene IDs? (y/n): ")
		userAns = tolower(substr(userAns, 1, 1))
		if (userAns != "y"){
			stop("Please re-run RSEM with all of your samples")
		} else {
			matched.transcripts = temp.transcripts[match(transcripts, temp.transcripts, nomatch=0)]
		}#end else
	} else {
		gene.counts.mat[,i] = temp.table$expected_count
		gene.tpm.mat[,i] = temp.table$TPM
	}#end else
}#end for (i in 1:length(sampleID))

if (length(matched.transcripts) != 0){
	count.mat = matrix(nrow=length(matched.transcripts), ncol=length(sampleID))
	colnames(count.mat) = sample.label
	
	for (i in 1:length(sampleID)){
		print(sampleID[i])
		temp.file = count.files[[i]]
		temp.table = read.table(temp.file, sep="\t", header=T)
		temp.transcripts = as.character(temp.table$Name)
		temp.counts = temp.table$expected_count
		temp.tpm = temp.table$TPM
		
		gene.counts.mat[,i] = temp.counts[match(matched.transcripts, temp.transcripts, nomatch=0)]
		gene.tpm.mat[,i] = temp.tpm[match(matched.transcripts, temp.transcripts, nomatch=0)]
	}#end for (i in 1:length(sampleID))
	
	transcripts = matched.transcripts
}#end if (length(matched.transcripts) != 0)

gene.counts.mat = round(gene.counts.mat)
gene.tpm.mat = round(gene.tpm.mat, digits=2)

gene.info = read.table(gene.annotation.file, header=T, sep="\t")

gene.symbol = as.character(levels(as.factor(as.character(gene.info$Gene.Name))))
geneID= gene.info$NCBI.GeneID[match(gene.symbol,as.character(gene.info$Gene.Name))]
Gene.Type = gene.info$Gene.Type[match(gene.symbol,as.character(gene.info$Gene.Name))]
gene.info = data.frame(Gene.Name = gene.symbol, geneID = geneID, Gene.Type = Gene.Type)

annotated.gene.tpm = data.frame(gene.info, gene.tpm.mat)
write.table(annotated.gene.tpm, file = gene.tpm.file, sep="\t", row.names=F, quote=T)

#result.file = paste(user.folder, gene.tpm.file, sep="/")
#write.table(annotated.gene.tpm, file=result.file, row.names=F, quote=F, sep="\t")

annotated.gene.counts = data.frame(gene.info, gene.counts.mat)
write.table(annotated.gene.counts, file = gene.counts.file, sep="\t", row.names=F, quote=T)

#result.file = paste(user.folder, gene.counts.mat, sep="/")
#write.table(annotated.gene.counts, file=result.file, row.names=F, quote=F, sep="\t")

gene.reads = apply(gene.counts.mat, 2, sum)

gene.trimmed.percent = apply(gene.counts.mat, 2, trimmed.counts, min.percent=0.3, max.percent=0.95)

expressed.gene.counts = apply(gene.tpm.mat, 2, count.defined.values, expr.cutoff = min.expression)
percent.expressed.genes = round( 100 * expressed.gene.counts / nrow(gene.tpm.mat), digits=1)

protein.coding.gene.counts = gene.counts.mat[gene.info$Gene.Type == "protein_coding", ]
protein.coding.gene.tpm = gene.tpm.mat[gene.info$Gene.Type == "protein_coding", ]

gene.trimmed.protein.coding.percent = apply(protein.coding.gene.counts, 2, trimmed.counts, min.percent=0.3, max.percent=0.95)

expressed.protein.coding.gene.counts = apply(protein.coding.gene.tpm, 2, count.defined.values, expr.cutoff = min.expression)
percent.protein.coding.expressed.genes = round( 100 * expressed.protein.coding.gene.counts / nrow(protein.coding.gene.tpm), digits=1)
coverage.table = data.frame(Sample = sample.label,
							gene.reads=gene.reads,
							Expressed.Genes = expressed.gene.counts, Percent.Expressed.Genes = percent.expressed.genes,
							gene.trimmed.percent=gene.trimmed.percent,
							
							expressed.protein.coding.genes = expressed.protein.coding.gene.counts,
							percent.expressed.protein.coding.genes = percent.protein.coding.expressed.genes,
							protein.coding.gene.trimmed.percent=gene.trimmed.protein.coding.percent)
write.table(coverage.table, file="gene_coverage_stats.txt", quote=F, row.names=F, sep="\t")
