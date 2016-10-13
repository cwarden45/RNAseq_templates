avgGroupExpression <- function (geneExpr, groups) {
	avg.expr = tapply(geneExpr, groups, mean)
	return(avg.expr)
}#end def avgGroupExpression

ratio2fc <- function(value)
{
	if(value >= 0){
		return(2^value)
	} else {
		return (-2^(-value))
	}
}#end def ratio2fc

param.table = read.table("parameters.txt", header=T, sep="\t")
comp.name=as.character(param.table$Value[param.table$Parameter == "comp_name"])
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
deg.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "deg_groups"]), split=","))
count.folder=as.character(param.table$Value[param.table$Parameter == "QoRTs_Merged_Folder"])
goseq.flag = as.character(param.table$Value[param.table$Parameter == "run_goseq"])
trt = as.character(param.table$Value[param.table$Parameter == "treatment_group"])
fc.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "fold_change_cutoff"]))
pvalue.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "pvalue_cutoff"]))
fdr.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "fdr_cutoff"]))
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
sample.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
gene.mapping.file = as.character(param.table$Value[param.table$Parameter == "GENCODE_Gene_Info"])

raw.output.folder = paste(count.folder,"/",comp.name,"/",sep="")

result.file = paste(raw.output.folder,"allGenes.results.txt",sep="")
counts.file = paste(raw.output.folder,"allGenes.expression.data.txt",sep="")
gene.file = paste(raw.output.folder,"sigGenes.genewiseResults.txt",sep="")
custom.result = paste(user.folder,"/DSG/feature_stats_",comp.name,".txt",sep="")
gene.result = paste(user.folder,"/DSG/sig_gene_stats_",comp.name,".txt",sep="")

result.table = read.table(result.file, head=T, sep="\t")
featureID = result.table$featureID

mapping.table = read.table(gene.mapping.file, head=T, sep="\t")
symbol = mapping.table$Gene.Name[match(result.table$transcripts, mapping.table$Ensembl.TranscriptID)]

temp.exonID = featureID[grep(":E\\d{3}$",featureID,perl=T)]
temp.exon.pvalue = result.table$pvalue[match(temp.exonID, featureID)]
temp.exon.fdr = result.table$padjust[match(temp.exonID, featureID)]

sample.table = read.table(sample.file, sep="\t", header=T)
sampleID = as.character(sample.table$userID)
countID = gsub("-",".",sampleID)
countID = paste("normCount_",countID,sep="")
group = as.factor(sample.table[,deg.groups][,1])
group.levels =levels(group)

count.table = read.table(counts.file, sep="\t", header=T)
count.feature = count.table$featureID
count.mat = count.table[,match(countID, names(count.table))]

group.counts = data.frame(t(apply(count.mat, 1, avgGroupExpression, groups = group)))
group.counts = log2(group.counts + 1)
colnames(group.counts) = paste("avg.log2.normCounts",levels(group),sep=".")


dsg.group = group
if (length(group.levels) > 2){
	dsg.group = as.character(dsg.group)
	dsg.group[dsg.group != trt] = "cntl"
	dsg.group = as.factor(as.character(dsg.group))
}
dsg.counts = data.frame(t(apply(count.mat, 1, avgGroupExpression, groups = dsg.group)))
colnames(dsg.counts) = levels(dsg.group)
dsg.counts = log2(dsg.counts + 1)
trt.expr = dsg.counts[,trt]
cntl.expr = dsg.counts[,levels(dsg.group)[levels(dsg.group) != trt]]
log2ratio = round(trt.expr - cntl.expr, digits = 2)
fc = round(sapply(log2ratio, ratio2fc), digits = 2)

custom.table = data.frame(featureID = featureID, geneName = symbol, transcripts = result.table$transcripts,
				feature.chr = result.table$chr, feature.start = result.table$start, feature.end = result.table$end,
				group.counts, log2ratio = log2ratio, fold.change=fc,
				feature.pvalue =result.table$pvalue, feature.fdr = result.table$padjust,
				gene.fdr = result.table$geneWisePadj)
write.table(custom.table, custom.result, sep="\t", row.names=F)

temp.exon.fc = fc[match(temp.exonID,count.feature)]

exonID = temp.exonID[!is.na(temp.exon.fdr) & !is.na(temp.exon.fc)]
exon.fdr = temp.exon.fdr[!is.na(temp.exon.fdr) & !is.na(temp.exon.fc)]
exon.pvalue = temp.exon.pvalue[!is.na(temp.exon.fdr) & !is.na(temp.exon.fc)]
exon.fc = temp.exon.fc[!is.na(temp.exon.fdr) & !is.na(temp.exon.fc)]

up.exons = as.character(exonID[(exon.fc > fc.cutoff) & (exon.pvalue < pvalue.cutoff)& (exon.fdr < fdr.cutoff)])
print(paste("Up-Regulated Exons:",length(up.exons),sep=""))
down.exons = as.character(exonID[(exon.fc < - fc.cutoff) & (exon.pvalue < pvalue.cutoff)& (exon.fdr < fdr.cutoff)])
print(paste("Down-Regulated Exons:",length(down.exons),sep=""))

gene.table = read.table(gene.file, head=T, sep="\t")

dsg = as.character(gene.table$geneID[gene.table$geneWisePadj < fdr.cutoff])
dsg = as.character(mapping.table$Gene.Name[match(dsg,as.character(mapping.table$Ensembl.GeneID))])
print(paste("Altered Genes:",length(dsg),sep=""))

gene.table = gene.table[gene.table$geneWisePadj < fdr.cutoff,]
write.table(gene.table, gene.result, sep="\t", row.names=F)

if(goseq.flag == "yes"){
	library(goseq)

	gene.symbol = as.character(levels(as.factor(as.character(symbol))))
	
	deg = rep(1, times = length(gene.symbol))
	deg[is.na(match(gene.symbol,dsg))]=0
	
	names(deg)=gene.symbol
	
	bias.file = paste(comp.name,"_goseq_gene_splicing_bias.png",sep="")
	png(bias.file)
	pwf=nullp(deg,genome,"geneSymbol")
	GO.wall=goseq(pwf,genome,"geneSymbol")
	GO.wall = data.frame(GO.wall, overenrichment.fdr=p.adjust(GO.wall$over_represented_pvalue))

	GO.wall = GO.wall[GO.wall$over_represented_pvalue < 0.05,]
	genes2go=getgo(gene.symbol[deg == 1],genome,"geneSymbol")
	go2genes=goseq:::reversemapping(genes2go)
	go.genes = rep("", times=nrow(GO.wall))
	dev.off()

	for( i in 1:length(go.genes)){
		go.category = GO.wall$category[i]
		result = paste(go2genes[[go.category]], sep="", collapse=",")
		go.genes[i] = result
	}#end for( i in 1:length(go.genes))

	go.table = data.frame(GO.wall, genes=go.genes)
	go.file = paste(user.folder,"/GO/",comp.name,"_goseq_gene_splicing.txt",sep="")
	write.table(go.table, file=go.file, row.names=F, quote=F, sep="\t")
}#end if(goseq.flag)
