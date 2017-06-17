gtf = "path/to/TxDb_[genome]_gene.gtf"
bed = "TxDb_[genome]_gene.bed"

### edit above this line ###


extract.id = function(char, id.type){
	char.info = unlist(strsplit(char,split=";"))
	info.text = char.info[grep(id.type, char.info)]
	info.text = gsub(id.type,"",info.text)
	info.text = gsub("\"","",info.text)
	info.text = gsub(" ","",info.text)
	return(info.text)
}#end def extract.id

gtf.table = read.table(gtf, head=F, sep="\t")

geneID.values = sapply(as.character(gtf.table$V9), extract.id, id.type="gene_id")
geneID = levels(as.factor(geneID.values))

gene.chr = gtf.table$V1[match(geneID, geneID.values)]
gene.strand = gtf.table$V7[match(geneID, geneID.values)]
gene.start = tapply(gtf.table$V4, geneID.values, min)
gene.stop = tapply(gtf.table$V4, geneID.values, max)

output.table = data.frame(gene.chr, gene.start, gene.stop,
						geneID, score=rep(0,length(geneID)), gene.strand)
write.table(output.table, bed, sep="\t", quote=F, col.names=F, row.names=F)
