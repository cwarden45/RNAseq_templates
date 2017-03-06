deg.file = ""
bdfunc.enrichment.file = "_custom_sig.txt"
sig.name = ""
up.status = "[trt] Up"
down.status = "[trt] Down"

deg.table = read.table(deg.file,head=T,sep="\t")

up.genes = paste(as.character(deg.table$symbol[deg.table$status == up.status]),collapse=",")
down.genes = paste(as.character(deg.table$symbol[deg.table$status == down.status]),collapse=",")

bd.func.table = data.frame(Pathway=sig.name, Positive.Genes=up.genes, Negative.Genes=down.genes)
write.table(bd.func.table, bdfunc.enrichment.file, quote=F, sep="\t", row.names=F)
