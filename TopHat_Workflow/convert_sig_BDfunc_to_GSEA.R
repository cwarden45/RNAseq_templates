bdfunc.sig = ".txt"
gsea.sig = ".gmt"
convert.upper=FALSE

bdfunc.table = read.table(bdfunc.sig, head=T,sep="\t")

gsea.lines = c()
for (i in 1:nrow(bdfunc.table)){
	bd.sig.name = bdfunc.table$Pathway[i]
	
	#up-regulated genes
	print("Activated Genes")
	up.genes = unlist(strsplit(as.character(bdfunc.table$Positive.Genes[i]), split=","))
	if(convert.upper)
		up.genes = toupper(up.genes)
	print(length(up.genes))
	new.list = c(paste(bd.sig.name,".UP",sep=""),"custom_signature",up.genes)
	gsea.txt = paste(new.list, collapse="\t")
	gsea.lines = c(gsea.lines, gsea.txt)
	
	#down-regulated genes
	print("Inhibited Genes")
	down.genes = unlist(strsplit(as.character(bdfunc.table$Negative.Genes[i]), split=","))
	if(convert.upper)
		down.genes = toupper(down.genes)
	print(length(down.genes))
	new.list = c(paste(bd.sig.name,".DOWN",sep=""),"custom_signature",down.genes)
	gsea.txt = paste(new.list, collapse="\t")
	gsea.lines = c(gsea.lines, gsea.txt)
}#end for (i in 1:nrow(bdfunc.table))

gsea.table = data.frame(gsea.lines)
write.table(gsea.table, gsea.sig, quote=F, row.names=F, col.names=F)