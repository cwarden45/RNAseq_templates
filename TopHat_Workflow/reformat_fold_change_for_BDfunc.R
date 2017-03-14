deg.file = ".txt"
out.file = "_fold_change.txt"
use_log2ratio = FALSE

deg.table = read.table(deg.file,head=T,sep="\t")

if(use_log2ratio){
	output.table = data.frame(gene=deg.table$symbol, fc=deg.table$log2ratio)
	write.table(output.table, out.file, sep="\t", row.names=F, quote=F)
}else{
	output.table = data.frame(gene=deg.table$symbol, fc=deg.table$fold.change)
	write.table(output.table, out.file, sep="\t", row.names=F, quote=F)
}

#manually replace header...
#BD-Func