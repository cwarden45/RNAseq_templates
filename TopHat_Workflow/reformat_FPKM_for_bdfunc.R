in.file = "log2_fpkm.txt"
out.file = "log2_fpkm_GROUP.txt"
meta.file = "sample_description.txt"
groupID = "Group"

meta.table = read.table(meta.file,head=T,sep="\t")
sampleID = meta.table$userID
group = meta.table[,groupID]

expr.table = read.table(in.file,head=T,sep="\t")
rpkm.mat = expr.table[,match(sampleID, names(expr.table))]
colnames(rpkm.mat)=group

output.table = data.frame(gene=expr.table$symbol, rpkm.mat)
colnames(output.table) = c("gene",as.character(group))
write.table(output.table, out.file, sep="\t", row.names=F, quote=F)

#manually change header...
#BD-Func
#groupA	groupB groupA groupB groupA groupB
