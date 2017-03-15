in.file = "log2_fpkm.txt"
gct.file = "../Result/GSEA/log2_fpkm.gct"
cls.file = "../Result/GSEA/log2_fpkm.cls"
meta.file = "sample_description.txt"
group = "Group"
convert.upper=TRUE

meta.table = read.table(meta.file,head=T,sep="\t")
sampleID = meta.table$userID
cls.group = meta.table[,group]

expr.table = read.table(in.file,head=T,sep="\t")
fpkm.mat = round(expr.table[,match(sampleID, names(expr.table))], digits=2)

if(convert.upper){
	#if working with other organisms (like mouse), convert symbol to upper-case
	gene = toupper(expr.table$symbol)
}#end if(convert.upper)

output.table = data.frame(NAME=gene,Description=gene,
						fpkm.mat)
write.table(output.table, gct.file, sep="\t", row.names=F, quote=F)

print(".gct file: manually change header...")
print("#1.2")
print(paste(nrow(fpkm.mat),ncol(fpkm.mat),sep=" "))

print(".cls file: manually change header...")
print(paste("",length(cls.group)," ",length(levels(cls.group))," 1",sep=""))
print(paste("# ",paste(levels(cls.group),collapse = " "),sep=""))
output.table = data.frame(t(data.frame(cls.group)))
write.table(output.table, cls.file, sep="\t", quote=F,
			row.names=F, col.names=F)
print("you might have to modify line #2 based upon order samples appear in table")
print("or, if continuous, modify as described in http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html?_Phenotype_Labels")
