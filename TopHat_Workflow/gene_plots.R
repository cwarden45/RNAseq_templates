expression.file = "log2_fpkm.txt"
meta.file = "sample_description.txt"
genes = c("ACTB","GAPDH")

#you'll probably still want to modify mtext parameters...otherwise, most modifications can be done above

fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black",colors())

meta.table = read.table(meta.file, sep="\t", head=T)
indexID = meta.table$userID
groupID = meta.table$Group
#groupID = factor(groupID, levels=c("Trt","Cntl"))

colorID = groupID

groupColors = rep("gray",times=length(colorID))
groups = levels(as.factor(as.character(colorID)))
for(i in 1:length(groups)){
	groupColors[colorID == groups[i]] = fixed.color.palatte[i]
}

expression.table = read.table(expression.file, sep="\t", head=T)
expression.mat = expression.table[,match(indexID, names(expression.table))]
expr.min = min(expression.mat, na.rm=T)
expr.max = max(expression.mat, na.rm=T)
all.genes = as.character(expression.table$symbol)

for (gene in genes){
	print(gene)
	gene.expr = as.numeric(expression.mat[all.genes == gene,])

	output.file = paste("aligned_FPKM_",gene,".png",sep="")
	png(output.file)
	par(mar=c(15,5,5,5))
	plot(as.numeric(groupID), gene.expr,
			ylim=c(expr.min, expr.max),
			xaxt="n", xlim=c(0,length(levels(groupID))+1),
			xlab="", ylab="log2(FPKM + 0.1) Expression", pch=16, col=groupColors,
			cex=1.2, main=gene)
	legend("bottom",groups,col=fixed.color.palatte[1:length(groups)],
			cex=0.8, pch=16, inset=-0.9, xpd=T, ncol=3)
	mtext(levels(groupID), side=1, at =1:length(levels(groupID)), las=2, line=2)
	dev.off()
	
	output.file = paste("aligned_FPKM_",gene,"_zoom.png",sep="")
	png(output.file)
	par(mar=c(15,5,5,5))
	plot(as.numeric(groupID), gene.expr,
			xaxt="n", xlim=c(0,length(levels(groupID))+1),
			xlab="", ylab="log2(FPKM + 0.1) Expression", pch=16, col=groupColors,
			cex=1.2, main=gene)
	legend("bottom",groups,col=fixed.color.palatte[1:length(groups)],
			cex=0.8, pch=16, inset=-0.9, xpd=T, ncol=3)
	mtext(levels(groupID), side=1, at =1:length(levels(groupID)), las=2, line=2)
	dev.off()
	
	output.file = paste("aligned_FPKM_",gene,"_zoom.pdf",sep="")
	pdf(output.file)
	par(mar=c(15,5,5,5))
	plot(as.numeric(groupID), gene.expr,
			xaxt="n", xlim=c(0,length(levels(groupID))+1),
			xlab="", ylab="log2(FPKM + 0.1) Expression", pch=16, col=groupColors,
			cex=1.2, main=gene)
	legend("bottom",groups,col=fixed.color.palatte[1:length(groups)],
			cex=0.8, pch=16, inset=-0.9, xpd=T, ncol=3)
	mtext(levels(groupID), side=1, at =1:length(levels(groupID)), las=2, line=2)
	dev.off()
}#end for (gene in genes)
