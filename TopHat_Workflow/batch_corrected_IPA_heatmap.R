IPA.prefix = ""
IPA.file = paste(IPA.prefix,".txt",sep="")
heatmap.file  = paste(IPA.prefix,".png",sep="")
metadata.table = "sample_description.txt"
expression.table = "log2_fpkm.txt"
cluster.distance="Pearson_Dissimilarity"
plot.groups = c("Group","Batch")
batch.group = "Batch"

###edit code above this line###

library(gplots)

standardize.arr = function(arr)
{
	center.arr = as.numeric(arr) - mean(as.numeric(arr), na.rm=T)
	norm.arr = center.arr / sd(center.arr, na.rm=T)
	return(norm.arr)
}#end def standardize.arr

batch.center = function(arr, batchID)
{
	batch.means = tapply(as.numeric(arr),as.character(batchID), mean, na.rm=T)
	for(batch in names(batch.means)){
		arr[batchID == batch] = arr[batchID == batch] - batch.means[batch]
	}#end for(batch in names(batch.means))
	return(arr)
}#end def batch.center

IPA.table = read.table(IPA.file, sep="\t", head=T)
IPA.genes = as.character(IPA.table$Genes.in.dataset)
IPA.annotation = as.character(IPA.table$Findings)
IPA.annotation[grep("Affected",IPA.annotation)]="Affected"
IPA.annotation[grep("Increases",IPA.annotation)]="Activated"
IPA.annotation[grep("Decreases",IPA.annotation)]="Inhibited"
evidence.color = rep("gray", length(IPA.annotation))
evidence.color[IPA.annotation=="Activated"]="red"
evidence.color[IPA.annotation=="Inhibited"]="blue"

fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

sample.description.table = read.table(metadata.table, sep="\t", head=T)
sample.label = sample.description.table$userID
batchID = as.character(sample.description.table[,batch.group])
	
expression.table = read.table(expression.table, sep="\t", head=T)
expression.genes = as.character(expression.table$symbol)
temp.rpkm = expression.table[,match(sample.label, names(expression.table))]
rownames(temp.rpkm)=expression.genes
	
matched.genes = expression.genes[match(IPA.genes, expression.genes, nomatch=0)]
temp.rpkm = temp.rpkm[match(matched.genes, expression.genes),]
evidence.color = evidence.color[match(matched.genes, IPA.genes)]
IPA.genes = IPA.genes[match(matched.genes, IPA.genes)]
print(dim(temp.rpkm))
temp.rpkm = t(apply(temp.rpkm, 1, batch.center, batchID = batchID))
print(dim(temp.rpkm))

	cor.dist = function(mat){
		cor.mat = cor(as.matrix(t(mat)))
		dis.mat = 1 - cor.mat
		return(as.dist(dis.mat))	
	}#end def cor.dist

	if (cluster.distance == "Pearson_Dissimilarity"){
		print("Using Pearson Dissimilarity as Distance in Heatmap...")
		dist.fun = cor.dist
	}else{
		dist.fun=dist
	}

	if(length(plot.groups) > 1){
		source("heatmap.3.R")
		grp1 = as.character(sample.description.table[,plot.groups[1]])
		grp2 = as.character(sample.description.table[,plot.groups[2]])
		group.levels = c(levels(as.factor(grp1)),levels(as.factor(grp2)))

		color.palette <- fixed.color.palatte[1:length(group.levels)]
		labelColors1 = rep("black",times=length(sample.label))
		for (i in 1:length(group.levels)){
			labelColors1[grp1 == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(group.levels))
		labelColors2 = rep("black",times=length(sample.label))
		for (i in 1:length(group.levels)){
			labelColors2[grp2 == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(group.levels))
		
		std.expr = apply(temp.rpkm, 1, standardize.arr)
		if(length(IPA.genes) < 25){
			colnames(std.expr) = IPA.genes
		} else {
			colnames(std.expr) = rep("", length(IPA.genes))
		}
		rownames(std.expr) = sample.label

		#have genes as rows
		std.expr = t(std.expr)
		
		column_annotation = data.frame(label2 = labelColors2, label1 = labelColors1)
		column_annotation = as.matrix(column_annotation)
		colnames(column_annotation) <- c(rev(plot.groups))

		row_annotation <- data.frame(evidence.color)
		row_annotation = as.matrix(t(row_annotation))
		rownames(row_annotation) <- c("IPA Prediction")

		png(file = heatmap.file)
		heatmap.3(std.expr,   distfun = dist.fun, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					RowSideColors=row_annotation, ColSideColors=column_annotation,
					RowSideColorsSize=2, ColSideColorsSize=2,
					trace="none", margins = c(10,13), dendrogram="both")
		legend("topright", legend=group.levels,
							col=color.palette,
							pch=15, cex=0.7)
		legend("bottomleft", c("Activated","Affected","Inhibited"), col=c("red","gray","blue"),
				pch=15, cex=0.7, xpd=T, inset=-0.1)
		dev.off()
	} else {
		grp = sample.description.table[,plot.groups]
		group.levels = levels(as.factor(grp))

		color.palette <- fixed.color.palatte[1:length(group.levels)]
		labelColors = rep("black",times=length(sample.label))
		for (i in 1:length(group.levels)){
			labelColors[grp == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(group.levels))

		std.expr = apply(temp.rpkm, 1, standardize.arr)
		if(length(IPA.genes) < 50){
			colnames(std.expr) = IPA.genes
		} else {
			colnames(std.expr) = rep("", length(IPA.genes))
		}
		rownames(std.expr) = sample.label
		
		#have genes as rows
		std.expr = t(std.expr)
		
		png(file = heatmap.file)
		heatmap.2(std.expr,   distfun = dist.fun, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					 ColSideColors=labelColors, RowSideColors=evidence.color,
					 trace="none", margins = c(10,15))
		legend("topright", legend=group.levels,	col=color.palette, pch=15, cex=0.7)
		legend("bottomleft", c("Activated","Affected","Inhibited"), col=c("red","gray","blue"), pch=15, cex=0.7)
		dev.off()
	}#end else
