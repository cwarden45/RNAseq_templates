deg.prefix = ""
deg.file = paste(deg.prefix,".txt",sep="")
heatmap.file  = paste("batch_corrected_",deg.prefix,".png",sep="")
metadata.table = "sample_description.txt"
expression.table = "log2_rpkm.txt"
plot.groups = c("Group")
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

deg.table = read.table(deg.file, sep="\t", head=T)
deg.genes = deg.table$symbol[deg.table$status != "No Change"]
	
fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

sample.description.table = read.table(metadata.table, sep="\t", head=T)
sample.label = sample.description.table$userID
batchID = as.character(sample.description.table[,batch.group])
	
expression.table = read.table(expression.table, sep="\t", head=T)
expression.genes = as.character(expression.table$symbol)
	
temp.rpkm = expression.table[match(deg.genes, expression.genes, nomatch=0),]
temp.rpkm = temp.rpkm[,match(sample.label, names(temp.rpkm))]
print(dim(temp.rpkm))
temp.rpkm = t(apply(temp.rpkm, 1, batch.center, batchID = batchID))
print(dim(temp.rpkm))
	
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
		if(length(deg.genes) < 25){
			colnames(std.expr) = deg.genes
		} else {
			colnames(std.expr) = rep("", length(deg.genes))
		}
		rownames(std.expr) = sample.label

		column_annotation <- as.matrix(deg.genes)
		colnames(column_annotation) <- c("")

		row_annotation <- data.frame(label1 = labelColors1, label2 = labelColors2)
		row_annotation = as.matrix(t(row_annotation))
		rownames(row_annotation) <- c(plot.groups)

		png(file = heatmap.file)
		heatmap.3(std.expr, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					RowSideColors=row_annotation, trace="none", margins = c(8,13),RowSideColorsSize=4, dendrogram="both")
		legend("topright", legend=group.levels,
							col=color.palette,
							pch=15, cex=0.7)
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
		if(length(deg.genes) < 25){
			colnames(std.expr) = deg.genes
		} else {
			colnames(std.expr) = rep("", length(deg.genes))
		}
		rownames(std.expr) = sample.label
		
		png(file = heatmap.file)
		heatmap.2(std.expr, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					 RowSideColors=labelColors, trace="none", margins = c(5,15))
		legend("topright", legend=group.levels,	col=color.palette, pch=15, cex=0.7)
		dev.off()
	}#end else
