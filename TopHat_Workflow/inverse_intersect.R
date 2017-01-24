#input info
normal.cat = "comp1"
normal.dir.file = ".txt"
normal.upID = "Increased Expression"
normal.downID = "Decreased Expression"
switched.cat = "comp2"
switched.dir.file = ".txt"
switched.upID = "Increased Expression"
switched.downID = "Decreased Expression"

#if you want to make a heatmap
heatmap.flag = FALSE
metadata.table = "sample_description.txt"
expression.table = "log2_rpkm.txt"
plot.groups = c("Group")

#if you want to run goseq
genome = "hg19"
goseq.flag = FALSE

#output names
overall.upID = "Increased comp1 Expression"
overall.downID = "Decreased comp1 Expression"
intersect.output.table = ".txt"
intersect.output.heatmap = ".png"
goseq.up.file = ".txt"
goseq.down.file = ".txt"

normal.dir.table = read.table(normal.dir.file, sep="\t", head=T)
switched.dir.table = read.table(switched.dir.file, sep="\t", head=T)

matched.genes = normal.dir.table$symbol[match(switched.dir.table$symbol, normal.dir.table$symbol, nomatch=0)]
normal.dir.table = normal.dir.table[match(matched.genes, normal.dir.table$symbol, nomatch=0), ]
switched.dir.table = switched.dir.table[match(matched.genes, switched.dir.table$symbol, nomatch=0), ]

status = rep("No Change", times=length(matched.genes))
up.genes = matched.genes[(normal.dir.table$status == normal.upID) & (switched.dir.table$status == switched.downID)]
status[(normal.dir.table$status == normal.upID) & (switched.dir.table$status == switched.downID)]=overall.upID
print(paste("Up-Regulated: ",length(up.genes),sep=""))
down.genes = matched.genes[(normal.dir.table$status == normal.downID) & (switched.dir.table$status == switched.upID)]
status[(normal.dir.table$status == normal.downID) & (switched.dir.table$status == switched.upID)]=overall.downID
print(paste("Down-Regulated: ",length(down.genes),sep=""))

gene.info = normal.dir.table[,1:2]
normal.stats = normal.dir.table[3:ncol(normal.dir.table)]
colnames(normal.stats) = paste(normal.cat, colnames(normal.stats),sep=".")
switched.stats = switched.dir.table[3:ncol(switched.dir.table)]
colnames(switched.stats) = paste(switched.cat, colnames(switched.stats),sep=".")

output.table = data.frame(gene.info, normal.stats, switched.stats, overall.status=status)
write.table(output.table, intersect.output.table, sep="\t", quote=F, row.names=F)

#heatmap
if(heatmap.flag){
	library(gplots)
	
	standardize.arr <- function(arr)
	{
		center.arr = as.numeric(arr) - mean(as.numeric(arr), na.rm=T)
		norm.arr = center.arr / sd(center.arr, na.rm=T)
		return(norm.arr)
	}#end def standardize.arr
	
	fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black",colors())
	sample.description.table = read.table(metadata.table, sep="\t", head=T)
	sample.label = sample.description.table$userID
	
	expression.table = read.table(expression.table, sep="\t", head=T)
	expression.genes = as.character(expression.table$symbol)
	
	deg.genes = as.character(matched.genes[status != "No Change"])
	heatmap.deg.table = expression.table[match(deg.genes, expression.genes, nomatch=0),]
	temp.rpkm = heatmap.deg.table[,match(sample.label, names(heatmap.deg.table))]
	
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

		png(file = intersect.output.heatmap)
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
		
		png(file = intersect.output.heatmap)
		heatmap.2(std.expr, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					 RowSideColors=labelColors, trace="none", margins = c(5,15))
		dev.off()
	}#end else
}#end heatmap code

#goseq
if (goseq.flag){
	library(goseq)
	
	deg = as.integer(status == overall.upID)
	deg = tapply(deg, as.character(matched.genes), sum)
	deg[deg >= 1] = 1
	gene.symbol = as.character(levels(as.factor(as.character(matched.genes))))
	names(deg)=gene.symbol
	
	bias.file = "goseq_up_bias.png"
	png(bias.file)
	pwf=nullp(deg,genome,"geneSymbol")
	GO.wall=goseq(pwf,genome,"geneSymbol")
	GO.wall = data.frame(GO.wall, overenrichment.fdr=p.adjust(GO.wall$over_represented_pvalue, "fdr"))

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
	write.table(go.table, file=goseq.up.file, row.names=F, quote=F, sep="\t")

	deg = as.integer(status == overall.downID)
	deg = tapply(deg, as.character(matched.genes), sum)
	deg[deg >= 1] = 1
	gene.symbol = as.character(levels(as.factor(as.character(matched.genes))))
	names(deg)=gene.symbol

	bias.file = "goseq_down_bias.png"
	png(bias.file)
	pwf=nullp(deg,genome,"geneSymbol")
	GO.wall=goseq(pwf,genome,"geneSymbol")
	GO.wall = data.frame(GO.wall, overenrichment.fdr=p.adjust(GO.wall$over_represented_pvalue, "fdr"))

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
	write.table(go.table, file=goseq.down.file, row.names=F, quote=F, sep="\t")
}#end goseq code
