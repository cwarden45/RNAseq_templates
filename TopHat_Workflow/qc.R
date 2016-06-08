colLab <- function(n, labelColors, clusMember) { 
   if(is.leaf(n)) { 
       a <- attributes(n) 
	   #print(a)
       # clusMember - vector of sample names (ordered to match label color.palette)
       # labelColors - a vector of color.palette for the above grouping 
       labCol <- labelColors[clusMember == a$label]
	   #print(labCol)
       attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol) 
   } 
   n 
}

count.defined.values <- function(arr)
{
	return(length(arr[!is.na(arr)]))
}#end def count.values

param.table = read.table("parameters.txt", header=T, sep="\t")
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
rpkm.file = as.character(param.table$Value[param.table$Parameter == "rpkm_file"])
min.expression = as.numeric(as.character(param.table$Value[param.table$Parameter == "rpkm_expression_cutoff"]))
plot.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "plot_groups"]), split=","))

fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black",colors())

sample.table = read.table(sample.description.file, sep="\t", header=T)
userID = as.character(sample.table$userID)

normalized.table = read.table(rpkm.file, sep="\t", header=T)
normalized.mat = normalized.table[,match(userID, names(normalized.table))]
normalized.mat[normalized.mat < log2(min.expression)] = NA
expr.max <- ceiling(max(normalized.mat, na.rm = T))
expr.min <- ceiling(min(normalized.mat, na.rm = T))

pca.values <- prcomp(na.omit(data.matrix(normalized.mat)))
print(dim(pca.values))
pc.values <- data.frame(pca.values$rotation)
variance.explained <- (pca.values $sdev)^2 / sum(pca.values $sdev^2)
pca.table <- data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))

pca.text.file = "pca_values.txt"
write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")

for (group in plot.groups){
	print(group)

	groups = levels(as.factor(sample.table[,group]))
	num.sample.types = length(groups)
	color.palette <- fixed.color.palatte[1:length(groups)]

	labelColors = rep("black",times=ncol(normalized.mat))
	for (i in 1:length(groups)){
		labelColors[sample.table[,group] == as.character(groups[i])] = color.palette[i]
	}#end for (i in 1:length(groups))

	pca.file = paste("pca_by_",group,".png",sep="")
	png(file=pca.file)
	plot(pc.values$PC1, pc.values$PC2, col = labelColors, xlab = paste("PC1 (",round(100* variance.explained[1] , digits = 2),"%)", sep = ""),
			ylab = paste("PC2 (",round(100* variance.explained[2] , digits = 2),"%)", sep = ""), pch=19)
	legend("topright",legend=groups,col=color.palette,  pch=19)
	dev.off()

	box.file = paste("box_plot_by_",group,".png",sep="")
	png(file=box.file)
	boxplot(normalized.mat, col=labelColors, xaxt='n')
	legend("topright",legend=groups,col=color.palette,  pch=19)
	dev.off()

	cluster.file = paste("cluster_by_",group,".png",sep="")
	dist1 <- dist(as.matrix(t(normalized.mat)))
	clusMember <- groups
	hc <- hclust(dist1)
	dend1 <- as.dendrogram(hc)
	png(file = cluster.file)
	dend1 <- dendrapply(dend1, colLab, labelColors=labelColors, clusMember=userID) 
	a <- attributes(dend1) 
	attr(dend1, "nodePar") <- c(a$nodePar, lab.col = labelColors) 
	op <- par(mar = par("mar") + c(0,0,0,10)) 
	plot(dend1, horiz=T)
	par(op) 
	dev.off()

	hist.file = paste("density_by_",group,".png",sep="")
	png(file = hist.file)
	for (i in 1:ncol(normalized.mat))
		{		
			data <- as.numeric(t(normalized.mat[,i]))
			
			if(i == 1)
				{
					den <- density(data, na.rm=T,from=expr.min, to=expr.max)
					expr <- den$x
					freq <- den$y
					plot(expr, freq, type="l", xlab = paste("Log2(RPKM > ",min.expression,") Expression",sep=""), ylab = "Density",
							xlim=c(expr.min,expr.max), ylim=c(0,0.2), col=labelColors[i])
					legend("topright",legend=groups,col=color.palette,  pch=19)
				}#end if(i == 1)
			else
				{
					den <- density(data, na.rm=T,from=expr.min, to=expr.max)
					expr <- den$x
					freq <- den$y
					lines(expr, freq, type = "l", col=labelColors[i])
				}#end else
		}#end for (i in 1:length(ncol(normalized.mat)))
	dev.off()	
}#end for (group in plot.groups)