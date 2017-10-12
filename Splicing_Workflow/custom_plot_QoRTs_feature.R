genes = c("")
geneIDs = c("")
#plot type can be "exon" or "junction"
plot.type = c("exon")
#separate multiple feature IDs per gene with commas (or leave plank if you don't have specific features to plot)
featureIDs = c("")

meta.file = "../sample_description.txt"
CPM.file = "../../Result/JunctionSeq/QoRTS_JunctionSeq_feature_CPM.txt"

fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black",colors())

meta.table = read.table(meta.file, sep="\t", head=T)
print(dim(meta.table))
indexID = meta.table$userID
groupID = meta.table$Group
#groupID = factor(groupID, levels=c("Ctrl","Trt"))

colorID = groupID

groupColors = rep("gray",times=length(colorID))
groups = levels(as.factor(as.character(colorID)))
for(i in 1:length(groups)){
	groupColors[colorID == groups[i]] = fixed.color.palatte[i]
}

CPM.table = read.table(CPM.file, sep="\t", head=T)
CPM.mat = CPM.table[,match(indexID, names(CPM.table))]

feature.info = CPM.table[,1:7]
rm(CPM.table)

for (i in 1:length(genes)){
	gene = genes[i]
	plot.type = plot.type
	target.features = featureIDs[i]
	gene.mat = CPM.mat[feature.info$geneName == gene,]
	gene.feature.info = feature.info[feature.info$geneName == gene,]
	if(nrow(gene.mat) == 0){
		print(paste("Could not find geneName '",gene,"'",sep=""))
		print("Moving on to the next gene...")
	}else{
		print(gene)
		if(plot.type[i] == "exon"){
			gene.mat = gene.mat[grep(":E\\d{3}$",gene.feature.info$featureID,perl=T),]
			gene.feature.info=gene.feature.info[grep(":E\\d{3}$",gene.feature.info$featureID,perl=T),]
		}else if(plot.type[i] == "junction"){
			gene.mat = gene.mat[grep(":J\\d{3}$",gene.feature.info$featureID,perl=T),]
			gene.feature.info=gene.feature.info[grep(":J\\d{3}$",gene.feature.info$featureID,perl=T),]
		}else{
			print(paste("Plot type '",plot.type[i],"' not recognized",sep=""))
			print("Please specify plot type as 'exon' or 'junction'")
			stop()
		}
		avg.feature.pos = (as.numeric(gene.feature.info$feature.start)+as.numeric(gene.feature.info$feature.end))/2
		
		gene.mat = gene.mat[order(avg.feature.pos),]
		gene.feature.info = gene.feature.info[order(avg.feature.pos),]
		gene.featureID = as.character(gene.feature.info$featureID)

		alt.color.features = unlist(strsplit(target.features,split=","))
		text.color = rep("black",length(gene.featureID))
		text.color[match(alt.color.features, gene.featureID)]="red"
		
		feature.min = min(gene.mat)
		feature.max = max(gene.mat)
		
		### Overall Plot ###
		pdf(paste("JunctionSeq_features_",gene,"_CPM_by_",plot.type[i],".pdf",sep=""))
		for (j in 1:ncol(gene.mat)){
			sample.CPM = as.numeric(gene.mat[,j])
			if(j == 1){
				par(mar=c(15,5,5,5))
				plot(1:length(gene.featureID), sample.CPM,
						ylim=c(feature.min, feature.max),
						xaxt="n", xlim=c(0,length(gene.featureID)+1),
						xlab="", ylab="Feature CPM", pch=16, col=groupColors[j],
						cex=0.4, main=gene)
				lines(1:length(gene.featureID), sample.CPM, col=groupColors[j])
				legend("top",groups,col=fixed.color.palatte[1:length(groups)],
						cex=0.8, pch=16, inset=-0.12, xpd=T, ncol=3)
				mtext(gene.featureID, side=1, at =1:length(gene.featureID), line=2,
					col=text.color, las=2, cex = 0.8)
			}else{
				points(1:length(gene.featureID), sample.CPM, pch=16, col=groupColors[j],
						cex=0.4)
				lines(1:length(gene.featureID), sample.CPM, col=groupColors[j])		
			}
		}#end for (j in 1:ncol(gene.mat))
		dev.off()
		
		### Dot-Plot for Specific Features ###
		for (j in 1:length(alt.color.features)){
			if(alt.color.features[j] %in% gene.feature.info$featureID){
				print(alt.color.features[j])
				
				feature.expr = as.numeric(gene.mat[gene.feature.info$featureID == alt.color.features[j],])
				
				output.file = paste("JunctionSeq_feature_",alt.color.features[j],"_Dot-Plot.pdf",sep="")
				output.file = gsub(":","_",output.file)
				pdf(output.file)
				par(mar=c(15,5,5,5))
				plot(as.numeric(groupID), feature.expr,
						xaxt="n", xlim=c(0,length(levels(groupID))+1),
						xlab="", ylab="Feature CPM", pch=16, col=groupColors,
						cex=0.8, main=alt.color.features[j])
				legend("bottom",groups,col=fixed.color.palatte[1:length(groups)],
						cex=0.8, pch=16, inset=-0.9, xpd=T, ncol=3)
				mtext(levels(groupID), side=1, at =1:length(levels(groupID)), las=2, line=2)
				dev.off()
			}else{
				print(paste("Could not find feature '",alt.color.features[j],"' for '",gene,"'",sep=""))
				print("Moving on to the next feature...")			
			}
		}#end for (j in 1:length(alt.color.features))
	}#end else

}#end for (gene in genes)