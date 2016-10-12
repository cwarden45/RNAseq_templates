#make sure you use ENCODE Gene ID, not normal gene symbol
genes = c()
geneIDs = c()

param.table = read.table("parameters.txt", header=T, sep="\t")
comp.name=as.character(param.table$Value[param.table$Parameter == "comp_name"])
count.folder=as.character(param.table$Value[param.table$Parameter == "QoRTs_Merged_Folder"])
lib.type=as.character(param.table$Value[param.table$Parameter == "pairing"])
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])

plot.folder = paste(user.folder,"/DSG/",comp.name,"/",sep="")

if (lib.type == "SE"){
	lib.type = "single-end"
}else if (lib.type == "PE"){
	lib.type = "paired-end"
}else{
	stop(paste("Need to map `pairing` to `single-end` or `paired-end`:",lib.type))
}

library("JunctionSeq")

jscsImage = paste(count.folder,"/",comp.name,"/jscs.RData",sep="")
load(jscsImage)

for (i in 1:length(genes)){
	gene = genes[i]
	geneID = geneIDs[i]
	print(gene)
	#buildAllPlotsForGene(jscs=jscs,geneID=geneID,sequencing.type = lib.type)
	output.name = paste(plot.folder,"/",gene,"_normCounts_TX.png",sep="")
	png(output.name)
	plotJunctionSeqResultsForGene(geneID=geneID, jscs=jscs,sequencing.type = lib.type, plot.type="normCounts", displayTranscripts=TRUE)
	dev.off()
}#end for (i in 1:length(genes))