normalizeTotalExpression = function (geneExpr, totalReads) {
	return(geneExpr / totalReads)
}#end def normalizeTotalExpression

param.table = read.table("parameters.txt", header=T, sep="\t")
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
count.folder=as.character(param.table$Value[param.table$Parameter == "QoRTs_Merged_Folder"])
comp.name=as.character(param.table$Value[param.table$Parameter == "comp_name"])
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])

count.file = paste(user.folder,"/QoRTS_JunctionSeq_feature_counts.txt",sep="")
CPM.file = paste(user.folder,"/QoRTS_JunctionSeq_feature_CPM.txt",sep="")

stat.file = paste(user.folder,"/DSG/feature_stats_",comp.name,".txt",sep="")
stat.table = read.table(stat.file, head=T, sep="\t")

feature.info = stat.table[,1:7]
rm(stat.table)

sample.description.table = read.table(sample.description.file,head=T, sep="\t")
sampleIDs = as.character(sample.description.table$sample.ID)
sample.label = as.character(sample.description.table$userID)
total.million.aligned.reads = as.numeric(sample.description.table$aligned.reads) / 1000000
countFiles = paste(count.folder,"/",sampleIDs,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt",sep="")

for (i in 1:length(countFiles)){
	if(!file.exists(countFiles[i])){
		gzFile = paste(countFiles[i],".gz",sep="")
		
		command = paste("gunzip -c ",gzFile," > ",countFiles[i],sep="")
		system(command)
	}#end if(!file.exists(countFiles[i]))
}#end for (i in 1:length(countFiles))


#some loss of overlapping features, but mostly covered
for (i in 1:length(countFiles)){		
		coverage.table = read.table(countFiles[i],head=F, sep="\t")
		feature.name = coverage.table[,1]
		feature.count = coverage.table[,2]
		
		sample.counts = feature.count[match(feature.info[,1], feature.name)]
		if(i == 1){
			count.mat = data.frame(sampleID = sample.counts)
		}else{
			count.mat = data.frame(count.mat, sample.counts)
		}
}#end for (i in 1:length(countFiles))
colnames(count.mat) = sample.label

annotated.count.table = data.frame(feature.info, count.mat)
write.table(annotated.count.table, count.file, quote=F, sep="\t", row.names=F)

CPM = round(t(apply(count.mat, 1, normalizeTotalExpression, totalReads = total.million.aligned.reads)), digits=1)
colnames(CPM) = sample.label

annotated.CPM.table = data.frame(feature.info, CPM)
write.table(annotated.CPM.table, CPM.file, quote=F, sep="\t", row.names=F)