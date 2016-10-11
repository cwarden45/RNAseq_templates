param.table = read.table("parameters.txt", header=T, sep="\t")
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
count.folder=as.character(param.table$Value[param.table$Parameter == "QoRTs_Merged_Folder"])
plot.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "plot_groups"]), split=","))


fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black",colors())
cov.max = log2(35000)

sample.description.table = read.table(sample.description.file,head=T, sep="\t")
sampleIDs = as.character(sample.description.table$sample.ID)
countFiles = paste(count.folder,"/",sampleIDs,"/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt",sep="")

for (i in 1:length(countFiles)){
	if(!file.exists(countFiles[i])){
		gzFile = paste(countFiles[i],".gz",sep="")
		
		command = paste("gunzip -c ",gzFile," > ",countFiles[i],sep="")
		system(command)
	}#end if(!file.exists(countFiles[i]))
}#end for (i in 1:length(countFiles))

for (group in plot.groups){
	print(group)

	qc.grp = sample.description.table[,group]
	qc.name = sampleIDs[!is.na(qc.grp)]
	qc.grp = qc.grp[!is.na(qc.grp)]
	
	splice.files = countFiles[match(qc.name, sampleIDs)]
	groups = levels(as.factor(as.character(qc.grp)))
	color.palette <- fixed.color.palatte[1:length(groups)]
	labelColors = rep("black",times=length(qc.name))
	for (i in 1:length(groups)){
		labelColors[qc.grp == as.character(groups[i])] = color.palette[i]
	}#end for (i in 1:length(groups))
	
	hist.file = paste("junction_dist_by_",group,".png",sep="")
	png(file = hist.file)
	for (i in 1:length(splice.files))
		{		
			coverage.table = read.table(splice.files[i],head=F, sep="\t")
			feature.name = coverage.table[,1]
			feature.count = coverage.table[,2]
			junction.names = feature.name[grep(":J\\d{3}$",feature.name,perl=T)]
			junction.counts = log2(feature.count[grep(":J\\d{3}$",feature.name,perl=T)]+1)
			
			if(i == 1)
				{
					den <- density(junction.counts, na.rm=T, from=0, to=cov.max)
					expr <- den$x
					freq <- den$y
					plot(expr, freq, type="l", xlab = "log2(Coverage+1)", ylab = "Density",
							xlim=c(0,cov.max), ylim=c(0,0.25), col=labelColors[i])
					legend("topright",legend=groups,col=color.palette,  pch=19)
				}#end if(i == 1)
			else
				{
					den <- density(junction.counts, na.rm=T, from=0, to=cov.max)
					expr <- den$x
					freq <- den$y
					lines(expr, freq, type = "l", col=labelColors[i])
				}#end else
		}#end for (i in 1:length(ncol(normalized.mat)))
	dev.off()
}#end for (group in plot.groups)