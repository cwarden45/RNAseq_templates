new.FDR.cutoff = 0.25

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

raw.output.folder = paste(count.folder,"/",comp.name,"/",sep="")
plot.folder = paste(user.folder,"/DSG/",comp.name,"/",sep="")

writeCompleteResults(jscs, outfile.prefix=raw.output.folder, save.jscs = FALSE,
						FDR.threshold = new.FDR.cutoff, gzip.output = FALSE)

if(FALSE){
	#you most likely don't want to output all genes in a more liberally defined list
	#...but you can change to TRUE if original criteria was too string
	buildAllPlots(jscs=jscs, outfile.prefix = plot.folder, sequencing.type = lib.type,
				use.plotting.device = "png",FDR.threshold = new.FDR.cutoff)
}#end if(FALSE)
