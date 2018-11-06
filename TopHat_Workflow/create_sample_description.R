forward.reads = list.files("../Reads")

forward.reads = forward.reads[grep("R1_001.fastq",forward.reads)]
forward.reads = gsub(".gz$","",forward.reads)

reg.obj = gregexpr("_(L\\d{3})_",forward.reads)
laneID = unlist(regmatches(forward.reads, reg.obj))
laneID = gsub("_","",laneID)

sampleID = gsub("_S\\d+_L\\d{3}_R1_001.fastq","",forward.reads)

extract.sampleID = function(char){
	sample.arr = unlist(strsplit(char,split="_"))
	sample.arr = sample.arr[2:length(sample.arr)]
	newID = paste(sample.arr,collapse=".")
	return(newID)
}

reg.obj = gregexpr("^\\d+",sampleID)
seqID = unlist(regmatches(sampleID, reg.obj))

print(sampleID)

userID=as.character(unlist(sapply(sampleID, extract.sampleID)))

Group = rep(NA, length(userID))
Batch = rep(NA, length(userID))
HTseq.file = paste(sampleID,"_gene_counts.txt",sep="")

output.table = data.frame(sampleID, seqID, userID,
							Group, HTseq.file)
write.table(output.table,"sample_description.txt", row.names=F, quote=F, sep="\t")