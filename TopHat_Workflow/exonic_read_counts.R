library(GenomicAlignments)

param.table = read.table("parameters.txt", header=T, sep="\t")
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
alignment.folder = as.character(param.table$Value[param.table$Parameter == "Alignment_Folder_PC"])
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
htseq.anno.folder = as.character(param.table$Value[param.table$Parameter == "HTseq_input_folder"])
total.reads.file = as.character(param.table$Value[param.table$Parameter == "total_counts_file"])
aligned.stats.file = as.character(param.table$Value[param.table$Parameter == "aligned_stats_file"])

setwd(output.folder)
length.file = paste(htseq.anno.folder,"\\",genome,"_chr_length.txt",sep="")
full.annotation.file = paste(htseq.anno.folder,"\\TxDb_",genome,"_exon_annotations.txt",sep="")

length.table = read.table(length.file, header=T, sep="\t")
chr_length = as.numeric(length.table$Length)
names(chr_length) = as.character(length.table$Chr)

total.reads.table = read.table(total.reads.file, header=T, sep="\t")
sampleIDs = as.character(total.reads.table$Sample)

exon.info = read.table(full.annotation.file, header=T, sep="\t")
    
#remove non-canonical chromosomes
nonCanonical <- grep("_", exon.info$chr)
 if (length(nonCanonical) > 0) {
     exon.info = exon.info[-nonCanonical, ]
}
chromosomes = as.character(levels(as.factor(as.character(exon.info$chr))))

aligned.reads = rep(0, times=length(sampleIDs))
exonic.reads = rep(0, times = length(sampleIDs))

bam.files = list.files(alignment.folder, pattern=".bam$")
if(length(grep("sort.bam",bam.files)) > 0){
	bam.files = bam.files[-grep("sort.bam",bam.files)]
}
sampleIDs = sub(".bam$","",bam.files)

for(i in 1:length(bam.files)){
	inputfile = paste(alignment.folder, bam.files[i], sep="/")
	print(inputfile)
	exon_reads <- list()
    total_reads <- list()
    
	for(chr in chromosomes){
    	print(chr)
        data <- readGAlignments(file = inputfile, use.names = TRUE, 
            					param = ScanBamParam(which = GRanges(chr, IRanges(1, chr_length[chr]))))
		total_reads[[chr]] <- names(data)
	}#end  for(chr in chromosomes)
	
	aligned.reads[i] = length(unique(unlist(total_reads))) 
	print(aligned.reads)
}#end for(i in 1:length(bam.files))

stat.table = data.frame(Sample = sampleIDs, aligned.reads = aligned.reads)
write.table(stat.table, aligned.stats.file, row.names=F, sep="\t", quote=F)
print(warnings())
