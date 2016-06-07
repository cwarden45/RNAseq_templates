library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)

setwd("R:\\path\\to\folder")

alignment.folder = "R:\\path\\to\folder\\"
genome="hg19"
length.file = paste("R:\\path\\to\\",genome,"_chr_length.txt",sep="")
full.annotation.file = paste("R:\\path\\to\\TxDb_",genome,"_exon_annotations.txt",sep="")
ignore.strand = T

length.table = read.table(length.file, header=T, sep="\t")
chr_length = as.numeric(length.table$Length)
names(chr_length) = as.character(length.table$Chr)

total.reads.file = "total_read_counts.txt"
total.reads.table = read.table(total.reads.file, header=T, sep="\t")
sampleIDs = as.character(total.reads.table$Sample)

exon.info = read.table(full.annotation.file, header=T, sep="\t")
    
#remove non-canonical chromosomes
bad <- grep("_", exon.info$chr)
 if (length(bad) > 0) {
     exon.info = exon.info[-bad, ]
}
chromosomes = as.character(levels(as.factor(as.character(exon.info$chr))))

aligned.reads = rep(0, times=length(sampleIDs))
exonic.reads = rep(0, times = length(sampleIDs))

bam.files <- list.files(alignment.folder, pattern=".bam$")
sampleIDs = sub(".bam$","",bam.files)

for(i in 1:length(bam.files)){
	inputfile = paste(alignment.folder, bam.files[i], sep="")
	print(inputfile)
	exon_reads <- list()
    total_reads <- list()
    
	for(chr in chromosomes){
    	print(chr)
        data <- readGAlignments(file = inputfile, use.names = TRUE, 
            					param = ScanBamParam(which = GRanges(chr, IRanges(1, chr_length[chr]))))
		total_reads[[chr]] <- names(data)

        chr.gene.table <- exon.info[which(exon.info$chr == i), ]
       	chr.gene.table <- exon.info[order(exon.info$start, exon.info$end), ]
        exon_range <- GRanges(seqnames = exon.info$chr, ranges=IRanges(exon.info$start, exon.info$end), strand = Rle(strand(exon.info$strand)))
        #findOverlaps will not distinguish between SE and PE reads
        hits <- findOverlaps(as(data, "GRanges"), exon_range, ignore.strand=ignore.strand)
        exon_reads[[chr]] <- names(data[(unique(as.matrix(hits)[, 1]))])
	}#end  for(chr in chromosomes)
	
	exonic.reads[i] = length(unique(unlist(exon_reads)))
	aligned.reads[i] = length(unique(unlist(total_reads))) 
	print(exonic.reads)
	print(aligned.reads)
}#end for(i in 1:length(bam.files))

percent.exonic = round(100 * exonic.reads / aligned.reads, digits = 1)
stat.table = data.frame(Sample = sampleIDs, aligned.reads = aligned.reads, exonic.reads = exonic.reads, percent.exonic=percent.exonic)
write.table(stat.table, "findOverlaps_exonic_stats.txt", row.names=F, sep="\t", quote=F)
print(warnings())