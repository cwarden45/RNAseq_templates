#if running TopHat with --fusion-search

sampleIDs = c("","","","")

alignmentFolder = "../hg19_TopHat_Alignment"

#use IGV coordinates

TMPRSS2.chr = "chr21"
TMPRSS2.start = 42834478
TMPRSS2.stop = 42880085

ERG.chr = "chr21"
ERG.start = 39737183
ERG.stop = 40035618

#try switching coordinates
ERG.start = 42834478
ERG.stop = 42880085

TMPRSS2.start = 39737183
TMPRSS2.stop = 40035618

fusion.strand = "rr"

min.span = 1

for (sampleID in sampleIDs){
	print(sampleID)
	fusion.junctions = paste(alignmentFolder,"/",sampleID,"/fusions.out",sep="")
	
	junction.table = read.table(fusion.junctions, head=F, sep="\t")
	print(dim(junction.table))
	junction.table = junction.table[junction.table$V5 >= min.span,]
	print(dim(junction.table))
	
	#look for chr with 5` to 3` strand
	junction.chr = paste(TMPRSS2.chr,ERG.chr,sep="-")
	#junction.table = junction.table[(junction.table$V1 == junction.chr)&(junction.table$V4 == fusion.strand), ]
	junction.table = junction.table[(junction.table$V1 == junction.chr), ]
	print(dim(junction.table))
	
	upstream.pos = junction.table$V2
	downstream.pos = junction.table$V3
	
	junction.table = junction.table[(TMPRSS2.start <= upstream.pos)&(upstream.pos <= TMPRSS2.stop), ]
	print(dim(junction.table))	

	upstream.pos = junction.table$V2
	downstream.pos = junction.table$V3
	
	junction.table = junction.table[(ERG.start <= downstream.pos)&(downstream.pos <= ERG.stop), ]
	print(dim(junction.table))	
	
	colnames(junction.table)[1]="fusion.chr"
	colnames(junction.table)[2]="upstream.pos"
	colnames(junction.table)[3]="downstream.pos"
	colnames(junction.table)[4]="fusion.orientation"
	colnames(junction.table)[5]="span.reads"
	colnames(junction.table)[6]="mate.pair"
	colnames(junction.table)[7]="mixed.support"
	colnames(junction.table)[8]="contradictory.reads"
	colnames(junction.table)[9]="cov.span.left"
	colnames(junction.table)[10]="cov.span.right"
	
	candidate.file = paste("TopHat2_TMPRSS2-ERG_",sampleID,".txt",sep="")
	write.table(junction.table, candidate.file, quote=F, sep="\t", row.names=F)
}#end for (sampleID in sampleIDs)

