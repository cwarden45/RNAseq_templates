#if running wtih --chimSegmentMin 12 --chimJunctionOverhangMin 12 

sampleIDs = c("","","","")

alignmentFolder = "../hg19_STAR_Alignment"

#use IGV coordinates

TMPRSS2.chr = "chr21"
TMPRSS2.start = 42834478
TMPRSS2.stop = 42880085
TMPRSS2.strand = "-"

ERG.chr = "chr21"
ERG.start = 39737183
ERG.stop = 40035618
ERG.strand = "-"

#ignore strand (although it would have to be switched for stranded library)

min.span = 3

for (sampleID in sampleIDs){
	print(sampleID)
	fusion.junctions = paste(alignmentFolder,"/",sampleID,"/",sampleID,"_Chimeric.out.junction",sep="")
	
	junction.table = read.table(fusion.junctions, head=F, sep="\t")
	print(dim(junction.table))
	
	upstream.pos = junction.table$V5
	downstream.pos = junction.table$V2
	
	expected.logical=((TMPRSS2.start <= upstream.pos)&(upstream.pos <= TMPRSS2.stop))&((ERG.start <= downstream.pos)&(downstream.pos <= ERG.stop))
	opposite.logical=((TMPRSS2.start <= downstream.pos)&(downstream.pos <= TMPRSS2.stop))&((ERG.start <= upstream.pos)&(upstream.pos <= ERG.stop))
	
	junction.table = junction.table[(expected.logical)|(opposite.logical), ]
	print(dim(junction.table))	
	
	junctionID = apply(junction.table[1:9],1,paste,collapse="\t")
	junction.counts = table(junctionID)
	
	print(length(junction.counts))
	junction.counts=junction.counts[junction.counts >= min.span]
	print(length(junction.counts))
	
	if(length(junction.counts)>0){
		junction.table=data.frame(junction.counts)
		if(length(junction.counts)==1){
			junction.table=data.frame(ID=names(junction.counts),junction.counts)
			
			colnames(junction.table)[1]="donor.chr\tdonor.pos\tdonor.strand\tacceptor.chr\tacceptor.pos\tacceptor.strand\tjunction.type\trepeat.left\trepeat.right"
			colnames(junction.table)[2]="junction.reads"
			
			candidate.file = paste("STAR_TMPRSS2-ERG_",sampleID,".txt",sep="")
			write.table(junction.table, candidate.file, quote=F, sep="\t", row.names=F, col.names=T)			
		}else{
			colnames(junction.table)[1]="donor.chr\tdonor.pos\tdonor.strand\tacceptor.chr\tacceptor.pos\tacceptor.strand\tjunction.type\trepeat.left\trepeat.right"
			colnames(junction.table)[2]="junction.reads"
			
			candidate.file = paste("STAR_TMPRSS2-ERG_",sampleID,".txt",sep="")
			write.table(junction.table, candidate.file, quote=F, sep="\t", row.names=F, col.names=T)
		}
	}#endif(length(junction.counts)>0)
}#end for (sampleID in sampleIDs)

