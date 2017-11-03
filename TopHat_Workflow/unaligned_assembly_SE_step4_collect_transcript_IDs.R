folder = "Oases_Assemblies"
num.top.transcripts = 2
transcript.list = paste(folder,"_top_",num.top.transcripts,"_transcript_IDs_for_BLAST.txt",sep="")

#folder = "Oases_Assemblies_Trim"
#num.top.transcripts = 2
#transcript.list = paste(folder,"_top_",num.top.transcripts,"_transcript_IDs_for_BLAST_Trim.txt",sep="")


subfolders = list.dirs(folder)

samples = c()
transcripts = c()
TPM = c()

for (subfolder in subfolders){
	TPM.file = file.path(subfolder, "results.xprs")
	if(file.exists(TPM.file)){
		sampleID = gsub(paste(folder,"/",sep=""),"",subfolder)
		print(sampleID)
		TPM.table = read.table(TPM.file, head=T, sep="\t")
		TPM.table=TPM.table[order(TPM.table$tpm, decreasing=T),]
		samples = c(samples, rep(sampleID,num.top.transcripts))
		transcripts = c(transcripts, as.character(TPM.table$target_id[1:num.top.transcripts]))
		TPM = c(TPM, TPM.table$tpm[1:num.top.transcripts])
	}
}#end for (subfolder in subfolders)

output.table = data.frame(Sample=samples, TranscriptID=transcripts, Transcript.TPM=TPM)
write.table(output.table, transcript.list, quote=F, sep="\t", row.names=F)