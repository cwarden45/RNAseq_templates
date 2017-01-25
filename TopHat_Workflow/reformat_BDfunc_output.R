#Useful if you want scores for individual samples (which means you didn't use group labels in header)

#You can use this if your enrichment file has one signature
#otherwise, you'll need to modify the code to extract a specific row

score.file = "[signal file]_[enrichment file]_sig.bdfunc"
histogram.file = "Score_Hist.pdf"
reformat.score.file = "Score_samples_in_rows.txt"

bdfunc.out = read.table(score.file, head=T, sep="\t")
scores = bdfunc.out[1,seq(2,ncol(bdfunc.out),4)]
samples = names(scores)
samples = gsub(".test.statistic","",samples)
samples = gsub("\\.","-",samples)
scores = as.numeric(scores)

pdf(histogram.file)
hist(scores, col="gray", main=paste("median= ",round(median(scores),digits=2),sep=""))
dev.off()

reformat.table = data.frame(sample=samples, score=scores)
write.table(reformat.table, reformat.score.file, quote=F, sep="\t", row.names=F)