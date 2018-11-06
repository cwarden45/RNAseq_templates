#new differentially expressed gene list
upID = " Up"
downID = " Down"

deg.file = ".txt"
deg.table = read.table(deg.file, head=T, sep="\t")
background.genes=as.character(deg.table$symbol)
deg.up = as.character(deg.table$symbol[deg.table$status == upID])
deg.down = as.character(deg.table$symbol[deg.table$status == downID])

#########################################################################################
#############					Earlier Gene Set					#####################
#########################################################################################

output.prefix = ""
set.name = ""

gs.upID = " Up"
gs.downID = " Down"

gs.file = "../../.txt"
gs.table = read.table(gs.file, head=T, sep="\t")

genes.up = as.character(gs.table$symbol[gs.table$status == gs.upID])
genes.down = as.character(gs.table$symbol[gs.table$status == gs.downID])

#### Ideally, modify above this point ####

##### Up-Regulated Genes ####
signature.name = c()
signature.pvalue = c()
overlap.genes = c()
num.overlap = c()
num.deg = c()
num.gene.set = c()


#up DEG, up GeneSet
matched.up = deg.up[match(genes.up, deg.up, nomatch=0)]
matched.background = intersect(genes.up, background.genes)
mat = matrix(c(length(matched.up), length(deg.up)-length(matched.up),
			length(matched.background),length(background.genes)-length(matched.background)),ncol=2)
result = fisher.test(mat, alternative="greater")
print(paste("Up DEG, ",set.name," Up --> ",result$p.value,sep=""))
print(length(matched.up))

signature.name = c(signature.name, paste(set.name,"_UP",sep=""))
signature.pvalue = c(signature.pvalue, result$p.value)
overlap.genes = c(overlap.genes, paste(matched.up,collapse=","))
num.overlap = c(num.overlap, length(matched.up))
num.deg = c(num.deg, length(deg.up))
num.gene.set = c(num.gene.set, length(genes.up))

#up DEG, down GeneSet
matched.up = deg.up[match(genes.down, deg.up, nomatch=0)]
matched.background = intersect(genes.down, background.genes)
mat = matrix(c(length(matched.up), length(deg.up)-length(matched.up),
			length(matched.background),length(background.genes)-length(matched.background)),ncol=2)
result = fisher.test(mat, alternative="greater")
print(paste("Up DEG, ",set.name," Down --> ",result$p.value,sep=""))
print(length(matched.up))

signature.name = c(signature.name, paste(set.name,"_DOWN",sep=""))
signature.pvalue = c(signature.pvalue, result$p.value)
overlap.genes = c(overlap.genes, paste(matched.up,collapse=","))
num.overlap = c(num.overlap, length(matched.up))
num.deg = c(num.deg, length(deg.up))
num.gene.set = c(num.gene.set, length(genes.up))

output.file = paste(output.prefix,"_UP.txt",sep="")
output.table = data.frame(signature.name, signature.pvalue,
							overlap.genes, num.overlap, num.deg, num.gene.set)
write.table(output.table, output.file, quote=F, sep="\t", row.names=F)

##### Down-Regulated Genes ####
signature.name = c()
signature.pvalue = c()
overlap.genes = c()
num.overlap = c()
num.deg = c()
num.gene.set = c()


#down DEG, up GeneSet
matched.down = deg.down[match(genes.up, deg.down, nomatch=0)]
matched.background = intersect(genes.up, background.genes)
mat = matrix(c(length(matched.down), length(deg.down)-length(matched.down),
			length(matched.background),length(background.genes)-length(matched.background)),ncol=2)
result = fisher.test(mat, alternative="greater")
print(paste("Down DEG, ",set.name," Up --> ",result$p.value,sep=""))
print(length(matched.down))

signature.name = c(signature.name, paste(set.name,"_UP",sep=""))
signature.pvalue = c(signature.pvalue, result$p.value)
overlap.genes = c(overlap.genes, paste(matched.down,collapse=","))
num.overlap = c(num.overlap, length(matched.down))
num.deg = c(num.deg, length(deg.down))
num.gene.set = c(num.gene.set, length(genes.down))

#down DEG, down GeneSet
matched.down = deg.down[match(genes.down, deg.down, nomatch=0)]
matched.background = intersect(genes.down, background.genes)
mat = matrix(c(length(matched.down), length(deg.down)-length(matched.down),
			length(matched.background),length(background.genes)-length(matched.background)),ncol=2)
result = fisher.test(mat, alternative="greater")
print(paste("Down DEG, ",set.name," Up --> ",result$p.value,sep=""))
print(length(matched.down))

signature.name = c(signature.name, paste(set.name,"_DOWN",sep=""))
signature.pvalue = c(signature.pvalue, result$p.value)
overlap.genes = c(overlap.genes, paste(matched.down,collapse=","))
num.overlap = c(num.overlap, length(matched.down))
num.deg = c(num.deg, length(deg.down))
num.gene.set = c(num.gene.set, length(genes.down))

output.file = paste(output.prefix,"_DOWN.txt",sep="")
output.table = data.frame(signature.name, signature.pvalue,
							overlap.genes, num.overlap, num.deg, num.gene.set)
write.table(output.table, output.file, quote=F, sep="\t", row.names=F)