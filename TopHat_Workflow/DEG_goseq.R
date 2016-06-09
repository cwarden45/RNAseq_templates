normalizeTotalExpression <- function (geneExpr, totalReads) {
	return(geneExpr / totalReads)
}#end def normalizeTotalExpression

avgGroupExpression <- function (geneExpr, groups) {
	avg.expr = tapply(geneExpr, groups, mean)
	return(avg.expr)
}#end def avgGroupExpression

standardize.arr <- function(arr)
{
	center.arr = as.numeric(arr) - mean(as.numeric(arr), na.rm=T)
	norm.arr = center.arr / sd(center.arr, na.rm=T)
	return(norm.arr)
}#end def standardize.arr

count.defined.values <- function(arr, expr.cutoff)
{
	return(length(arr[arr > expr.cutoff]))
}#end def count.values

count.na.values <- function(arr)
{
	return(length(arr[is.na(arr)]))
}#end def count.values

ratio2fc <- function(value)
{
	if(value >= 0){
		return(2^value)
	} else {
		return (-2^(-value))
	}
}#end def ratio2fc

gene.lm <- function(arr, var1, var2=c(), var3=c())
{	
	if (length(var2) == 0){
		fit = lm(as.numeric(arr) ~ var1)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else if (length(var3) == 0){
		fit = lm(as.numeric(arr) ~ var1 + var2)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else {
		fit = lm(as.numeric(arr) ~ var2*var3 + var2 + var3)
		result = summary(fit)
		pvalue = result$coefficients[4,4]
	}
	return(pvalue)
}#end def gene.lm

gene.aov <- function(arr, var1, var2=c(), var3=c())
{	
	if (length(var2) == 0){
		fit = aov(as.numeric(arr) ~ var1)
		result = summary(fit)
		aov.pvalue = result[[1]][['Pr(>F)']][1]
	} else if (length(var3) == 0){
		fit = aov(as.numeric(arr) ~ var1 + var2)
		result = summary(fit)
		aov.pvalue = result[[1]][['Pr(>F)']][1]
	} else {
		fit = aov(as.numeric(arr) ~ var2*var3 + var2 + var3)
		result = summary(fit)
		aov.pvalue = result[[1]][['Pr(>F)']][3]
	}
	return(pvalue)
}#end def gene.aov

library(gplots)
fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black",colors())

param.table = read.table("parameters.txt", header=T, sep="\t")
comp.name=as.character(param.table$Value[param.table$Parameter == "comp_name"])
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
min.expression = as.numeric(as.character(param.table$Value[param.table$Parameter == "rpkm_expression_cutoff"]))
min.fraction.expressed = as.numeric(as.character(param.table$Value[param.table$Parameter == "minimum_fraction_expressed"]))
fc.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "fold_change_cutoff"]))
pvalue.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "pvalue_cutoff"]))
fdr.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "fdr_cutoff"]))
fc.cutoff2 = as.numeric(as.character(param.table$Value[param.table$Parameter == "sec_fold_change_cutoff"]))
pvalue.cutoff2 = as.numeric(as.character(param.table$Value[param.table$Parameter == "sec_pvalue_cutoff"]))
fdr.cutoff2 = as.numeric(as.character(param.table$Value[param.table$Parameter == "sec_fdr_cutoff"]))
gene.length.cutoff.kb = as.numeric(as.character(param.table$Value[param.table$Parameter == "min_length_kb"]))
deg.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "deg_groups"]), split=","))
plot.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "plot_groups"]), split=","))
trt.group = as.character(param.table$Value[param.table$Parameter == "treatment_group"])
trt.group2 = as.character(param.table$Value[param.table$Parameter == "secondary_trt"])
interaction.flag = as.character(param.table$Value[param.table$Parameter == "interaction"])
pvalue.method = as.character(param.table$Value[param.table$Parameter == "pvalue_method"])
goseq.flag = as.character(param.table$Value[param.table$Parameter == "run_goseq"])
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
counts.file = as.character(param.table$Value[param.table$Parameter == "counts_file"])
aligned.stats.file = as.character(param.table$Value[param.table$Parameter == "aligned_stats_file"])

setwd(output.folder)

sample.description.table = read.table(sample.description.file, sep="\t", header=T)
longID = sample.description.table$sampleID
sample.label = sample.description.table$userID

deg.group.table = sample.description.table[,deg.groups]
if (length(deg.groups) == 1){
	deg.meta = sample.description.table[!is.na(deg.group.table),]
} else {
	deg.grp.na.counts = apply(deg.group.table, 1, count.na.values)
	deg.meta = sample.description.table[deg.grp.na.counts == 0,]
}

counts.table = read.table(counts.file, head=T, sep="\t")
counts= counts.table[,match(sample.label,names(counts.table))]
counts = matrix(as.numeric(unlist(counts)), ncol=length(sample.label))

counts =counts[rowSums(!is.na(counts)) == length(sample.label),]
counts.table = counts.table[rowSums(!is.na(counts)) == length(sample.label),]
print(dim(counts.table))
print(dim(counts))

counts =counts[as.numeric(counts.table$length.kb) > gene.length.cutoff.kb,]
counts.table = counts.table[as.numeric(counts.table$length.kb) > gene.length.cutoff.kb,]
print(dim(counts.table))
print(dim(counts))

genes = counts.table$symbol
gene.length.kb = as.numeric(counts.table$length.kb)

exonic.stat.table = read.table(aligned.stats.file, header=T, sep="\t")
aligned.reads = as.numeric(exonic.stat.table$aligned.reads[match(longID, exonic.stat.table$Sample)])

total.million.aligned.reads = aligned.reads / 1000000
print(total.million.aligned.reads)

rpk = matrix(ncol=ncol(counts), nrow=nrow(counts))
for (i in 1:ncol(counts)){
	temp.counts = as.numeric(counts[,i])
	temp.rpk = temp.counts / gene.length.kb
	rpk[,i] = temp.rpk 
}
RPKM = log2(round(t(apply(rpk, 1, normalizeTotalExpression, totalReads = total.million.aligned.reads)), digits=2) + min.expression)
colnames(RPKM) = as.character(sample.label)

RPKM = RPKM[,match(sample.label, colnames(RPKM))]

print(dim(counts))
expressed.sample.count = apply(RPKM, 1, count.defined.values, expr.cutoff = log2(min.expression))
counts.table = counts.table[expressed.sample.count > round(min.fraction.expressed * ncol(counts)),]
counts = counts[expressed.sample.count > round(min.fraction.expressed * ncol(counts)),]
RPKM = RPKM[expressed.sample.count > round(min.fraction.expressed * ncol(counts)),]
print(dim(counts))

genes = counts.table$symbol
gene.length.kb = as.numeric(counts.table$length.kb)

if(length(plot.groups) == 1){
	print("Averaging Expression for One Variable (for plot.groups)")
	grp = sample.description.table[,plot.groups]
} else if ((length(plot.groups) == 2)&(interaction.flag == "no")){
	print("Averaging Expression for First Variable (for plot.groups)")
	grp = sample.description.table[,plot.groups[1]]
} else if (length(plot.groups) == 2){
	print("Averaging Expression for Interaction Variable (for plot.groups)")
	grp = paste(sample.description.table[,plot.groups[1]],sample.description.table[,plot.groups[2]],sep=":")
} else {
	stop("Code only compatible with 2 variables (with or without a 3rd interaction variable")
}

groupIDs = as.character(levels(as.factor(grp)))
average.rpkm = data.frame(t(apply(RPKM, 1, avgGroupExpression, groups = grp)))
average.rpkm = average.rpkm
colnames(average.rpkm) = paste("avg.log2.rpkm", sub("-",".",groupIDs), sep=".")

if(length(deg.groups) == 1){
	print("Averaging Expression for One Variable (for deg.groups)")
	contrast.grp = sample.description.table[,deg.groups]
} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
	print("Averaging Expression for First Variable (for deg.groups)")
	contrast.grp = sample.description.table[,deg.groups[1]]
} else if (length(deg.groups) == 2){
	print("Averaging Expression for Interaction Variable (for deg.groups)")
	contrast.grp = paste(sample.description.table[,deg.groups[1]],sample.description.table[,deg.groups[2]],sep=":")
} else {
	stop("Code only compatible with 2 variables (with or without a 3rd interaction variable")
}

groupIDs = as.character(levels(as.factor(contrast.grp)))
contrast.rpkm = data.frame(t(apply(RPKM, 1, avgGroupExpression, groups = contrast.grp)))
colnames(contrast.rpkm) = paste("avg.log2.rpkm", sub("-",".",groupIDs), sep=".")

if(interaction.flag == "no"){
trt.expr = contrast.rpkm[,paste("avg.log2.rpkm", sub("-",".",trt.group), sep=".")]
cntl.expr = contrast.rpkm[,paste("avg.log2.rpkm", sub("-",".",groupIDs[groupIDs != trt.group]), sep=".")]

log2ratio = round(trt.expr - cntl.expr, digits = 2)
fc = round(sapply(log2ratio, ratio2fc), digits = 2)
fc.table = data.frame(log2ratio=log2ratio, fold.change=fc)
} else if ((interaction.flag == "model")|(interaction.flag == "filter-overlap")){
prim.groups = as.character(levels(as.factor(sample.description.table[,deg.groups[1]])))
prim.trt = paste(trt.group,trt.group2,sep=":")
prim.cntl = paste(prim.groups[prim.groups != trt.group],trt.group2,sep=":")
prim.trt.expr = contrast.rpkm[,paste("avg.log2.rpkm", sub("-",".",prim.trt), sep=".")]
prim.cntl.expr = contrast.rpkm[,paste("avg.log2.rpkm", sub("-",".",prim.cntl), sep=".")]

prim.log2ratio = round(prim.trt.expr - prim.cntl.expr, digits = 2)
prim.fc = round(sapply(prim.log2ratio, ratio2fc), digits = 2)

sec.groups = as.character(levels(as.factor(sample.description.table[,deg.groups[2]])))
sec.trt = paste(trt.group, sec.groups[sec.groups != trt.group2], sep=":")
sec.cntl = paste(prim.groups[prim.groups != trt.group], sec.groups[sec.groups != trt.group2], sep=":")
sec.trt.expr = contrast.rpkm[,paste("avg.log2.rpkm", sub("-",".",sec.trt), sep=".")]
sec.cntl.expr = contrast.rpkm[,paste("avg.log2.rpkm", sub("-",".",sec.cntl), sep=".")]

sec.log2ratio = round(sec.trt.expr - sec.cntl.expr, digits = 2)
sec.fc = round(sapply(sec.log2ratio, ratio2fc), digits = 2)

overall.log2ratio = prim.log2ratio - sec.log2ratio
overall.fc = round(sapply(overall.log2ratio, ratio2fc), digits = 2)

fc.table = data.frame(fc1 = prim.fc, fc2=sec.fc, fc3=overall.fc)
colnames(fc.table) = c(paste("fold.change",trt.group,":",trt.group2,sep="."),
						paste("fold.change",trt.group,":",sec.groups[sec.groups != trt.group2], sep="."),
						"overall.fold.change")
}else{
		stop("interaction must be \"no\", \"model\", or \"filter-overlap\"")
}#end else

rep.check = 1
for (deg.group in deg.groups){
	deg.group.values = deg.meta[,deg.group]
	min.reps = min(table(deg.group.values))
	if (min.reps < 2){
		rep.check=0
		print("There are not at least 2 samples per-group.  P-value will not be calculated.")
		print("In the future, please make sure you at least have duplicate samples.")
	}#end if (min.reps < 2)
}#end for (deg.group in deg.groups)

if(rep.check == 1){
	if (pvalue.method == "edgeR"){
		library(edgeR)
		
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("edgeR, Two-Step Analysis")
			prim.counts = counts[,sample.description.table[,deg.groups[2]]==sec.groups[sec.groups==trt.group2]]
			prim.edgeR.grp = sample.description.table[,deg.groups[1]][sample.description.table[,deg.groups[2]]==sec.groups[sec.groups==trt.group2]]

			y <- DGEList(counts=prim.counts, genes=genes)
			y = estimateCommonDisp(y)
			design <- model.matrix(~prim.edgeR.grp)
			fit <- glmFit(y, design)
			lrt <- glmLRT(fit, coef=2)

			prim.pvalue = lrt$table$PValue
			
			sec.counts = counts[,sample.description.table[,deg.groups[2]]==sec.groups[sec.groups!=trt.group2]]
			sec.edgeR.grp = sample.description.table[,deg.groups[1]][sample.description.table[,deg.groups[2]]!=sec.groups[sec.groups==trt.group2]]

			y <- DGEList(counts=sec.counts, genes=genes)
			y = estimateCommonDisp(y)
			design <- model.matrix(~sec.edgeR.grp)
			fit <- glmFit(y, design)
			lrt <- glmLRT(fit, coef=2)

			sec.pvalue = lrt$table$PValue
			
		} else {
			y <- DGEList(counts=counts, genes=genes)
			y = estimateCommonDisp(y)
			if (length(deg.groups) == 1){
				print("edgeR with 1 variable")
				var1 = sample.description.table[,deg.groups]
				design <- model.matrix(~var1)
				fit <- glmFit(y, design)
				lrt <- glmLRT(fit, coef=2)
				test.pvalue = lrt$table$PValue
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("edgeR with 2 variables")
				var1 = sample.description.table[,deg.groups[1]]
				var2 = sample.description.table[,deg.groups[2]]
				design <- model.matrix(~var1 + var2)
				fit <- glmFit(y, design)
				lrt <- glmLRT(fit, coef=2)
				test.pvalue = lrt$table$PValue
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("edgeR with 2 variables plus interaction")
				var1 = sample.description.table[,deg.groups[1]]
				var2 = sample.description.table[,deg.groups[2]]
				design <- model.matrix(~var1*var2 + var1 + var2)
				fit <- glmFit(y, design)
				lrt <- glmLRT(fit, coef=4)
				test.pvalue = lrt$table$PValue
			}
		}#end else
	} else if (pvalue.method == "lm"){
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("RPKM linear regression, Two-Step Analysis")
			prim.RPKM = RPKM[,sample.description.table[,deg.groups[2]]==sec.groups[sec.groups==trt.group2]]
			prim.lm.grp = sample.description.table[,deg.groups[1]][sample.description.table[,deg.groups[2]]==sec.groups[sec.groups==trt.group2]]

			prim.pvalue = apply(prim.RPKM, 1, gene.lm, var1=prim.lm.grp)
			
			sec.RPKM = RPKM[,sample.description.table[,deg.groups[2]]==sec.groups[sec.groups!=trt.group2]]
			sec.lm.grp = sample.description.table[,deg.groups[1]][sample.description.table[,deg.groups[2]]!=sec.groups[sec.groups==trt.group2]]
			
			sec.pvalue = apply(sec.RPKM, 1, gene.lm, var1=sec.lm.grp)
		} else {
			if (length(deg.groups) == 1){
				print("RPKM linear regression with 1 variable")
				var1 = sample.description.table[,deg.groups]
				test.pvalue = apply(RPKM, 1, gene.lm, var1=var1)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("RPKM linear regression with 2 variables")
				var1 = sample.description.table[,deg.groups[1]]
				var2 = sample.description.table[,deg.groups[2]]
				test.pvalue = apply(RPKM, 1, gene.lm, var1=var1, var2=var2)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("RPKM linear regression with 2 variables plus interaction")
				var1 = as.factor(paste(sample.description.table[,deg.groups[1]],sample.description.table[,deg.groups[2]],sep=":"))
				var2 = sample.description.table[,deg.groups[1]]
				var3 = sample.description.table[,deg.groups[2]]
				test.pvalue = apply(RPKM, 1, gene.lm, var1=var1, var2=var2, var3=var3)
			}
		}#end else
	} else if (pvalue.method == "ANOVA"){
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("RPKM ANOVA, Two-Step Analysis")
			prim.RPKM = RPKM[,sample.description.table[,deg.groups[2]]==sec.groups[sec.groups==trt.group2]]
			prim.aov.grp = sample.description.table[,deg.groups[1]][sample.description.table[,deg.groups[2]]==sec.groups[sec.groups==trt.group2]]

			prim.pvalue = apply(prim.RPKM, 1, gene.aov, var1=prim.aov.grp)
			
			sec.RPKM = RPKM[,sample.description.table[,deg.groups[2]]==sec.groups[sec.groups!=trt.group2]]
			sec.aov.grp = sample.description.table[,deg.groups[1]][sample.description.table[,deg.groups[2]]!=sec.groups[sec.groups==trt.group2]]
			
			sec.pvalue = apply(sec.RPKM, 1, gene.aov, var1=sec.aov.grp)
		} else {
			if (length(deg.groups) == 1){
				print("RPKM ANOVA with 1 variable")
				var1 = sample.description.table[,deg.groups]
				test.pvalue = apply(RPKM, 1, gene.aov, var1=var1)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("RPKM ANOVA with 2 variables")
				var1 = sample.description.table[,deg.groups[1]]
				var2 = sample.description.table[,deg.groups[2]]
				test.pvalue = apply(RPKM, 1, gene.aov, var1=var1, var2=var2)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("RPKM ANOVA with 2 variables plus interaction")
				var1 = as.factor(paste(sample.description.table[,deg.groups[1]],sample.description.table[,deg.groups[2]],sep=":"))
				var2 = sample.description.table[,deg.groups[1]]
				var3 = sample.description.table[,deg.groups[2]]
				test.pvalue = apply(RPKM, 1, gene.aov, var1=var1, var2=var2, var3=var3)
			}
		}#end else
	} else{
		stop("pvalue_method must be \"edgeR\", \"lm\", or \"ANOVA\"")
	}
} else{
	test.pvalue = rep(1,times=length(genes))
	prim.pvalue = rep(1,times=length(genes))
	sec.pvalue = rep(1,times=length(genes))
}#end else

if (interaction.flag == "no"){
	fdr = p.adjust(test.pvalue, "fdr")
	status = rep("No Change", times=length(fdr))
	status[(fc >= fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Up",sep="")
	status[(fc <= -fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Down",sep="")
} else{
	trt.group = prim.trt
	if(interaction.flag == "model"){
		fdr = p.adjust(test.pvalue, "fdr")
		status = rep("No Change", times=length(fdr))
		status[(overall.fc >= fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Up",sep="")
		status[(overall.fc <= -fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Down",sep="")
		pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
	} else if (interaction.flag == "filter-overlap"){
		fdr = p.adjust(prim.pvalue, "fdr")
		pass1.status = rep("No Change", times=length(fdr))
		pass1.status[(prim.fc >= fc.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Up",sep="")
		pass1.status[(prim.fc <= -fc.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Down",sep="")

		print(paste("Primary Up-Regulated: ",length(pass1.status[pass1.status == paste(trt.group," Up",sep="")]),sep=""))
		print(paste("Primary Down-Regulated: ",length(pass1.status[pass1.status == paste(trt.group," Down",sep="")]),sep=""))
			
		sec.fdr = p.adjust(sec.pvalue, "fdr")
		pass2.status = rep("No Change", times=length(fdr))
		pass2.status[(sec.fc >= fc.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Up",sep="")
		pass2.status[(sec.fc <= -fc.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Down",sep="")

		print(paste("Secondary Up-Regulated: ",length(pass2.status[pass2.status == paste(trt.group," Up",sep="")]),sep=""))
		print(paste("Secondary Down-Regulated: ",length(pass2.status[pass2.status == paste(trt.group," Down",sep="")]),sep=""))
			
		pvalue.table = data.frame(prim.pvalue = prim.pvalue, prim.FDR = fdr,
										sec.pvalue=sec.pvalue, sec.fdr=sec.fdr)
			
		status = rep("No Change", times=length(fdr))
		status[(pass1.status == paste(trt.group," Up",sep="")) & (pass2.status == "No Change")] = paste(trt.group," Up",sep="")
		status[(pass1.status == paste(trt.group," Down",sep="")) & (pass2.status == "No Change")] = paste(trt.group," Down",sep="")
	} else{
		stop("interaction must be \"no\", \"model\", or \"filter-overlap\"")
	}#end else
}#end else

print(paste("Up-Regulated: ",length(status[status == paste(trt.group," Up",sep="")]),sep=""))
print(paste("Down-Regulated: ",length(status[status == paste(trt.group," Down",sep="")]),sep=""))

if (interaction.flag == "filter-overlap"){
	pvalue.method = paste(pvalue.method,"two-step_filtered",sep="_")
}

if(rep.check == 1){
	deg.table = data.frame(symbol = genes, gene.length.kb=gene.length.kb,
							average.rpkm, fc.table,
							pvalue.table, status = status)
} else {
	deg.table = data.frame(symbol = genes, gene.length.kb=gene.length.kb,
							average.rpkm, fc.table, status = status)	
}#end else
deg.file = paste(comp.name,"_",pvalue.method,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".txt",sep="")
deg.file = gsub(":",".",deg.file)
write.table(deg.table, file=deg.file, row.names=F, quote=F, sep="\t")

final.deg.file = paste(user.folder,"/DEG/",comp.name,"_DEG_stats.txt",sep="")
write.table(deg.table, file=final.deg.file, row.names=F, quote=F, sep="\t")

temp.rpkm = RPKM
temp.rpkm = temp.rpkm[status != "No Change", ]
deg.genes = genes[status != "No Change"]

if(length(plot.groups) > 1){
	source("heatmap.3.R")
	grp1 = as.character(sample.description.table[,plot.groups[1]])
	grp2 = as.character(sample.description.table[,plot.groups[2]])
	group.levels = c(levels(as.factor(grp1)),levels(as.factor(grp2)))

	color.palette <- fixed.color.palatte[1:length(group.levels)]
	labelColors1 = rep("black",times=length(sample.label))
	for (i in 1:length(group.levels)){
		labelColors1[grp1 == as.character(group.levels[i])] = color.palette[i]
	}#end for (i in 1:length(group.levels))
	labelColors2 = rep("black",times=length(sample.label))
	for (i in 1:length(group.levels)){
		labelColors2[grp2 == as.character(group.levels[i])] = color.palette[i]
	}#end for (i in 1:length(group.levels))
	
	std.expr = apply(temp.rpkm, 1, standardize.arr)
	if(length(deg.genes) < 25){
		colnames(std.expr) = deg.genes
	} else {
		colnames(std.expr) = rep("", length(deg.genes))
	}
	rownames(std.expr) = sample.label

	column_annotation <- as.matrix(deg.genes)
	colnames(column_annotation) <- c("")

	row_annotation <- data.frame(label1 = labelColors1, label2 = labelColors2)
	row_annotation = as.matrix(t(row_annotation))
	rownames(row_annotation) <- c(plot.groups)

	heatmap.file <- paste(comp.name,"_",pvalue.method,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
	heatmap.file = gsub(":",".",heatmap.file)
	png(file = heatmap.file)
	heatmap.3(std.expr, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
				RowSideColors=row_annotation, trace="none", margins = c(8,13),RowSideColorsSize=4, dendrogram="both")
	legend("topright", legend=group.levels,
						col=color.palette,
						pch=15, cex=0.7)
	dev.off()
	
	if(interaction.flag != "none"){
		temp.fc.table = as.matrix(fc.table[,-ncol(fc.table)])
		temp.fc.table = temp.fc.table[status != "No Change", ]
		if(length(deg.genes) < 25){
			rownames(temp.fc.table) = deg.genes
		} else {
			rownames(temp.fc.table) = rep("",times=length(deg.genes))
		}
		colnames(temp.fc.table) = gsub(".:.",":",gsub("fold.change.","",colnames(temp.fc.table)))
	
		temp.fc.table[temp.fc.table < -10] = -10
		temp.fc.table[temp.fc.table > 10] = 10
	
		heatmap.file <- paste("fold_change_",comp.name,"_",pvalue.method,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
		heatmap.file = gsub(":",".",heatmap.file)
		png(file = heatmap.file)
		heatmap.2(temp.fc.table, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					trace="none", margins = c(20,5), cexCol=1.5)
		dev.off()
	}#end if(interaction.flag != "none")
	
} else {
	group.levels = levels(as.factor(sample.description.table[,plot.groups]))

	color.palette <- fixed.color.palatte[1:length(group.levels)]
	labelColors = rep("black",times=length(sample.label))
	for (i in 1:length(group.levels)){
		labelColors[grp == as.character(group.levels[i])] = color.palette[i]
	}#end for (i in 1:length(group.levels))

	std.expr = apply(temp.rpkm, 1, standardize.arr)
	if(length(deg.genes) < 25){
		colnames(std.expr) = deg.genes
	} else {
		colnames(std.expr) = rep("", length(deg.genes))
	}
	rownames(std.expr) = sample.label
	
	heatmap.file <- paste(comp.name,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
	heatmap.file = gsub(":",".",heatmap.file)
	png(file = heatmap.file)
	heatmap.2(std.expr, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
				 RowSideColors=labelColors, trace="none", margins = c(5,15))
	dev.off()
}#end else

#goseq
if (goseq.flag == "yes"){
	library(goseq)
	
	deg = as.integer(status == paste(trt.group," Up",sep=""))
	deg = tapply(deg, as.character(genes), sum)
	deg[deg >= 1] = 1
	gene.symbol = as.character(levels(as.factor(as.character(genes))))
	names(deg)=gene.symbol
	
	bias.file = paste(comp.name,"_goseq_up_bias.png",sep="")
	png(bias.file)
	pwf=nullp(deg,genome,"geneSymbol")
	GO.wall=goseq(pwf,genome,"geneSymbol")
	GO.wall = data.frame(GO.wall, overenrichment.fdr=p.adjust(GO.wall$over_represented_pvalue))

	GO.wall = GO.wall[GO.wall$over_represented_pvalue < 0.05,]
	genes2go=getgo(gene.symbol[deg == 1],genome,"geneSymbol")
	go2genes=goseq:::reversemapping(genes2go)
	go.genes = rep("", times=nrow(GO.wall))
	dev.off()

	for( i in 1:length(go.genes)){
		go.category = GO.wall$category[i]
		result = paste(go2genes[[go.category]], sep="", collapse=",")
		go.genes[i] = result
	}#end for( i in 1:length(go.genes))

	go.table = data.frame(GO.wall, genes=go.genes)
	go.file = paste(user.folder,"/GO/",comp.name,"_goseq_UP.txt",sep="")
	write.table(go.table, file=go.file, row.names=F, quote=F, sep="\t")

	deg = as.integer(status == paste(trt.group," Down",sep=""))
	deg = tapply(deg, as.character(genes), sum)
	deg[deg >= 1] = 1
	gene.symbol = as.character(levels(as.factor(as.character(genes))))
	names(deg)=gene.symbol

	bias.file = paste(comp.name,"_goseq_down_bias.png",sep="")
	png(bias.file)
	pwf=nullp(deg,genome,"geneSymbol")
	GO.wall=goseq(pwf,genome,"geneSymbol")
	GO.wall = data.frame(GO.wall, overenrichment.fdr=p.adjust(GO.wall$over_represented_pvalue))

	GO.wall = GO.wall[GO.wall$over_represented_pvalue < 0.05,]
	genes2go=getgo(gene.symbol[deg == 1],genome,"geneSymbol")
	go2genes=goseq:::reversemapping(genes2go)
	go.genes = rep("", times=nrow(GO.wall))
	dev.off()

	for( i in 1:length(go.genes)){
		go.category = GO.wall$category[i]
		result = paste(go2genes[[go.category]], sep="", collapse=",")
		go.genes[i] = result
	}#end for( i in 1:length(go.genes))

	go.table = data.frame(GO.wall, genes=go.genes)
	go.file = paste(user.folder,"/GO/",comp.name,"_goseq_DOWN.txt",sep="")
	write.table(go.table, file=go.file, row.names=F, quote=F, sep="\t")
}#end if (goseq.flag == "yes")
