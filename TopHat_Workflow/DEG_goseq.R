normalizeTotalExpression = function (geneExpr, totalReads) {
	return(geneExpr / totalReads)
}#end def normalizeTotalExpression

avgGroupExpression = function (geneExpr, groups) {
	avg.expr = tapply(geneExpr, groups, mean)
	return(avg.expr)
}#end def avgGroupExpression

standardize.arr = function(arr)
{
	center.arr = as.numeric(arr) - mean(as.numeric(arr), na.rm=T)
	norm.arr = center.arr / sd(center.arr, na.rm=T)
	return(norm.arr)
}#end def standardize.arr

center.arr = function(arr)
{
	center.arr = as.numeric(arr) - mean(as.numeric(arr), na.rm=T)
	return(center.arr)
}#end def center.arr

count.defined.values = function(arr, expr.cutoff)
{
	sig.figures = 1
	if (expr.cutoff >= 0)
		sig.figures = 0
	expr.cutoff = round(expr.cutoff, digits=sig.figures)
	arr = round(arr, digits=sig.figures)
	return(length(arr[arr > expr.cutoff]))
}#end def count.values

count.na.values = function(arr)
{
	return(length(arr[is.na(arr)]))
}#end def count.values

ratio2fc = function(value)
{
	if(value >= 0){
		return(2^value)
	} else {
		return (-2^(-value))
	}
}#end def ratio2fc

gene.lm = function(arr, var1, var2=c(), var3=c())
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

gene.aov = function(arr, var1, var2=c(), var3=c())
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
	return(aov.pvalue)
}#end def gene.aov

calc.gene.cor = function(arr, indep.var)
{	
	na.count = length(arr[!is.na(arr)])
	if((na.count >= 3) & (sd(arr) != 0)){
		gene.cor.coef = cor(arr,indep.var)
	} else {
		gene.cor.coef = NA
	}
	return(gene.cor.coef)
}#end def calc.gene.cor

cor.dist = function(mat){
	cor.mat = cor(as.matrix(t(mat)))
	dis.mat = 1 - cor.mat
	return(as.dist(dis.mat))	
}#end def cor.dist
	
library(gplots)
fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

param.table = read.table("parameters.txt", header=T, sep="\t")
comp.name=as.character(param.table$Value[param.table$Parameter == "comp_name"])
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
min.expression = as.numeric(as.character(param.table$Value[param.table$Parameter == "fpkm_expression_cutoff"]))
aligned.type = as.character(param.table$Value[param.table$Parameter == "FPKM_norm"])
min.fraction.expressed = as.numeric(as.character(param.table$Value[param.table$Parameter == "minimum_fraction_expressed"]))
fc.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "fold_change_cutoff"]))
cor.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "cor_cutoff"]))
cor.cutoff2 = as.numeric(as.character(param.table$Value[param.table$Parameter == "sec_cor_cutoff"]))
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
interaction.flag[interaction.flag == "none"]="no"
pvalue.method = as.character(param.table$Value[param.table$Parameter == "pvalue_method"])
fdr.method = as.character(param.table$Value[param.table$Parameter == "fdr_method"])
goseq.flag = as.character(param.table$Value[param.table$Parameter == "run_goseq"])
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
counts.file = as.character(param.table$Value[param.table$Parameter == "counts_file"])
aligned.stats.file = as.character(param.table$Value[param.table$Parameter == "aligned_stats_file"])
cluster.distance = as.character(param.table$Value[param.table$Parameter == "cluster_distance"])

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

if((aligned.type == "aligned")|(aligned.type =="quantified")){
	if(aligned.type =="aligned"){
		exonic.stat.table = read.table(aligned.stats.file, header=T, sep="\t")
		aligned.reads = as.numeric(exonic.stat.table$aligned.reads[match(longID, exonic.stat.table$Sample)])
	}else if(aligned.type =="quantified"){
		aligned.reads=apply(counts, 2, sum)
	}else{
		stop("Print RPKM_norm must be either 'aligned' or 'quantified'")
	}#end else

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
}else if(aligned.type == "TMM"){
	library(edgeR)
	
	colnames(counts)=sample.label
	y = DGEList(counts=counts, genes=genes)
	y = calcNormFactors(y, method="TMM")
	
	RPKM = round(log2(rpkm(y, gene.length = 1000 * as.numeric(gene.length.kb))+min.expression), digits=2)
}else{
	stop("Print RPKM_norm must be either 'aligned', 'quantified', or 'TMM'")
}#end else

RPKM = RPKM[,match(sample.label, colnames(RPKM))]

print(dim(counts))
expressed.sample.count = apply(RPKM, 1, count.defined.values, expr.cutoff = log2(min.expression))
counts.table = counts.table[expressed.sample.count >= round(min.fraction.expressed * ncol(counts)),]
counts = counts[expressed.sample.count >= round(min.fraction.expressed * ncol(counts)),]
RPKM = RPKM[expressed.sample.count >= round(min.fraction.expressed * ncol(counts)),]
print(dim(counts))

genes = counts.table$symbol
gene.length.kb = as.numeric(counts.table$length.kb)
colnames(counts) = as.character(sample.label)
rownames(counts) = as.character(genes)

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
if(length(groupIDs) == 1){
	average.rpkm = t(average.rpkm)
} else {
	average.rpkm = average.rpkm
}
colnames(average.rpkm) = paste("avg.log2.rpkm", sub("-",".",groupIDs), sep=".")

#remove undefined group IDs (so, you can visualize samples that aren't in your comparison)
if(length(deg.groups) == 1){
	var1 = sample.description.table[,deg.groups]
	deg.counts = counts[,!is.na(var1)]
	deg.RPKM = RPKM[,!is.na(var1)]
	var1 = var1[!is.na(var1)]
	if (trt.group != "continuous"){
		var1 = as.factor(as.character(var1[!is.na(var1)]))
	}
} else if (length(deg.groups) == 2){
	if(interaction.flag == "filter-overlap"){
		var1 = sample.description.table[,deg.groups[1]]
		var2 = sample.description.table[,deg.groups[2]]
		deg.samples = !is.na(var1)&!is.na(var2)
		deg.counts = counts[,deg.samples]
		var1 = var1[deg.samples]
		var2 = var2[deg.samples]
			
		prim.deg.grp = var1[var2 == trt.group2]
		prim.counts = counts[,var2 == trt.group2]
		prim.RPKM = RPKM[,var2 == trt.group2]

		sec.deg.grp = var1[var2 != trt.group2]
		sec.counts = counts[,var2 != trt.group2]
		sec.RPKM = RPKM[,var2 != trt.group2]

	} else{
		var1 = sample.description.table[,deg.groups[1]]
		var2 = sample.description.table[,deg.groups[2]]
		deg.samples = !is.na(var1)&!is.na(var2)
		deg.counts = counts[,deg.samples]
		deg.RPKM = RPKM[,deg.samples]
		var1 = var1[deg.samples]
		if (trt.group != "continuous"){
			var1 = as.factor(as.character(var1[!is.na(var1)]))
		}
		var2 = var2[deg.samples]
		if (trt.group2 != "continuous"){
			var2 = as.factor(as.character(var2[!is.na(var2)]))
		}
	}
} else {
	stop("Code currently doesn't support more than 2 group model for DEG (with or without interaction)")
}

if(length(deg.groups) == 1){
	print("Averaging Expression for One Variable (for deg.groups)")
	contrast.grp = var1
} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
	print("Averaging Expression for First Variable (for deg.groups)")
	contrast.grp = var1
} else if (length(deg.groups) == 2){
	print("Averaging Expression for Interaction Variable (for deg.groups)")
	contrast.grp = paste(var1,var2,sep=":")
} else {
	stop("Code only compatible with 2 variables (with or without a 3rd interaction variable")
}

if (trt.group == "continuous"){
	contrast.grp = as.numeric(contrast.grp)
	
	gene.cor = apply(deg.RPKM, 1, calc.gene.cor, indep.var=contrast.grp)
	
	fc.table = data.frame(cor=gene.cor)
} else {
	groupIDs = as.character(levels(as.factor(contrast.grp)))
	contrast.rpkm = data.frame(t(apply(deg.RPKM, 1, avgGroupExpression, groups = contrast.grp)))
	colnames(contrast.rpkm) = paste("avg.log2.rpkm", sub("-",".",groupIDs), sep=".")
}#end else

if((interaction.flag == "no") & (trt.group != "continuous")){
	print("Calculating fold-change for primary variable")
	trt.expr = contrast.rpkm[,paste("avg.log2.rpkm", sub("-",".",trt.group), sep=".")]
	cntl.expr = contrast.rpkm[,paste("avg.log2.rpkm", sub("-",".",groupIDs[groupIDs != trt.group]), sep=".")]

	log2ratio = round(trt.expr - cntl.expr, digits = 2)
	fc = round(sapply(log2ratio, ratio2fc), digits = 2)
	fc.table = data.frame(log2ratio=log2ratio, fold.change=fc)
} else if ((interaction.flag == "model")|(interaction.flag == "filter-overlap")){
	if ((trt.group == "continuous")&(trt.group2 == "continuous")){
		print("Calculating  correlation for secondary variable")
		sec.contrast.grp = as.numeric(var2)
		
		gene.cor2 = apply(deg.RPKM, 1, calc.gene.cor, indep.var=sec.contrast.grp)
		
		fc.table = data.frame(prim.cor=gene.cor, sec.cor = gene.cor2)
	} else if (trt.group == "continuous"){
		print("Fold-change / correlation cutoff not used for mixed variable analysis")
		print("NOTE: 'Up-Regulated' R output refers to genes that vary with FDR and p-value cutoffs")
		print("However, fold-change / correlation values for each separate variable are still provided")

		sec.groupIDs = var2
		sec.groups = as.character(levels(as.factor(sec.groupIDs)))
		sec.contrast.rpkm = data.frame(t(apply(deg.RPKM, 1, avgGroupExpression, groups = sec.groupIDs)))
		colnames(sec.contrast.rpkm) = paste("avg.log2.rpkm", sub("-",".",groupIDs), sep=".")
		sec.trt.expr = sec.contrast.rpkm[,paste("avg.log2.rpkm", sub("-",".",trt.group2), sep=".")]
		sec.cntl.expr = sec.contrast.rpkm[,paste("avg.log2.rpkm", sub("-",".",sec.groups[sec.groups != trt.group2]), sep=".")]

		sec.log2ratio = round(sec.trt.expr - sec.cntl.expr, digits = 2)
		sec.fc = round(sapply(sec.log2ratio, ratio2fc), digits = 2)
		
		fc.table = data.frame(prim.cor=gene.cor, sec.fc = sec.fc)
	} else if (trt.group2 == "continuous"){	
		print("Fold-change / correlation cutoff not used for mixed variable analysis")
		print("NOTE: 'Up-Regulated' R output refers to genes that vary with FDR and p-value cutoffs")
		print("However, fold-change / correlation values for each separate variable are still provided")

		prim.groupIDs = var1
		prim.groups = as.character(levels(as.factor(prim.groupIDs)))
		prim.contrast.rpkm = data.frame(t(apply(deg.RPKM, 1, avgGroupExpression, groups = prim.groupIDs)))
		colnames(prim.contrast.rpkm) = paste("avg.log2.rpkm", sub("-",".",prim.groups), sep=".")
		prim.trt = trt.group
		prim.cntl = prim.groups[prim.groups != trt.group]
		prim.trt.expr = prim.contrast.rpkm[,paste("avg.log2.rpkm", sub("-",".",prim.trt), sep=".")]
		prim.cntl.expr = prim.contrast.rpkm[,paste("avg.log2.rpkm", sub("-",".",prim.cntl), sep=".")]

		prim.log2ratio = round(prim.trt.expr - prim.cntl.expr, digits = 2)
		prim.fc = round(sapply(prim.log2ratio, ratio2fc), digits = 2)
		
		sec.contrast.grp = as.numeric(var2)
		gene.cor2 = apply(deg.RPKM, 1, calc.gene.cor, indep.var=sec.contrast.grp)
		
		fc.table = data.frame(prim.fc=prim.fc, sec.cor = gene.cor2)
	} else {
		print("Calculating fold-change table for primary variables (within subsets of secondary variable)")
		prim.groups = paste(var1,var2,sep=":")
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
	}#end else
}else if(trt.group == "continuous"){
	print("Skipping fold-change calculation for continuous variable")
}else{
		stop("interaction must be \"no\", \"model\", or \"filter-overlap\"")
}#end else

rep.check = 1
for (i in 1:length(deg.groups)){
	deg.group = deg.groups[i]
	
	if((i == 1) & (trt.group != "continuous")){
		deg.group.values = as.factor(as.character(deg.meta[,deg.group]))
		min.reps = min(table(deg.group.values))
		if (min.reps < 2){
			rep.check=0
			print("There are not at least 2 samples per-group in order to calculate p-value.")
			print("In the future, please make sure you at least have duplicate samples.")
		}#end if (min.reps < 2)
	} else if ((i == 2) & (trt.group2 != "continuous")){
		deg.group.values = as.factor(as.character(deg.meta[,deg.group]))
		min.reps = min(table(deg.group.values))
		if (min.reps < 2){
			rep.check=0
			print("There are not at least 2 samples per-group in order to calculate p-value.")
			print("In the future, please make sure you at least have duplicate samples.")
		}#end if (min.reps < 2)
	} else if (i > 2){
		stop("Workflow currently doesn't support use of more than 2 variables")
	}
}#end for (deg.group in deg.groups)

if(rep.check == 1){
	#start p-value calculation
	if (pvalue.method == "edgeR"){
		library(edgeR)
		
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("edgeR, Two-Step Analysis")
				
			if (trt.group == "continuous"){
				prim.deg.grp = as.numeric(prim.deg.grp)
			}
			
			y = DGEList(counts=prim.counts, genes=genes)
			if(aligned.type == "TMM"){
				y = calcNormFactors(y, method="TMM")
			}#end if(aligned.type == "TMM")
			y = estimateCommonDisp(y)
			design = model.matrix(~prim.deg.grp)
			fit = glmFit(y, design)
			lrt = glmLRT(fit, coef=2)

			prim.pvalue = lrt$table$PValue

			if (trt.group2 == "continuous"){
				sec.deg.grp = as.numeric(prim.deg.grp)
			}
			
			y = DGEList(counts=sec.counts, genes=genes)
			if(aligned.type == "TMM"){
				y = calcNormFactors(y, method="TMM")
			}#end if(aligned.type == "TMM")
			y = estimateCommonDisp(y)
			design = model.matrix(~sec.edgeR.grp)
			fit = glmFit(y, design)
			lrt = glmLRT(fit, coef=2)

			sec.pvalue = lrt$table$PValue
			
		} else {
			y = DGEList(counts=deg.counts, genes=genes)
			if(aligned.type == "TMM"){
				y = calcNormFactors(y, method="TMM")
			}#end if(aligned.type == "TMM")
			y = estimateCommonDisp(y)
			
			if (length(deg.groups) == 1){
				print("edgeR with 1 variable")
				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				design = model.matrix(~var1)
				fit = glmFit(y, design)
				lrt = glmLRT(fit, coef=2)
				test.pvalue = lrt$table$PValue
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("edgeR with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				} else{
					var1 = as.factor(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				} else{
					var2 = as.factor(var2)
				}
				design = model.matrix(~var1 + var2)
				fit = glmFit(y, design)
				lrt = glmLRT(fit, coef=2)
				test.pvalue = lrt$table$PValue
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("edgeR with 2 variables plus interaction")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				design = model.matrix(~var1*var2 + var1 + var2)
				fit = glmFit(y, design)
				lrt = glmLRT(fit, coef=4)
				test.pvalue = lrt$table$PValue
			}
		}#end else
	}else if (pvalue.method == "edgeR-robust"){
		library(edgeR)
		
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("edgeR-robust, Two-Step Analysis")
				
			if (trt.group == "continuous"){
				prim.deg.grp = as.numeric(prim.deg.grp)
			}
			
			y = DGEList(counts=prim.counts, genes=genes)
			if(aligned.type == "TMM"){
				y = calcNormFactors(y, method="TMM")
			}#end if(aligned.type == "TMM")
			design = model.matrix(~prim.deg.grp)
			y = estimateGLMRobustDisp(y, design)
			#fit = glmQLFit(y, design, robust=TRUE)
			fit = glmFit(y, design)
			lrt = glmLRT(fit, coef=2)

			prim.pvalue = lrt$table$PValue

			if (trt.group2 == "continuous"){
				sec.deg.grp = as.numeric(prim.deg.grp)
			}
			
			y = DGEList(counts=sec.counts, genes=genes)
			if(aligned.type == "TMM"){
				y = calcNormFactors(y, method="TMM")
			}#end if(aligned.type == "TMM")
			y = estimateGLMRobustDisp(y)
			design = model.matrix(~sec.edgeR.grp)
			y = estimateGLMRobustDisp(y, design)
			#fit = glmQLFit(y, design, robust=TRUE)
			fit = glmFit(y, design)
			lrt = glmLRT(fit, coef=2)

			sec.pvalue = lrt$table$PValue
			
		} else {
			y = DGEList(counts=deg.counts, genes=genes)
			if(aligned.type == "TMM"){
				y = calcNormFactors(y, method="TMM")
			}#end if(aligned.type == "TMM")
			
			if (length(deg.groups) == 1){
				print("edgeR-robust with 1 variable")
				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				design = model.matrix(~var1)
				y = estimateGLMRobustDisp(y, design)
				
				#fit = glmQLFit(y, design, robust=TRUE)
				fit = glmFit(y, design)
				lrt = glmLRT(fit, coef=2)
				test.pvalue = lrt$table$PValue
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("edgeR-robust with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				} else{
					var1 = as.factor(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				} else{
					var2 = as.factor(var2)
				}
				design = model.matrix(~var1 + var2)
				y = estimateGLMRobustDisp(y, design)
				
				#fit = glmQLFit(y, design, robust=TRUE)
				fit = glmFit(y, design)
				lrt = glmLRT(fit, coef=2)
				test.pvalue = lrt$table$PValue
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("edgeR-robust with 2 variables plus interaction")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				design = model.matrix(~var1*var2 + var1 + var2)
				y = estimateGLMRobustDisp(y, design)
				#fit = glmQLFit(y, design, robust=TRUE)
				fit = glmFit(y, design)
				lrt = glmLRT(fit, coef=4)
				test.pvalue = lrt$table$PValue
			}
		}#end else
	}else if (pvalue.method == "DESeq"){
		library(DESeq)
		
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("DESeq, Two-Step Analysis")
			
			print("Add Code...")
			stop()
			
		} else {
			if (length(deg.groups) == 1){
				print("DESeq with 1 variable (negative binomial)")
				if (trt.group == "continuous"){
					print("Won't work?")
					stop()
				}

				colData = data.frame(var1=var1)
				rownames(colData) = colnames(deg.counts)
				dds = newCountDataSet(deg.counts, var1)
				#dds = newCountDataSet(deg.counts, colData)
				dds = estimateSizeFactors(dds)
				dds = estimateDispersions(dds)
				png(paste(comp.name,"DESeq_dispersion.png",sep="_"))
				plotDispEsts(dds)
				dev.off()
				cats = levels(var1)
				res = nbinomTest(dds, trt.group, cats[cats != trt.group])
				test.pvalue = res$pval
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("DESeq with 2 variables (GLM)")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				} else{
					var1 = as.factor(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				} else{
					var2 = as.factor(var2)
				}

				colData = data.frame(var1=var1, var2=var2)
				rownames(colData) = colnames(deg.counts)

				dds = newCountDataSet(deg.counts, colData)
				dds = estimateSizeFactors(dds)
				#may need to add pooled-CR, per-condition, or blind
				dds = estimateDispersions(dds, method = "pooled")
				png(paste(comp.name,"DESeq_dispersion.png",sep="_"))
				plotDispEsts(dds)
				dev.off()
				fit1 = fitNbinomGLMs(dds, count ~ var2 + var1)
				fit0 = fitNbinomGLMs(dds, count ~ var2 )
				test.pvalue = nbinomGLMTest( fit1, fit0 )
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("DESeq with 2 variables plus interaction (GLM)")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}

				colData = data.frame(var1=var1, var2=var2)
				dds = newCountDataSet(deg.counts, colData)
				dds = estimateSizeFactors(dds)
				#may need to add pooled-CR, per-condition, or blind
				dds = estimateDispersions(dds, method = "pooled")
				png(paste(comp.name,"DESeq_dispersion.png",sep="_"))
				plotDispEsts(dds)
				dev.off()
				fit1 = fitNbinomGLMs(dds, count ~ var2 + var1 + var1*var2)
				fit0 = fitNbinomGLMs(dds, count ~ var2 + var1)
				test.pvalue = nbinomGLMTest( fit1, fit0 )
				print("test result...")
				stop()
			}
		}#end else
	} else if (pvalue.method == "DESeq2"){
		library(DESeq2)
		
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("DESeq2, Two-Step Analysis")
			
			if (trt.group == "continuous"){
				prim.deg.grp = as.numeric(prim.deg.grp)
			}
			
			colData = data.frame(var1=prim.deg.grp)
			rownames(colData) = colnames(prim.counts)
			dds <- DESeqDataSetFromMatrix(countData = prim.counts,
							colData = colData,
							design = ~ var1)
			dds <- DESeq(dds)
			res <- results(dds)
			prim.pvalue = res$pvalue
			

			if (trt.group2 == "continuous"){
				sec.edgeR.grp = as.numeric(prim.edgeR.grp)
			}
			
			colData = data.frame(var1=sec.edgeR.grp)
			rownames(colData) = colnames(sec.counts)
			dds <- DESeqDataSetFromMatrix(countData = sec.counts,
							colData = colData,
							design = ~ var1)
			dds <- DESeq(dds)
			res <- results(dds)
			sec.pvalue = res$pvalue
			
		} else {
			if (length(deg.groups) == 1){
				print("DESeq2 with 1 variable")
				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				colData = data.frame(var1=var1)
				rownames(colData) = colnames(deg.counts)
				dds <- DESeqDataSetFromMatrix(countData = deg.counts,
								colData = colData,
								design = ~ var1)
				dds <- DESeq(dds)
				res <- results(dds)
				test.pvalue = res$pvalue
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("DESeq2 (LRT) with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				} else{
					var1 = as.factor(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				} else{
					var2 = as.factor(var2)
				}

				colData = data.frame(var1=var1, var2=var2)
				rownames(colData) = colnames(deg.counts)
				dds <- DESeqDataSetFromMatrix(countData = deg.counts,
								colData = colData,
								design = ~ var1 + var2)
				Wald.flag = TRUE
				if (Wald.flag){
					dds <- DESeq(dds)
					if (trt.group == "continuous"){
						res <- results(dds, name = "var1")
					} else{
						other.groups = as.character(levels(as.factor(as.character(var1[var1 != trt.group]))))
						if(length(other.groups) > 1){
							print("DESeq2 Wald-test will look at differences between two groups.")
							print("You can manually switch code to use LRT instead of Wald test (set Wald.flag = FALSE).")
							#The paired sample design in the DESeq2 manual uses the Wald test
							stop("Or, please consider using an interaction model with LRT if your primary variable has more than two groups.")
						}
						res <- results(dds, contrast = c("var1", trt.group, other.groups))
					}
				} else {
					dds <- DESeq(dds, test="LRT", reduced = ~ var2)
					res = results(dds)
					#print(head(res))
				}
				test.pvalue = res$pvalue
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("DESeq2 (LRT) with 2 variables plus interaction")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}

				colData = data.frame(var1=var1, var2=var2)
				rownames(colData) = colnames(deg.counts)
				dds <- DESeqDataSetFromMatrix(countData = deg.counts,
								colData = colData,
								design = ~ var1*var2 + var1 + var2)
				dds <- DESeq(dds, test="LRT", reduced = ~ var1 + var2)
				res <- results(dds)
				test.pvalue = res$pvalue
			}
		}#end else
	} else if (pvalue.method == "limma-voom"){
		library(edgeR)
		library(limma)
		
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("limma-voom, Two-Step Analysis")

			if (trt.group == "continuous"){
				prim.deg.grp = as.numeric(prim.deg.grp)
			}
			
			y <- DGEList(counts=prim.counts, genes=genes)
			png(paste(comp.name,"prim_voom_plot.png",sep="_"))
			design <- model.matrix(~prim.deg.grp)
			v <- voom(y,design,plot=TRUE)
			dev.off()
			fit <- lmFit(v,design)
			fit <- eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			prim.pvalue = pvalue.mat[,2]
			

			if (trt.group2 == "continuous"){
				sec.deg.grp = as.numeric(sec.deg.grp)
			}
			
			y <- DGEList(counts=sec.counts, genes=genes)
			design <- model.matrix(~sec.deg.grp)
			png(paste(comp.name,"sec_voom_plot.png",sep="_"))
			v <- voom(y,design,plot=TRUE)
			dev.off()
			fit <- lmFit(v,design)
			fit <- eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			sec.pvalue = pvalue.mat[,2]
			
		} else {
			y <- DGEList(counts=deg.counts, genes=genes)
			if (length(deg.groups) == 1){
				print("limma-voom with 1 variable")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				design <- model.matrix(~var1)
				png(paste(comp.name,"voom_plot.png",sep="_"))
				v <- voom(y,design,plot=TRUE)
				dev.off()
				fit <- lmFit(v,design)
				fit <- eBayes(fit)
				pvalue.mat = data.frame(fit$p.value)
				test.pvalue = pvalue.mat[,2]
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("limma-voom with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				} else{
					var1 = as.factor(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				} else{
					var2 = as.factor(var2)
				}
				design <- model.matrix(~var1 + var2)
				png(paste(comp.name,"voom_plot.png",sep="_"))
				v <- voom(y,design,plot=TRUE)
				dev.off()
				fit <- lmFit(v,design)
				fit <- eBayes(fit)
				pvalue.mat = data.frame(fit$p.value)
				test.pvalue = pvalue.mat[,2]
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("limma-voom with 2 variables plus interaction")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				design <- model.matrix(~var1*var2 + var1 + var2)
				png(paste(comp.name,"voom_plot.png",sep="_"))
				v <- voom(y,design,plot=TRUE)
				dev.off()
				fit <- lmFit(v,design)
				fit <- eBayes(fit)
				pvalue.mat = data.frame(fit$p.value)
				test.pvalue = pvalue.mat[,4]
			}
		}#end else
	}else if (pvalue.method == "DSS"){
		library(DSS)
		
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("DSS, Two-Step Analysis")
			
			print("Add Code...")
			stop()
			
		} else {
			if (length(deg.groups) == 1){
				print("DSS with 1 variable")
				if (trt.group == "continuous"){
					print("Won't work?")
					stop()
				}

				cats = levels(var1)
				var1=as.character(var1)
				names(var1)=colnames(deg.counts)
				seqData=newSeqCountSet(as.matrix(deg.counts), var1)
				seqData=estNormFactors(seqData)
				seqData=estDispersion(seqData)
				res=waldTest(seqData, cats[cats != trt.group], trt.group)
				
				genes = rownames(deg.counts)
				test.pvalue = res$pval[match(genes,rownames(res))]
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("DSS+edgeR with 2 variables")
				library(edgeR)

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				} else{
					var1 = as.factor(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				} else{
					var2 = as.factor(var2)
				}

				colData = data.frame(var1=var1, var2=var2)
				rownames(colData) = colnames(deg.counts)

				X=model.matrix(~var1+var2, data=colData)
				seqData=newSeqCountSet(deg.counts, as.data.frame(X))
				seqData=estNormFactors(seqData)
				seqData=estDispersion(seqData)
				
				fit.edgeR = glmFit(deg.counts, X, dispersion=dispersion(seqData))
				lrt.edgeR = glmLRT(glmfit=fit.edgeR, coef=2)
				test.pvalue = lrt.edgeR$table$PValue
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("DSS+edgeR with 2 variables plus interaction=")
				library(edgeR)
				
				print("Add Code...")
				stop()
			}
		}#end else
	} else if (pvalue.method == "lm"){
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("RPKM linear regression, Two-Step Analysis")

			if (trt.group == "continuous"){
				prim.deg.grp = as.numeric(prim.deg.grp)
			}
			
			prim.pvalue = apply(prim.RPKM, 1, gene.lm, var1=prim.deg.grp)
			

			if (trt.group2 == "continuous"){
				sec.deg.grp = as.numeric(sec.deg.grp)
			}
			
			sec.pvalue = apply(sec.RPKM, 1, gene.lm, var1=sec.deg.grp)
		} else {
			if (length(deg.groups) == 1){
				print("RPKM linear regression with 1 variable")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				test.pvalue = apply(deg.RPKM, 1, gene.lm, var1=var1)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("RPKM linear regression with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				test.pvalue = apply(deg.RPKM, 1, gene.lm, var1=var1, var2=var2)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("RPKM linear regression with 2 variables plus interaction")
				var3 = as.factor(paste(var1,var2,sep=":"))

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group == "continuous"){
					var2 = as.numeric(var2)
				}
				test.pvalue = apply(deg.RPKM, 1, gene.lm, var1=var3, var2=var1, var3=var2)
			}
		}#end else
	} else if (pvalue.method == "ANOVA"){
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("RPKM ANOVA, Two-Step Analysis")

			if (trt.group == "continuous"){
				prim.deg.grp = as.numeric(prim.deg.grp)
			}
			
			prim.pvalue = apply(prim.RPKM, 1, gene.aov, var1=prim.deg.grp)
			

			if (trt.group2 == "continuous"){
				sec.deg.grp = as.numeric(sec.deg.grp)
			}
			
			sec.pvalue = apply(sec.RPKM, 1, gene.aov, var1=sec.deg.grp)
		} else {
			if (length(deg.groups) == 1){
				print("RPKM ANOVA with 1 variable")
				
				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				test.pvalue = apply(deg.RPKM, 1, gene.aov, var1=var1)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("RPKM ANOVA with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				test.pvalue = apply(deg.RPKM, 1, gene.aov, var1=var1, var2=var2)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("RPKM ANOVA with 2 variables plus interaction")
				var3 = as.factor(paste(var1,var2,sep=":"))

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				test.pvalue = apply(deg.RPKM, 1, gene.aov, var1=var3, var2=var1, var3=var2)
			}
		}#end else
	} else{
		stop("pvalue_method must be \"edgeR\", \"limma-voom\", \"DESeq\", \"DESeq2\", \"DSS\", \"lm\", or \"ANOVA\"")
	}
} else{
	test.pvalue = rep(1,times=length(genes))
	prim.pvalue = rep(1,times=length(genes))
	sec.pvalue = rep(1,times=length(genes))
}#end else

if (trt.group == "continuous"){
	upID = "Increased Expression"
	downID = "Decreased Expression"
} else {
	upID = paste(trt.group," Up",sep="")
	downID = paste(trt.group," Down",sep="")	
}


if (interaction.flag == "no"){
	if (fdr.method == "BH"){
		fdr = p.adjust(test.pvalue, "fdr")
	} else if (fdr.method == "q-value"){
		library(qvalue)
		qobj <- qvalue(p = test.pvalue)
		fdr = qobj$qvalue
		png(paste(comp.name,"_",pvalue.method,"_qvalue_plot.png",sep=""))
		qHist = hist(qobj)
		print(qHist)
		dev.off()
	} else if (fdr.method == "q-lfdr"){
		library(qvalue)
		qobj <- qvalue(p = test.pvalue)
		fdr = qobj$lfdr
		png(paste(comp.name,"_",pvalue.method,"_qvalue_plot.png",sep=""))
		qHist = hist(qobj)
		print(qHist)
		dev.off()
	} else {
		stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
	}
	status = rep("No Change", times=length(fdr))
	if (trt.group == "continuous"){
		status[(gene.cor >= cor.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
		status[(gene.cor <= -cor.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
		pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
	} else{
		status[(fc >= fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
		status[(fc <= -fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
		pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
	}#end else
} else{
	trt.group = prim.trt
	if(interaction.flag == "model"){
		if (fdr.method == "BH"){
			fdr = p.adjust(test.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj <- qvalue(p = test.pvalue)
			fdr = qobj$qvalue
			png(paste(comp.name,"_",pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj <- qvalue(p = test.pvalue)
			fdr = qobj$lfdr
			png(paste(comp.name,"_",pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}
		status = rep("No Change", times=length(fdr))
		if ((trt.group == "continuous")&(trt.group2 == "continuous")){
			status[(gene.cor.int >= cor.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			status[(gene.cor.int <= -cor.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
			pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
		} else if ((trt.group != "continuous")&(trt.group2 != "continuous")){
			status[(overall.fc >= fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			status[(overall.fc <= -fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
			pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
		} else {
			upID = "Variable Expression"
			status[(test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
		}#end else
	} else if (interaction.flag == "filter-overlap"){
		if (fdr.method == "BH"){
			fdr = p.adjust(prim.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj <- qvalue(p = prim.pvalue)
			fdr = qobj$qvalue
			png(paste(comp.name,"_",pvalue.method,"_prim_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj <- qvalue(p = prim.pvalue)
			fdr = qobj$lfdr
			png(paste(comp.name,"_",pvalue.method,"_prim_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}
		pass1.status = rep("No Change", times=length(fdr))
		if (trt.group == "continuous"){
			pass1.status[(gene.cor >= cor.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Up",sep="")
			pass1.status[(gene.cor <= -cor.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Down",sep="")
		} else{
			pass1.status[(prim.fc >= fc.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Up",sep="")
			pass1.status[(prim.fc <= -fc.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Down",sep="")
		}#end else

		print(paste("Primary Up-Regulated: ",length(pass1.status[pass1.status == paste(trt.group," Up",sep="")]),sep=""))
		print(paste("Primary Down-Regulated: ",length(pass1.status[pass1.status == paste(trt.group," Down",sep="")]),sep=""))

		if (fdr.method == "BH"){
			sec.fdr = p.adjust(sec.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj <- qvalue(p = sec.pvalue)
			sec.fdr = qobj$qvalue
			png(paste(comp.name,"_",pvalue.method,"_sec_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj <- qvalue(p = sec.pvalue)
			sec.fdr = qobj$lfdr
			png(paste(comp.name,"_",pvalue.method,"_sec_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}		

		pass2.status = rep("No Change", times=length(fdr))
		if (trt.group2 == "continuous"){
			pass2.status[(gene.cor2 >= cor.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Up",sep="")
			pass2.status[(gene.cor2 <= -cor.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Down",sep="")
		} else{
			pass2.status[(sec.fc >= fc.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Up",sep="")
			pass2.status[(sec.fc <= -fc.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Down",sep="")
		}#end else

		print(paste("Secondary Up-Regulated: ",length(pass2.status[pass2.status == paste(trt.group," Up",sep="")]),sep=""))
		print(paste("Secondary Down-Regulated: ",length(pass2.status[pass2.status == paste(trt.group," Down",sep="")]),sep=""))
			
		pvalue.table = data.frame(prim.pvalue = prim.pvalue, prim.FDR = fdr,
										sec.pvalue=sec.pvalue, sec.fdr=sec.fdr)
			
		status = rep("No Change", times=length(fdr))
		status[(pass1.status == paste(trt.group," Up",sep="")) & (pass2.status == "No Change")] = upID
		status[(pass1.status == paste(trt.group," Down",sep="")) & (pass2.status == "No Change")] = downID
	} else{
		stop("interaction must be \"no\", \"model\", or \"filter-overlap\"")
	}#end else
}#end else

print(paste("Up-Regulated: ",length(status[status == upID]),sep=""))
print(paste("Down-Regulated: ",length(status[status == downID]),sep=""))

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

if(length(deg.genes) > 1){
	if (cluster.distance == "Pearson_Dissimilarity"){
		print("Using Pearson Dissimilarity as Distance in Heatmap...")
		dist.fun = cor.dist
	}else{
		dist.fun=dist
	}
	
	if(length(plot.groups) > 1){
		source("heatmap.3.R")
		if((trt.group != "continuous")&(trt.group2 != "continuous")){
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
		}else if((trt.group == "continuous")&(trt.group2 == "continuous")){
			stop("Add code for two continuous variables")
		}else if(trt.group == "continuous"){
			grp1 = as.numeric(sample.description.table[,plot.groups[1]])
			grp2 = as.character(sample.description.table[,plot.groups[2]])
		
			labelColors1 = rep("black",times=length(sample.label))
			library("RColorBrewer")
			continuous.color.breaks = 10
				
			plot.var = as.numeric(grp1)
			plot.var.min = min(plot.var, na.rm=T)
			plot.var.max = max(plot.var, na.rm=T)
				
			plot.var.range = plot.var.max - plot.var.min
			plot.var.interval = plot.var.range / continuous.color.breaks
				
			color.range = colorRampPalette(c("green","black","orange"))(n = continuous.color.breaks)
			plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
			for (j in 1:continuous.color.breaks){
				#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
				labelColors1[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
			}#end for (j in 1:continuous.color.breaks)
			
			group.levels = c(levels(as.factor(grp2)))
			color.palette <- fixed.color.palatte[3:(2+length(group.levels))]
			labelColors2 = rep("black",times=length(sample.label))
			for (i in 1:length(group.levels)){
				labelColors2[grp2 == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
		}else{
			grp1 = as.character(sample.description.table[,plot.groups[1]])
			grp2 = as.numeric(sample.description.table[,plot.groups[2]])

			group.levels = c(levels(as.factor(grp1)))
			color.palette <- fixed.color.palatte[1:(length(group.levels))]
			labelColors1 = rep("black",times=length(sample.label))
			for (i in 1:length(group.levels)){
				labelColors1[grp1 == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
			
			labelColors2 = rep("black",times=length(sample.label))
			library("RColorBrewer")
			continuous.color.breaks = 10
				
			plot.var = as.numeric(grp2)
			plot.var.min = min(plot.var, na.rm=T)
			plot.var.max = max(plot.var, na.rm=T)
				
			plot.var.range = plot.var.max - plot.var.min
			plot.var.interval = plot.var.range / continuous.color.breaks
				
			color.range = colorRampPalette(c("purple","black","cyan"))(n = continuous.color.breaks)
			plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
			for (j in 1:continuous.color.breaks){
				#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
				labelColors2[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
			}#end for (j in 1:continuous.color.breaks)
		}
		
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
		heatmap.3(std.expr,   distfun = dist.fun, hclustfun = hclust,
			  col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
			RowSideColors=row_annotation, trace="none", margins = c(10,15),RowSideColorsSize=4, dendrogram="both")
		if((trt.group != "continuous")&(trt.group2 != "continuous")){
					legend("topright", legend=group.levels,
							col=color.palette,
							pch=15, cex=0.7)
		}else if((trt.group == "continuous")&(trt.group2 == "continuous")){
			stop("Add code for two continuous variables")
		}else if(trt.group == "continuous"){
			legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
			legend("topright", legend=group.levels, col=color.palette, pch=15, cex=0.7)
		}else{
			legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
			legend("topright", legend=group.levels, col=color.palette, pch=15, cex=0.7)
		}
		dev.off()
			
		if(interaction.flag != "no"){
			temp.fc.table = as.matrix(fc.table)
			if (((trt.group == "continuous") & (trt.group2 == "continuous")) | ((trt.group != "continuous") & (trt.group2 != "continuous"))){
				temp.fc.table = temp.fc.table[,-ncol(temp.fc.table)]
			}
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
			heatmap.2(temp.fc.table,  distfun = dist.fun, hclustfun = hclust,
				  col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
				  trace="none", margins = c(20,5), cexCol=1.5)
			dev.off()
		}#end if(interaction.flag != "no")
		
	} else {
		labelColors = rep("black",times=length(sample.label))
		if(trt.group == "continuous"){
			library("RColorBrewer")
			continuous.color.breaks = 10
			
			plot.var = as.numeric(sample.description.table[,plot.groups])
			plot.var.min = min(plot.var, na.rm=T)
			plot.var.max = max(plot.var, na.rm=T)
			
			plot.var.range = plot.var.max - plot.var.min
			plot.var.interval = plot.var.range / continuous.color.breaks
			
			color.range = colorRampPalette(c("green","black","orange"))(n = continuous.color.breaks)
			plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
			for (j in 1:continuous.color.breaks){
				#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
				labelColors[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
			}#end for (j in 1:continuous.color.breaks)
		}else{
			group.levels = levels(as.factor(sample.description.table[,plot.groups]))
			color.palette = fixed.color.palatte[1:length(group.levels)]
			for (i in 1:length(group.levels)){
				labelColors[grp == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
		}

		std.expr = apply(temp.rpkm, 1, standardize.arr)
		if(length(deg.genes) < 25){
			colnames(std.expr) = deg.genes
		} else {
			colnames(std.expr) = rep("", length(deg.genes))
		}
		rownames(std.expr) = sample.label
			
		heatmap.file = paste(comp.name,"_",pvalue.method,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
		heatmap.file = gsub(":",".",heatmap.file)
		png(file = heatmap.file)
		heatmap.2(std.expr, distfun = dist.fun, hclustfun = hclust,
			  col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
			 RowSideColors=labelColors, trace="none", margins = c(5,15))

		if(trt.group == "continuous"){
			legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
		}else{
			legend("topright", group.levels, col=color.palette, pch=15)
		}
		dev.off()
	}#end else
}#end if(length(deg.genes) > 1)

#goseq
if (goseq.flag == "yes"){
	library(goseq)
	
	deg = as.integer(status == upID)
	deg = tapply(deg, as.character(genes), sum)
	deg[deg >= 1] = 1
	gene.symbol = as.character(levels(as.factor(as.character(genes))))
	names(deg)=gene.symbol
	
	bias.file = paste(comp.name,"_goseq_up_bias.png",sep="")
	png(bias.file)
	pwf=nullp(deg,genome,"geneSymbol")
	GO.wall=goseq(pwf,genome,"geneSymbol")
	GO.wall = data.frame(GO.wall, overenrichment.fdr=p.adjust(GO.wall$over_represented_pvalue,"fdr"))

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
	go.file = paste(comp.name,"_goseq_UP.txt",sep="")
	write.table(go.table, file=go.file, row.names=F, quote=F, sep="\t")
	
	go.file = paste(user.folder,"/GO/",comp.name,"_goseq_UP.txt",sep="")
	write.table(go.table, file=go.file, row.names=F, quote=F, sep="\t")

	deg = as.integer(status == downID)
	deg = tapply(deg, as.character(genes), sum)
	deg[deg >= 1] = 1
	gene.symbol = as.character(levels(as.factor(as.character(genes))))
	names(deg)=gene.symbol

	bias.file = paste(comp.name,"_goseq_down_bias.png",sep="")
	png(bias.file)
	pwf=nullp(deg,genome,"geneSymbol")
	GO.wall=goseq(pwf,genome,"geneSymbol")
	GO.wall = data.frame(GO.wall, overenrichment.fdr=p.adjust(GO.wall$over_represented_pvalue, "fdr"))

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
	go.file = paste(comp.name,"_goseq_DOWN.txt",sep="")
	write.table(go.table, file=go.file, row.names=F, quote=F, sep="\t")
	
	go.file = paste(user.folder,"/GO/",comp.name,"_goseq_DOWN.txt",sep="")
	write.table(go.table, file=go.file, row.names=F, quote=F, sep="\t")
}#end if (goseq.flag == "yes")
