### Order to Run Scripts ###

0) It is recommended that you filter the GENCODE .gtf file to match the gene annotation set using `filter_gtf.py`

1) se_alignment.py or cluster_se_alignment.py (from [TopHat_Workflow](https://github.com/cwardn45/RNAseq_templates/tree/master/TopHat_Workflow))

2) run_QoRTs.py

3) run_JunctionSeq.R

### Dependencies (some optional) ###

Most Python scripts can be run using this [Docker image](https://hub.docker.com/r/cwarden45/rnaseq-dependencies/)

*Alignment*

TopHat: https://ccb.jhu.edu/software/tophat/tutorial.shtml


*Read / Junction Counting*

QoRTs: http://hartleys.github.io/QoRTs/index.html

*Differential Splicing*

JunctionSeq: https://bioconductor.org/packages/release/bioc/html/JunctionSeq.html
MATS: http://rnaseq-mats.sourceforge.net/

### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name	| Name of differential expression comparison (used to name output file)
|plot_groups | Names of columns in *sample_description_file* to be plotted in QC and differential expression plots.  Use commas to plot multiple groups|
|deg_groups|Names of columns in *sample_description_file* to be plotted in QC and differential expression plots.  Use commas to include multiple variables (for multivariate model or gene list filtering)|
|treatment_group|Treatment group for primary variable; enter *continuous* for a continuous variable and a correlation will be provided instead of a fold-change value.|
|Raw_Code_PC|Path to output folder for most results|
|Result_Folder|Path to output folder for selected, final results|
|Alignment_Folder|Path to TopHat Alignments|
|fdr_method|Method to Calculate FDR.  Can be *BH* (Benjamini and Hochberg),*q-value*, or *q-lfdr*|
|genome|Name of genome build|
|Java_Mem|Java memory allocation for QoRTs|
|Threads|Number of Threads for JunctionSeq Analysis|
|GENCODE_GTF|Path to GENCODE GTF file for QoRTs quantification.  Recommend using filtered gene annotations from [Genome_Ref_Code](https://github.com/cwarden45/RNAseq_templates/tree/master/Genome_Ref_Code)|
|GENCODE_Gene_Info|Path to GENCODE gene annotation file (created using [Genome_Ref_Code](https://github.com/cwarden45/RNAseq_templates/tree/master/Genome_Ref_Code))|
|sample_description_file|Name of Sample Description File|
|pvalue_cutoff|Maximum p-value to consider a gene differenitally expressed|
|fdr_cutoff|Maximum FDR to consider a gene differentially expressed|
|sec_fold_change_cutoff|If comparing two gene lists, fold-change threshold for list you want to filter out|
|sec_cor_cutoff|If comparing two gene lists and using a continuous variable, minimum absolute correlation to consider a gene differentially expressed|
|sec_pvalue_cutoff|If comparing two gene lists, p-value threshold for list you want to filter out|
|sec_fdr_cutoff|If comparing two gene lists, FDR threshold for list you want to filter out|
|strand|Library type.  Can be *no*, *yes*, or *reverse* (*reverse* typically used for Illumina stranded libraries)|
|pairing|*SE* for single-end reads, *PE* for paired-end reads|
|run_goseq| Run goseq?  It is useful to leave this as 'no' initially, and then switch to 'yes' after optimizing differential expression parameters|
|interaction| Method for comparing an interaction of two variables.  Can be *model*, *filter-overlap*, or *no*|
|secondary_trt| If comparing two gene lists, this is treatment group for the list that you want to filter out; enter *continuous* for a continuous variable and a correlation will be provided instead of a fold-change value (also converts second variable from factor to numeric, even if interaction is set to *no*)|
