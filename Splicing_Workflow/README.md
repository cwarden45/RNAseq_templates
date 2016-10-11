### Order to Run Scripts ###

0) It is recommended that you filter the GENCODE .gtf file to match the gene annotation set using `filter_gtf.py`

1) se_alignment.py or cluster_se_alignment.py (from [TopHat_Workflow](https://github.com/cwardn45/RNAseq_templates/tree/master/TopHat_Workflow))

2) run_QoRTs.py

3) junction_hist.R

4) run_JunctionSeq.R (plots created with FDR cutoff, but you can usually detect something with FDR < 0.05)

5) create_stat_table.R (this is the step where you would run goseq)

### Dependencies (some optional) ###

Most Python scripts can be run using this [Docker image](https://hub.docker.com/r/cwarden45/rnaseq-dependencies/)

*Alignment*

TopHat: https://ccb.jhu.edu/software/tophat/tutorial.shtml


*Read / Junction Counting*

QoRTs: http://hartleys.github.io/QoRTs/index.html

*Differential Splicing*

JunctionSeq: https://bioconductor.org/packages/release/bioc/html/JunctionSeq.html

SGSeq: https://bioconductor.org/packages/release/bioc/html/SGSeq.html

MATS: http://rnaseq-mats.sourceforge.net/


*Functional Enrichment*

goseq: http://bioconductor.org/packages/release/bioc/html/goseq.html

### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name	| Name of differential expression comparison (used to name output file)
|plot_groups | Names of columns in *sample_description_file* to be plotted in QC plots.  Use commas to plot multiple groups|
|deg_groups|Names of columns in *sample_description_file* to be used for differential splicing analysis.  Use commas to include multiple variables (for multivariate model).  However, primary variable will be called "condition" for JunctionSeq, which is only designed for categorical variables.|
|Raw_Code_PC|Path to output folder for most results|
|Result_Folder|Path to output folder for selected, final results|
|Alignment_Folder|Path to TopHat Alignments|
|QoRTs_Count_Folder| Per-Sample QoRTs QC/counts|
|QoRTs_Merged_Folder| Folder with QoRTs merged counts formatted for JunctionSeq|
|fdr_method|Method to Calculate FDR.  Can be *BH* (Benjamini and Hochberg),*q-value*, or *q-lfdr*|
|genome|Name of genome build|
|Java_Mem|Java memory allocation for QoRTs|
|Threads|Number of Threads for JunctionSeq Analysis (although it will switch to single-core for Windows users)|
|GENCODE_GTF|Path to GENCODE GTF file for QoRTs quantification.  Recommend using filtered gene annotations from [Genome_Ref_Code](https://github.com/cwarden45/RNAseq_templates/tree/master/Genome_Ref_Code)|
|GENCODE_Gene_Info|Path to GENCODE gene annotation file (created using [Genome_Ref_Code](https://github.com/cwarden45/RNAseq_templates/tree/master/Genome_Ref_Code))|
|sample_description_file|Name of Sample Description File.  "unique.ID" and "sample.ID" must match the QoRTs counts folder names, but you can use *userID* to change sample labels in result|
|pvalue_cutoff|Maximum p-value to consider a gene differenitally spliced|
|fdr_cutoff|Maximum FDR to consider a gene differentially spliced|
|strand|Library type.  Can be *no*, *yes*, or *reverse* (*reverse* typically used for Illumina stranded libraries)|
|pairing|*SE* for single-end reads, *PE* for paired-end reads|
|run_goseq| Run goseq? *yes* or *no*|
