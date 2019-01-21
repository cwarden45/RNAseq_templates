### Acknowledgements ###

Code developed when providing analysis support from the City of Hope Integrative Genomics Core with requests / suggestions from *Ya-Huei Kuo* (1st successful testing with QoRTs / JunctionSeq), *[Internal Collaborator B] (testing / revision with QoRTs/JunctionSeq/SGSeq; currently, did not want to be explicitly acknowledged here)*, *Joo Song + Alex Herrera / Vanessa Jonsson* (testing with QoRTs/JunctionSeq/SGSeq), and *Marcia Miller* (validation with QoRTs junctions; for project with Ronald Goto and Jibin Zhang).  I am very grateful for their discussions.

In general, some projects have been analyzed or re-analyzed with different staff members (sometimes in different labs/institutes, sometimes within the IGC).  So, the use of data with modified versions of these templates may not be what is eventually used in the associated publication.


### \~Temporary Note\~ ###
**I apologize for the confusion, but I would like to emphasize that these are called “templates” because I almost always have to modify the code for each project (beyond the parameter files), meaning they will be more difficult for other people to use in the same way.**  This was unfortunately not immediately clear to me when I created the templates.

I also believe that the process of writing the scripts for analysis (such as the templates) is very important for the learning process, and it is very important that you understand all the steps for analysis before presenting them in a paper.

So, I will post an update when more specific guidance / suggestions can be provided.  **In the meantime, I apologize, but I cannot provide user support for the templates.**


### Order to Run Scripts ###

0) It is recommended that you filter the GENCODE .gtf file to match the gene annotation set using `filter_gtf.py`, which also creates the necessary gene info file (which is not identical to what is created from .fasta file in transcriptome workflow)

1) `cluster_se_alignment.py` or `cluster_pe_alignment.py` (or se_alignment.py) (from [TopHat_Workflow](https://github.com/cwardn45/RNAseq_templates/tree/master/TopHat_Workflow))

2) `run_QoRTs.py`

3) `junction_hist.R`

4) `run_JunctionSeq.R` (plots created with FDR cutoff of 0.001, but the R image is saved if you want to make new plots)

5) `create_stat_table.R` (this is the step where you would run goseq)

6) `CPM_feature_table.R` (CPM values for exons and junctions; assumes there is a column called "aligned.reads" in the sample description file)

7) `run_SGSeq.R` (optional: used to predict novel exons - uses different annotation database and is computationally intenstive, so try to make your number of differentially spliced genes below 500)

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

*Gene Fusions*

JAFFA: https://github.com/Oshlack/JAFFA/wiki

*Functional Enrichment*

goseq: http://bioconductor.org/packages/release/bioc/html/goseq.html

### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name	| Name of differential expression comparison (used to name output file)
|plot_groups | Names of columns in *sample_description_file* to be plotted in QC plots.  Use commas to plot multiple groups|
|plot_type | Are lables for Junction QC plots "discrete" or "continous"?  Use commas to describe multiple variables.  If continous, orange=high, green=low|
|deg_groups|Names of columns in *sample_description_file* to be used for differential splicing analysis.  Use commas to include multiple variables (for multivariate model).  However, primary variable will be called "condition" for JunctionSeq, which is only designed for categorical variables.|
|Raw_Code_PC|Path to output folder for most results|
|Result_Folder|Path to output folder for selected, final results|
|Alignment_Folder|Path to TopHat Alignments|
|QoRTs_Count_Folder| Per-Sample QoRTs QC/counts|
|QoRTs_Merged_Folder| Folder with QoRTs merged counts formatted for JunctionSeq|
|genome|Name of genome build|
|Java_Mem|Java memory allocation for QoRTs|
|Threads|Number of Threads for JunctionSeq Analysis (although it will switch to single-core for Windows users)|
|GENCODE_GTF|Path to GENCODE GTF file for QoRTs quantification.  Recommend using filtered gene annotations from [Genome_Ref_Code](https://github.com/cwarden45/RNAseq_templates/tree/master/Genome_Ref_Code)|
|GENCODE_Gene_Info|Path to GENCODE gene annotation file (created using [Genome_Ref_Code](https://github.com/cwarden45/RNAseq_templates/tree/master/Genome_Ref_Code))|
|sample_description_file|Name of Sample Description File.  "unique.ID" and "sample.ID" must match the QoRTs counts folder names, but you can use *userID* to change sample labels in result|
|treatment_group|Treatment group for primary variable (used when calculating fold-change values for exons from normalized counts).  JunctionSeq not designed to work with continuous variable.|
|fold_change_cutoff|Minimum fold-change difference to consider an exon differentially expressed|
|pvalue_cutoff|Maximum p-value to consider a gene differenitally spliced|
|fdr_cutoff|Maximum FDR to consider a gene differentially spliced|
|strand|Library type.  Can be *no*, *yes*, or *reverse* (*reverse* typically used for Illumina stranded libraries)|
|pairing|*SE* for single-end reads, *PE* for paired-end reads|
|run_goseq| Run goseq? *yes* or *no*|
