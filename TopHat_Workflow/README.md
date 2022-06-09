### Acknowledgements ###

Code developed when providing analysis support from the City of Hope Integrative Genomics Core with requests / suggestions from **Sanjay Awasthi** (formerly at City of Hope, currently at [TTUHSC](http://www.ttuhsc.edu/about/profile.aspx?eRaiderUserName=sawasthi); includes a collaboration with Sharad Singhal), **[Internal / External Collaborator A]** (previously at City of Hope, quantification test for re-analysis, initially with exons, influencing use of htseq-count in new scripts; *not able to get a response for acknowledgement approval*), **Rebecca Deegan** ([UNMC](https://www.unmc.edu/biochemistry/faculty/oberley-deeggan.html), first conversion from the IGC RPKB scripts; Arpita Chatterjee was the lab contact), **Jeremy Jones** (early parameterization), **Marc Baum** ([OCIS](http://www.oak-crest.org/); early RPKB and modified template re-analysis),  **Shiuan Chen** (application to 10+ runs; including [Petrossian et al. 2018](https://www.ncbi.nlm.nih.gov/pubmed/29963233); also includes some samples that are a collaboration with *Susan Neuhausen*, such as those in [Kanaya et al. 2019](https://www.ncbi.nlm.nih.gov/pubmed/30796839)), **Debbie Thurmond**  (application to 10+ runs, including [Oh et al. 2018](https://www.ncbi.nlm.nih.gov/pubmed/30305365) and [Merz et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35222279/)), **[External Collaborator B]** (application to 10+ runs; *currently, did not want to be explicitly acknowledged here*), and **Yapeng Su** (Caltech, includes collaboration with **Wei Wei from UCLA/Crump Institute**; application to 10+ runs, including content within [Su et al. 2019](https://www.biorxiv.org/content/10.1101/724740v1.abstract) preprint).

*De Novo* assembly (usually of unaligned reads) was developed / tested for **Marc Baum** ([OCIS](http://www.oak-crest.org/)), **Debbie Thurmond** (with Karla Merz and Miwon Ahn), and **[Internal Collaborator B]** (*currently, did not want to be explicitly acknowledged here*).

I would also like to thank the following users that gave us permission to acknowledge these templates (or scripts similar to these templates) were used for analysis of their data: **David Ann**, **Behnam Badie** (with processed data for project with Jacob Berlin in [GSE102294](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102294), as well as other runs), **Yuan Chen**, **Xiao Chen / James Wang** ([StemBios Technologies, Inc.](http://www.stembios.com/)), **Hsun “Teresa” Ku**, **Marcia Miller** (with lab members Ronald Goto and Jibin Zhang), **John Rossi** (Sorah Yoon’s project), **Paul Salvaterra** (including [Ubina et al. 2018][https://www.biorxiv.org/content/early/2018/09/17/419978](https://www.frontiersin.org/articles/10.3389/fnins.2019.01007/full)), **Sharad Singhal** (including [Nagaprashantha et al. 2018]( https://www.ncbi.nlm.nih.gov/pubmed/29719590)), **Joo Song + Alex Herrera / Vanessa Jonsson**, **Zuoming Sun**,  **De-Fu Zeng**, **Jianhua Yu** (Lei Tian’s project), *Yuan Chun Ding* (for limited QC of a joint project with **Susan Neuhausen** and **Sunita Patel**), and **Yiling Hong** ([Western University of Health Sciences](https://www.westernu.edu/veterinary/faculty-spotlight-dr-yiling-hong/), “initial” analysis for Xu Dong’s project; includes brief test of enrichR R package code; referenced in [Hong et al. 2021](https://www.biorxiv.org/content/10.1101/2021.08.06.455467v1.abstract) preprint).  There are also **2** internal collaborators with whom I did not get a clear confirmation for approval.

There was also testing of MiTranscriptome lncRNA annotations for **Shiuan Chen** (and one other lab) for ribosome-depleted RNA-Seq libraries, but those results may not be used in the resulting publication.

Parts of the *exonic_reads_counts.R* script were based upon code for earlier gene quantification (so, I became familiar with the GenomicRanges functions from earlier IGC code).  So, that name reflects the influence on the code that I wrote, but it is not a precise description for it’s current use: quantifying chromosomal aligned read counts (for desired chromosomes).

My understanding is that *Min-Hsuan “Michael” Chen* is also currently using portions of this code (chunks within some of the early steps for htseq-count reformatting, but otherwise re-written from scratch) for IGC analysis of some RNA-Seq samples.

In general, some projects have been analyzed or re-analyzed with different staff members (sometimes in different labs/institutes, sometimes within the IGC).  So, the use of data with modified versions of these templates may not be what is eventually used in the associated publication.

*Xiwei Wu* (the Director of the Integrative Genomics Core) was kind enough to review this acknowledgement, and he approves of the content.


### \~Temporary Note\~ ###
**I apologize for the confusion, but I would like to emphasize that these are called “templates” because I almost always have to modify the code for each project (beyond the parameter files), meaning they will be more difficult for other people to use in the same way.**  This was unfortunately not immediately clear to me when I created the templates.

I also believe that the process of writing the scripts for analysis (such as the templates) is very important for the learning process, and it is very important that you understand all the steps for analysis before presenting them in a paper.

So, I will post an update when more specific guidance / suggestions can be provided.  **In the meantime, I apologize, but I cannot provide user support for the templates.**

**Update (10/11/2019):** It is a work in progress (and doesn't represent everything from my experiences that I want to show), but I have started to try and show the need to have testing for every project using public data [here](https://sourceforge.net/projects/rnaseq-deg-methodlimit/) (which also includes a public log / "lab notebook" for this project)

### Order to Run Scripts ###

1) `cluster_se_alignment.py` or `cluster_pe_alignment.py`

2) `extract_total_reads.py` (from TopHat alignment) or `extract_total_reads_FastQC.py` (if using STAR)

3) `exonic_read_counts.R`

4) `cluster_HTseq_counts.py`

4b) Start `create_TDFs.py'

4c) `run_RSeQC.py` + `collect_RSeQC_stats.py` / `housekeeping_RSeQC_read_dist_barplot.py`
--> Need to have sample description (possibly created by `create_sample_description.R`)

5) `HTseq_reformat.R`

6) `qc.R`

7) `DEG_goseq.R`

7b) Can create additional files with `web_enrichment_input.R`

### Dependencies (some optional) ###

Most Python scripts can be run using this [Docker image](https://hub.docker.com/r/cwarden45/rnaseq-dependencies/)

*Alignment*

TopHat2: https://ccb.jhu.edu/software/tophat/tutorial.shtml

Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

 - TopHat2 alignments run against bowtie2 indexed reference, created using `bowtie2-build ref_name.fa ref_name`

GenomicAlignments: https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html

*Read Counting / FPKM*

HTseq: http://www-huber.embl.de/HTSeq/doc/install.html#install

GenomicRanges: https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html

*Differential Expression*

edgeR: https://bioconductor.org/packages/release/bioc/html/edgeR.html

limma-voom: https://bioconductor.org/packages/release/bioc/html/limma.html

DESeq: https://bioconductor.org/packages/release/bioc/html/DESeq.html

DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html

DSS: http://bioconductor.org/packages/release/bioc/html/DSS.html

qvalue: https://bioconductor.org/packages/release/bioc/html/qvalue.html

*Visualization*

gplots: https://cran.r-project.org/web/packages/gplots/index.html

RColorBrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html

heatmap.3: https://github.com/obigriffith/biostar-tutorials/blob/master/Heatmaps/heatmap.3.R

heatmap.3 example: https://www.biostars.org/p/18211/

IGV / igvtools: http://software.broadinstitute.org/software/igv/

*Gene Set Enrichment*

goseq: http://bioconductor.org/packages/release/bioc/html/goseq.html

BD-Func: https://sourceforge.net/projects/bdfunc/

GSEA: http://software.broadinstitute.org/gsea/index.jsp

IPA: https://www.qiagenbioinformatics.com/products/ingenuity-pathway-analysis/

*QC*

RSeQC: http://rseqc.sourceforge.net/

(use with `run_RSeQC.py` and `collect_RSeQC_metrics.py`)

### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name	| Name of differential expression comparison (used to name output file)
|plot_groups | Names of columns in *sample_description_file* to be plotted in QC and differential expression plots.  Use commas to plot multiple groups|
|deg_groups|Names of columns in *sample_description_file* to be plotted in QC and differential expression plots.  Use commas to include multiple variables (for multivariate model or gene list filtering)|
|treatment_group|Treatment group for primary variable; enter *continuous* for a continuous variable and a correlation will be provided instead of a fold-change value.|
|Raw_Code_PC|Path to output folder for most results|
|Result_Folder|Path to output folder for selected, final results|
|Alignment_Folder_MAC|Path to TopHat Alignments|
|Reads_Folder_MAC|Path to Reads for TopHat Alignment|
|Cluster_Email|If running alignment on a cluster, e-mail for notifications|
|pvalue_method|Method to Calculate P-value.  Can be *edgeR*, *edgeR-robust*, *limma-voom*, *DESeq*, *DESeq2*, *DSS*, *lm* (linear regression), or *aov* (ANOVA)|
|fdr_method|Method to Calculate FDR.  Can be *BH* (Benjamini and Hochberg),*q-value*, or *q-lfdr*|
|genome|Name of genome build|
|Bowtie2_Ref| Path to Bowtie2 ref|
|Threads|Number of Threads for TopHat Alignment|
|txGTF_MAC|Path to GTF file for HTseq quantification|
|lncRNA_GTF_MAC|Path to lncRNA GTF file for HTseq quantification|
|lncRNA_ann_name|Game of extra lncRNA database (*GENCODE* or *MiTranscriptome*)|
|RSeQC_bed_MAC|Path to .bed file with housekeeping gene annotations (for RSeQC)|
|HTseq_input_folder|Path to folder containing GTF gene information|
|sample_description_file|Name of Sample Description File|
|total_counts_file|Name of File to Contain Total Read Counts|
|aligned_stats_file|Name of File to Contain Aligned Read Counts|
|cluster_distance| Distance metric for dendrogram.  Can be *Euclidean* or *Pearson_Dissimilarity*|
|fpkm_file|Name of File to Contain log2(FPKM + *rpkm_expression_cutoff*) Expression Values|
|FPKM_norm|How to count number of aligned reads: *aligned*, *quantified*, or *TMM*|
|counts_file|Name of File to Contain Read Counts Per Gene|
|fpkm_expression_cutoff|Rounding Value for FPKM (minimum reliable expression level)|
|minimum_fraction_expressed|Minimum fraction of samples with expression above *fpkm_expression_cutoff*. Filter for differential expression analysis.|
|fold_change_cutoff|Minimum fold-change difference to consider a gene differentially expressed|
|cor_cutoff|If using a continuous variable, minimum absolute correlation to consider a gene differentially expressed|
|pvalue_cutoff|Maximum p-value to consider a gene differenitally expressed|
|fdr_cutoff|Maximum FDR to consider a gene differentially expressed|
|sec_fold_change_cutoff|If comparing two gene lists, fold-change threshold for list you want to filter out|
|sec_cor_cutoff|If comparing two gene lists and using a continuous variable, minimum absolute correlation to consider a gene differentially expressed|
|sec_pvalue_cutoff|If comparing two gene lists, p-value threshold for list you want to filter out|
|sec_fdr_cutoff|If comparing two gene lists, FDR threshold for list you want to filter out|
|min_length_kb|Minimum length of gene to consider for differential expression anlysis (in kilobases)|
|strand|Library type.  Can be *no*, *yes*, or *reverse* (*reverse* typically used for Illumina stranded libraries)|
|run_goseq| Run goseq?  It is useful to leave this as 'no' initially, and then switch to 'yes' after optimizing differential expression parameters|
|interaction| Method for comparing an interaction of two variables.  Can be *model*, *filter-overlap*, or *no*|
|secondary_trt| If comparing two gene lists, this is treatment group for the list that you want to filter out; enter *continuous* for a continuous variable and a correlation will be provided instead of a fold-change value (also converts second variable from factor to numeric, even if interaction is set to *no*)|
|Read_Pairing|Are you using single-end (*SE*, typical for gene expession analysis) or paired-end (*PE*) reads?|
