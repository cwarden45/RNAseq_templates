### Order to Run Scripts ###

1) `cluster_pe_STAR_alignment.py` (or `cluster_se_STAR_alignment.py` + `STAR_post_processing.py`)

2) `GATK_variant_calls.py`

3) `filter_variants.py`

4) `annotate_variants.py`

5) `variant_summary.R`


### Dependencies (some optional) ###

**STAR**: https://github.com/alexdobin/STAR

**Picard**: https://broadinstitute.github.io/picard/

**GATK**: https://software.broadinstitute.org/gatk/

**VarScan**: http://varscan.sourceforge.net/

**bedtools**: http://bedtools.readthedocs.io/en/latest/

**ANNOVAR**: http://annovar.openbioinformatics.org/en/latest/

**UCSC Table Browser (used to download RepeatMasker for any genome)**: https://genome.ucsc.edu/cgi-bin/hgTables

-also used to download human GWAS Catalog annotations

**ORegAnno (for human genome; includes miRNA targets, and you might get some other regulatory mutations with variable UTRs)**: http://www.oreganno.org/dump/
-link to temporary website annotation

**REDIportal (hg19 RNA-editing events)**: http://srv00.recas.ba.infn.it/atlas/


### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name	| Name of differential expression comparison (used to name output file)
|Result_Folder|Path to output folder for selected, final results|
|Alignment_Folder|Path to TopHat Alignments|
|Java_Mem|Java memory allocation for Picard/GATK|
|Threads|Number of Threads for STAR|
|strand|Library type.  Can be *no*, *yes*, or *reverse* (*reverse* typically used for Illumina stranded libraries)|
|Cluster_Email|If running STAR on cluster, e-mail address for updates|
|Max_Alt_Alleles|Maximum alternative alleles for GATK HaplotypeCaller (default = 6)|
|Path_to_STAR|Full path to STAR binary|
|ANNOVAR_Path|Path to ANNOVAR script folder|
|genome|Name of genome build|
|RepeatMasker_BED|BED file of repetitive elements to exclude variants|
|RNA_Editing_BED|BED file of RNA editing elements to exclude variants|

