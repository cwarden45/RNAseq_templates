### Acknowledgements ###
Code developed when providing analysis support from the City of Hope Integrative Genomics Core with requests / suggestions from *[Internal Collaborator B]* (*currently, did not want to be explicitly acknowledged here*), *Yapeng Su* (Caltech; also, paired Xenograft samples with Exomes for collaboration with *Wei Wei from UCLA/Crump Institute*), *Joo Song + Alex Herrera / Vanessa Jonsson* (Access RNA-Seq, human), and *David Ann* (mouse samples).  I am very grateful for their discussions.

In general, some projects have been analyzed or re-analyzed with different staff members (sometimes in different labs/institutes, sometimes within the IGC).  So, the use of data with modified versions of these templates may not be what is eventually used in the associated publication.


### \~Temporary Note\~ ###
**I apologize for the confusion, but I would like to emphasize that these are called “templates” because I almost always have to modify the code for each project (beyond the parameter files), meaning they will be more difficult for other people to use in the same way.**  This was unfortunately not immediately clear to me when I created the templates.

I also believe that the process of writing the scripts for analysis (such as the templates) is very important for the learning process, and it is very important that you understand all the steps for analysis before presenting them in a paper.

So, I will post an update when more specific guidance / suggestions can be provided.  **In the meantime, I apologize, but I cannot provide user support for the templates.**


### Order to Run Scripts ###

1) `cluster_pe_STAR_alignment.py`

2) `cluster_GATK_joint_variant_calls.py` or `GATK_joint_variant_calls.py`

3) `filter_variants.py` or `filter_variants_joint.py`

4) `annotate_variants.py` or `annotate_variants_joint.py` or `cluster_annotate_variants_joint.py` (with `bed_to_ANNOVAR.py`)

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

