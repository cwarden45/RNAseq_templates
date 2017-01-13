# RNAseq_templates
template scripts for RNA-Seq analysis

Genome_Ref_Code is used to prepare annotations for TopHat_Workflow

Transcriptome_Ref_Code is used to prepare sequences / annotations for Transcriptome_Workflow

Splicing_Workflow and Mutation_Workflow use genomic alignments

Splicing_Workflow uses TopHat alignment from TopHat_Workflow, but provides new quantification for exons / splice junctions

Mutation_Workflow uses 2-pass STAR alignment for GATK variant calling

Please see the README file in each section.
