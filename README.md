### \~Temporary Note\~ ###
**I apologize for the confusion, but I would like to emphasize that these are called “templates” because I almost always have to modify the code for each project (beyond the parameter files), meaning they will be more difficult for other people to use in the same way.**  This was unfortunately not immediately clear to me when I created the templates.

I also believe that the process of writing the scripts for analysis (such as the templates) is very important for the learning process, and it is very important that you understand all the steps for analysis before presenting them in a paper.

So, I will post an update when more specific guidance / suggestions can be provided.  **In the meantime, I apologize, but I cannot provide user support for the templates.**

# RNAseq_templates
template scripts for RNA-Seq analysis

Genome_Ref_Code is used to prepare annotations for TopHat_Workflow

Transcriptome_Ref_Code is used to prepare sequences / annotations for Transcriptome_Workflow

Splicing_Workflow and Mutation_Workflow use genomic alignments

Splicing_Workflow uses TopHat alignment from TopHat_Workflow, but provides new quantification for exons / splice junctions

Mutation_Workflow uses 2-pass STAR alignment for GATK variant calling

Please see the README file in each section.
