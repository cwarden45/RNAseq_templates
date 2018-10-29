### \~Temporary Note\~ ###
**I apologize for the confusion, but I would like to emphasize that these are called “templates” because I almost always have to modify the code for each project (beyond the parameter files), meaning they will be more difficult for other people to use in the same way.**  This was unfortunately not immediately clear to me when I created the templates.

I also believe that the process of writing the scripts for analysis (such as the templates) is very important for the learning process, and it is very important that you understand all the steps for analysis before presenting them in a paper.

So, I will post an update when more specific guidance / suggestions can be provided.

### Refence Notes for Transcriptome_Workflow ###

**For both RSEM and Sailfish/Salmon Workflows**: Use RSEM_gene_map.py to create a filtered transcriptome from *GENCODE* annotations, and use MiTranscriptome_annotations.py to create annotation file for *MiTranscriptome* sequences.

This also creates gene to transcript mapping file for RSEM and annotation file for all transcriptome-based workflows

NOTE: The MiTranscriptome annotation file has a column called "Ensembl.TranscriptID".  This is not really an Ensembl ID - it is just the transcript name repeated so that the annoation file will work with the GENCODE reference scripts.  This is also why NCBI.GeneID is "NA" for all values.

**Salmon Index Command**: `/opt/salmon/bin/salmon index -t [RSEM_gene_map.py trimmed rsemFa].fa -i salmon_index`
- Recommend using default setting of *--type quasi*
- *-k* must be less than the length of the read and must be odd number (manual recommends less than 1/2 read length).  Using shorter reads, use a lower value of k (such as 19).

### Annotation Resources ###

    GENCODE (scripts use .fasta sequences and annotations in sequence name)

Human:
http://www.gencodegenes.org/releases/reference_releases.html

Mouse:
http://www.gencodegenes.org/mouse_releases/reference_releases.html

    MiTranscriptome

Human:
http://www.mitranscriptome.org/

*use Table Browser in Genome Browser to download sequences for cancer-specific genes for reference build of interest
