import os

mapFile = "/path/to/RSEM_Bowtie_index/transcript_to_gene.txt"
refFa = "/path/to/RSEM_Bowtie_index/gencode.v24.transcriptID.fa"
rsemIndex = "/path/to/RSEM_Bowtie_index/Bowtie_gencode_v24"
command = "/opt/RSEM/rsem-prepare-reference --bowtie --bowtie-path /opt/bowtie --transcript-to-gene-map " + mapFile + " " + refFa + " " + rsemIndex
os.system(command)
