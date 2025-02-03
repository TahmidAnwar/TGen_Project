# Fragmentation Pattern Analysis of Circulating Tumor DNA in DIPG cell lines

# Description

The research aimed to analyze how ctDNA fragmentation patterns are affected by drug induced tumor cell death. DIPG cell line PBT29 was treated with four chemotherapeutic drugs: GB13, Quisinostat, Pevonedistat, and Temozolomide (TMZ). ctDNA was extracted from cell culture supernatants and sequenced using low-pass whole genome sequencing. 

# Features
  - Low pass whole genome sequencing using illumina.
  - **Slurm file:** contains the DNA sequencing pipelines, including the trimmed fastq reads, sequence alignment using BWA, filtering bam files using samtools, producing fragment size metrices using picard.
  - **R files:** contains the script for fragment pattern analysis using R, including different types of plots for fragment size distributions of PBT29 for the four drugs.
  
**The codes can be found in the master branch of this repository**
