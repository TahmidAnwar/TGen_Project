#!/bin/bash

#SBATCH --job-name=test_bwa   # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mstampar@tgen.org     # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=25gb                     # Job memory request
#SBATCH --time=07:45:00               # Time limit hrs:min:sec
#SBATCH --output=Test_bwa.log      # Standard output and error log
#SBATCH -p compute 

uname -a

# Load the needed modules

module load BWA/0.7.17-GCCcore-5.4.0
module load SAMtools 
module load picard


# Define variables
ref_genome=/home/mstampar/tools/ref_genome/hg38/genome.fa
read1=/home/mstampar/BWA_Test_For_Tahmid/001_TP1_R1.fastq.gz
read2=/home/mstampar/BWA_Test_For_Tahmid/001_TP1_R2.fastq.gz
num_threads=10

# Summary overview paper of fragmentomics
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9903826/pdf/main.pdf

# Very basic video of how illumina sequencing works and the structure of the read coming from the sequencer:
# https://www.youtube.com/watch?v=fCd6B5HRaZ8
# basic structure of read1: adapter + index + fragment from sample
# basic structure of read2: fragment from sample + index + adapter

# 1. Trim the reads (currently I don't see the fastp module on gemini, but these reads have already been trimmed i believe)

# When read length exceeds DNA insert size, a run can sequence beyond the DNA insert and read bases from
# the sequencing adapter. To prevent these bases from appearing in FASTQ files, the adapter sequence is
# trimmed from the 3′ends of reads. Trimming the adapter sequence improves alignment accuracy and
# performance in Illumina FASTQ generation pipelines.

#fastp --in1 tp001_R1.fastq.gz --in2 tp001_R2.fastq.gz --out1 tp001_R1_clipped.fq.gz --out2 tp001_R2_clipped.fq.gz -j tp001_clipped.json -h tp001_clipped.html -Q --thread 4 -z 9

# 2. Align fastq files to produce bam file (manual with most commonly used flags in bwa commands: https://bio-bwa.sourceforge.net/bwa.shtml)

## A little more background on what is in sam format file header https://medium.com/@shilparaopradeep/samtools-guide-learning-how-to-filter-and-manipulate-with-sam-bam-files-2c28b25d29e8 

# Next 3 steps can be done in one liner:

# bwa mem -t 10 /home/mstampar/tools/ref_genome/hg38/genome.fa   /home/mstampar/BWA_Test_For_Tahmid/001_TP1_R1.fastq.gz /home/mstampar/BWA_Test_For_Tahmid/001_TP1_R2.fastq.gz
  | # samtools view -bSh - | samtools sort -O bam -T tempTP01 -o /scratch/mstampar/BWA_Test_For_Tahmid/TP01_test.bam

## Use bwa tool to align 70bp-1Mbp paired end reads (R1 and R2 fastq.gz) agains reference fasta file (hg30 genome.fa).
## You may want to use the option of multiple threads to make this step quicker.
## The output file will be .sam format file

bwa mem -t 10  $ref_genome  $read1  $read2 > /scratch/mstampar/BWA_Test_For_Tahmid/TP01.sam

## Use samtools (manual with the samtools commands: https://www.htslib.org/doc/samtools.html and more specific samtools view options :https://www.htslib.org/doc/samtools-view.html )
## to convert the sam file into bam file format
## Input should be the above created sam and output should be the bam file. You will want to use options that will probusce BAM format with header.
samtools view -bh /scratch/mstampar/BWA_Test_For_Tahmid/TP01.sam > /scratch/mstampar/BWA_Test_For_Tahmid/TP01.bam

## Use samtools (manual with the samtools commands: https://www.htslib.org/doc/samtools.html and more specific samtools sort options :https://www.htslib.org/doc/samtools-sort.html)
## to sort the above created bam file into sorted .bam file
## You will want to use the options to specify the nanme os the output and the format of the output
samtools sort -O bam -T tempTP01 -o /scratch/mstampar/BWA_Test_For_Tahmid/TP01.bam /scratch/mstampar/BWA_Test_For_Tahmid/TP01.bam

## Use samtools to index the bam file so you can view it in IGV 
samtools index /scratch/mstampar/BWA_Test_For_Tahmid/TP01.bam

# 3. Filtering the bam files.
## Use samtools to get FLAG statistics (calculator of SAM flags: https://broadinstitute.github.io/picard/explain-flags.html)

samtools flagstat  /scratch/mstampar/BWA_Test_For_Tahmid/TP01.bam 

## Use samtools or awk to extract the insert size data using these criteria: (a) reads have the “proper pair” flag,
## (b) mapping quality >= 30, and (c) CIGAR string doesn’t include “I” or “D”.
## https://www.biostars.org/p/253774/

## for (a) and (b)
samtools view -b -f 2 -q 30 /scratch/mstampar/BWA_Test_For_Tahmid/TP01.bam > /scratch/mstampar/BWA_Test_For_Tahmid/TP01_pp_30.bam 
samtools flagstat  /scratch/mstampar/BWA_Test_For_Tahmid/TP01_pp_30.bam 
samtools index /scratch/mstampar/BWA_Test_For_Tahmid/TP01_pp_30.bam 

## for (c)
samtools view -h /scratch/mstampar/BWA_Test_For_Tahmid/TP01_pp_30.bam | awk '$1 ~ "^@" || $6 !~ "I|D"' | samtools view -b - > /scratch/mstampar/BWA_Test_For_Tahmid/TP01_filtered.bam 
samtools flagstat  /scratch/mstampar/BWA_Test_For_Tahmid/TP01_filtered.bam
samtools index /scratch/mstampar/BWA_Test_For_Tahmid/TP01_filtered.bam 

# 4. Use picard tool as described in fragmentomics R library tutorial to produce the insert size metric .txt file from bam file produced in previous step.
picard CollectInsertSizeMetrics I=/scratch/mstampar/BWA_Test_For_Tahmid/TP01_filtered.bam O=/scratch/mstampar/BWA_Test_For_Tahmid/insert_size_metrics_TP01.txt H=/scratch/mstampar/BWA_Test_For_Tahmid/TP01_insert_size_histogram.pdf M=0.5
