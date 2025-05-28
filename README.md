............NGS-Data-Analysis-GATK.......................

🔬 Project Overview

This repository contains the complete pipeline and scripts used for Next-Generation Sequencing.
g (NGS) data analysis using GATK (Genome Analysis Toolkit) to perform alignment and mutation calling between tumor and normal samples.
The goal of this project is to identify somatic variants by comparing sequencing data from tumor tissue with matched normal tissue.

NGS-Data-Analysis-GATK/
│
├── data/                       # Input FASTQ or BAM files

├── reference/                 # Reference genome (e.g., hg38 or hg19) and index files

├── scripts/                   # Bash and Python scripts for analysis

├── results/                   # Output files: BAMs, VCFs, logs

├── README.md                  # This file

└── pipeline.sh                # Main pipeline script

Tools & Technologies

GATK (v4.x)

BWA-MEM – For sequence alignment

Samtools – For BAM file processing

Picard Tools – For marking duplicates

GATK Mutect2 – For somatic variant calling

FastQC – Quality check

MultiQC – Aggregated report

bcftools – For VCF filtering and stats

🚀 Workflow Steps

Quality Control
Raw FASTQ files are checked using FastQC and summarized with MultiQC.

Reference Genome Preparation
Indexed reference genome files are prepared using BWA and Samtools.

Read Alignment
FASTQ reads are aligned to the reference genome using BWA-MEM.

Post-processing

Sorting and indexing of BAM files

Marking duplicates using Picard

Base Quality Recalibration (BQSR)
Performed using GATK BaseRecalibrator and ApplyBQSR.

Somatic Variant Calling
Mutect2 is used to identify somatic mutations by comparing tumor vs. normal BAM files.

Variant Filtering & Annotation

GATK FilterMutectCalls

Optional: Variant annotation using VEP or ANNOVAR

📊 Output Files
BAM files (sorted, duplicates marked, indexed)

Recalibration tables

Raw and filtered VCFs

Summary reports from FastQC and MultiQC

📌 Prerequisites

Linux/Unix-based system

Conda or Docker for environment management (optional but recommended)

Reference genome (e.g., hg38.fa with BWA index, .dict, and .fai files)

GATK, BWA, Samtools, Picard, FastQC installed

🧪 Sample Commands

# Run BWA alignment
bwa mem -t 8 reference/hg38.fa data/tumor_R1.fastq.gz data/tumor_R2.fastq.gz > results/tumor.sam

# Convert, sort, and index BAM
samtools view -Sb results/tumor.sam | samtools sort -o results/tumor_sorted.bam
samtools index results/tumor_sorted.bam

# Mark duplicates
picard MarkDuplicates I=results/tumor_sorted.bam O=results/tumor_dedup.bam M=results/metrics.txt

# Mutect2 variant calling
gatk Mutect2 -R reference/hg38.fa -I results/tumor_dedup.bam -tumor TumorSample \
             -I results/normal_dedup.bam -normal NormalSample \
             -O results/unfiltered.vcf

# Filtering
gatk FilterMutectCalls -V results/unfiltered.vcf -R reference/hg38.fa -O results/filtered.vcf

📌 Citation

If you use this pipeline, please cite:

GATK: https://gatk.broadinstitute.org

BWA: Li H. (2009). Bioinformatics

Samtools: Li H. et al. (2009). Bioinformatics
