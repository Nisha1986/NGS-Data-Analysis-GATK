NGS Alignment and Mutation Calling
This repository contains the NGS Data Analysis Project, which focuses on Task 2 of a bioinformatics challenge. The primary objective is to process and analyze raw sequencing data (FASTQ files) from both normal and tumor samples using various bioinformatics tools. The workflows documented here include quality control (via FastQC), alignment (e.g., using BWA), somatic variant calling (Mutect2, VarScan2, etc.), and calculation of a background mutation level from the normal sample to ensure high-confidence identification of truly somatic mutations.

Below is a concise objectives list for your NGS data analysis task, rewritten in bullet points for clarity:

Objective 1: Perform quality checks on raw FASTQ files (tumor and normal) to assess read quality, adapter contamination, and duplication.

Objective 2: Align paired-end reads to the human reference genome (e.g., hg38) using alignment tools such as Bowtie2 or BWA.

                  Objective 3: Identify somatic mutations present in the tumor sample but absent (or at much lower frequency) in the normal sample.

                  Benchmark Software: Employing established tools (e.g., Mutect2, Strelka2, VarScan2).Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2

Custom Scripts: Leveraging command-line utilities (Samtools, bcftools) and Python/R libraries for mutation detection.

Objective 4: Estimate the background mutation level using the normal sample, accounting for sequencing errors or biases, and determine the minimum reads-per-million threshold for confident variant calls.

Data Description

- **Samples**:
  - **Normal Tissue**: 
    - `PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz`
    - `PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz`
  - **Cancer Tissue**:
    - `PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz`
    - `PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz`
- **Platform**: Illumina, Paired-End (2×100 bp)
- **Reference Genome**: Typically `hg38` or `hg19` (indexed)

## Sub-Tasks

1. **Quality Control**  
   Perform checks using FastQC and summarize quality metrics (e.g., sequence counts, per-base quality, read duplication levels).

2. **Alignment and Mutation Calling**  
   a. Align the samples to the human genome.  
   b. Identify somatic mutations in the tumor sample vs. normal tissue.  
      - **Benchmark Software**: Use tools like Mutect2, Strelka2, or VarScan2.  
      - **Custom Code**: Scripts leveraging Samtools, bcftools, Python/R.  
   c. Calculate the **median background mutation level**using the normal tissue.

---

## Quality Control

1. **Tool**: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. **Example Commands**:
   ```bash
   # If necessary, decompress the FASTQ first:
   gunzip -c Tumor_R1.fastq.gz > Tumor_R1.fastq

   # Then run FastQC:
   fastqc Tumor_R1.fastq --outdir fastqc_results/

**Align the samples to the human genome** using tools like **Bowtie2** or **BWA**.
   - Examples:
     ```bash
     # Using BWA
     bwa mem -t 8 hg38.fa sample_R1.fastq.gz sample_R2.fastq.gz > sample.sam
     samtools view -bS sample.sam | samtools sort -o sample.sorted.bam
     samtools index sample.sorted.bam

     # OR Using Bowtie2
     bowtie2 -x hg38 -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -S sample.sam
     ...
     ```

2. **Identify somatic mutations** present in the cancer sample but absent in the normal tissue.
3. Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
   **i. Benchmark Software**: Use established tools such as **Mutect2**, **Strelka2**, or **VarScan2** for somatic mutation identification and background mutation estimation.
   ```bash
   gatk Mutect2 \
       -R hg38.fa \
       -I tumor.bam -tumor TUMOR_SAMPLE \
       -I normal.bam -normal NORMAL_SAMPLE \
       -O raw_variants.vcf
References
FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
BWA: http://bio-bwa.sourceforge.net/
GATK: https://gatk.broadinstitute.org/
bcftools: http://samtools.github.io/bcftools/   
