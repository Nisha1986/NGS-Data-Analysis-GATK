#!/bin/bash
################################################################################
# Task 2: NGS Data Analysis
# Objective: Process and analyze raw sequencing data from normal and tumor 
# samples (paired-end, 2×100 bp, Illumina).
# Steps: QC, trimming, alignment, variant calling (Mutect2), filtering, 
# and extracting variants for normal/tumor.
# Reference: GATK4 best practices 
# (https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)
################################################################################

# Data files:
# Tumor:   PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz 
#          PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz
# Normal:  PA221MH-lib09-P19-Norm_S1_L0es

01_R1_001.fastq.gz 
#          PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz

################################################################################
# Step 0: Decompress FASTQ files if necessary
################################################################################
# NOTE: It's often fine to keep FASTQ files compressed for tools like fastqc or
# fastp. If a tool requires uncompressed files, you can run gunzip.
# However, if you want to keep them compressed, remove or comment out this line.

gunzip *.gz

echo "FASTQ files decompressed."

################################################################################
# Step 1: Prep files (Download and index reference) - only once
################################################################################

echo "Preparing reference files..."

# 1a. Download reference genome (hg38) from UCSC:
mkdir -p ~/Desktop/demo/supporting_files/hg38
wget -P ~/Desktop/demo/supporting_files/hg38/ \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gunzip ~/Desktop/demo/supporting_files/hg38/hg38.fa.gz

# 1b. Index reference with samtools
samtools faidx ~/Desktop/demo/supporting_files/hg38/hg38.fa

# 1c. Create sequence dictionary with GATK
gatk CreateSequenceDictionary \
  -R ~/Desktop/demo/supporting_files/hg38/hg38.fa \
  -O ~/Desktop/demo/supporting_files/hg38/hg38.dict

# 1d. Download known sites (e.g., dbsnp for BQSR or annotation) 
# from GATK resource bundle:
wget -P ~/Desktop/demo/supporting_files/hg38/ \
  https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

wget -P ~/Desktop/demo/supporting_files/hg38/ \
  https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

echo "Reference download and indexing complete."

################################################################################
# Step 2: Quality Control with FastQC
################################################################################

echo "STEP 2: Quality Control (FastQC)"

# Install FastQC if not installed (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install -y fastqc

# Run FastQC on all .fastq files in the directory
fastqc *.fastq

################################################################################
# Step 3: Trimming (fastp)
################################################################################

echo "STEP 3: Trimming with fastp"

# Install fastp if not installed (Ubuntu/Debian)
sudo apt-get install -y fastp

# Example trimming command for tumor R1/R2
fastp \
  -i PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq \
  -I PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq \
  -o Trim_PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq \
  -O Trim_PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq \
  --adapter_fasta adapter1.fasta \
  --cut_mean_quality 30

################################################################################
# Step 4: Align to reference using BWA-MEM
################################################################################

echo "STEP 4: Align to reference using BWA-MEM"

cd ~/Desktop/demo/supporting_files/hg38/
bwa index hg38.fa
cd -

# BWA alignment for Normal Sample
bwa mem -t 8 ~/Desktop/demo/supporting_files/hg38/hg38.fa \
  Trim_PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq \
  Trim_PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq \
  > Normal_Sample.sam

# BWA alignment for Tumor Sample
bwa mem -t 8 ~/Desktop/demo/supporting_files/hg38/hg38.fa \
  Trim_PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq \
  Trim_PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq \
  > Tumor_Sample.sam

################################################################################
# Step 5: Convert SAM to BAM and sort
################################################################################

echo "STEP 5: Convert and sort SAM files"

samtools view -@ 8 -bS Tumor_Sample.sam \
  | samtools sort -@ 8 -o Tumor_Sample.sorted.bam

samtools view -@ 8 -bS Normal_Sample.sam \
  | samtools sort -@ 8 -o Normal_Sample.sorted.bam

################################################################################
# Step 6: Index the BAM files
################################################################################

echo "STEP 6: Index the BAM files"

samtools index Tumor_Sample.sorted.bam
samtools index Normal_Sample.sorted.bam

################################################################################
# Step 7: Verify Alignments with samtools flagstat
################################################################################

echo "STEP 7: Verify alignment metrics with samtools flagstat"

samtools flagstat Tumor_Sample.sorted.bam
samtools flagstat Normal_Sample.sorted.bam


###############################################################################
# STEP 8: Add or Replace Read Groups (Picard) for Normal and Tumor
###############################################################################

echo "=== STEP 8: Adding/Replacing Read Groups for Normal and Tumor ==="

# For Normal Sample:
java -jar ./picard.jar AddOrReplaceReadGroups \
    I=/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/Normal_Sample.sorted.bam \
    O=/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/Normal_Sample_RG.bam \
    RGID=group1 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=Normal_Sample

echo "Read groups added to Normal_Sample.sorted.bam, output: Normal_Sample_with_RG.bam"

# For Tumor Sample:
java -jar ./picard.jar AddOrReplaceReadGroups \
    I=/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/Tumor_Sample.sorted.bam \
    O=/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/Tumor_Sample_RG.bam \
    RGID=group2 \
    RGLB=lib2 \
    RGPL=illumina \
    RGPU=unit2 \
    RGSM=Tumor_Sample

echo "Read groups added to Tumor_Sample.sorted.bam, output: Tumor_Sample_with_RG.bam"

################################################################################
# Step 9: Somatic Variant Calling with GATK Mutect2
################################################################################

echo "STEP 9: Somatic Variant Calling - GATK Mutect2"

gatk Mutect2 \
  -R ~/Desktop/demo/supporting_files/hg38/hg38.fa \
  -I Tumor_Sample.sorted.bam \
  -I Normal_Sample.sorted.bam \
  --tumor-sample TUMOR_SAMPLE \
  --normal-sample NORMAL_SAMPLE \
  -O raw_variants.vcf

# Output: raw_variants.vcf, .vcf.idx, .vcf.stats

################################################################################
# Step 10: Filter the raw variants
################################################################################

echo "STEP 10: Filter the raw variants with FilterMutectCalls"

gatk FilterMutectCalls \
  -R ~/Desktop/demo/supporting_files/hg38/hg38.fa \
  -V raw_variants.vcf \
  -O filtered_variants.vcf

# Compress the resulting VCF
bgzip filtered_variants.vcf
tabix -p vcf filtered_variants.vcf.gz

################################################################################
# Step 11: Extract normal/tumor subsets (optional)
################################################################################

echo "Extract normal/tumor variants from multi-sample VCF"

bcftools view \
  --samples NORMAL_SAMPLE \
  --output Normal_variant.vcf.gz \
  --output-type z \
  filtered_variants.vcf.gz

bcftools view \
  --samples TUMOR_SAMPLE \
  --output Tumor_variant.vcf.gz \
  --output-type z \
  filtered_variants.vcf.gz

# Index these new VCFs
bcftools index Normal_variant.vcf.gz
bcftools index Tumor_variant.vcf.gz

################################################################################
# Step 12: Filter the new normal/tumor VCFs further if needed
################################################################################

echo "STEP 12: Filter the new normal/tumor VCFs (if necessary)"

gatk FilterMutectCalls \
  -R ~/Desktop/demo/supporting_files/hg38/hg38.fa \
  -V Normal_variant.vcf.gz \
  -O Normal_filtered_variants.vcf.gz

gatk FilterMutectCalls \
  -R ~/Desktop/demo/supporting_files/hg38/hg38.fa \
  -V Tumor_variant.vcf.gz \
  -O Tumor_filtered_variants.vcf.gz
################################################################################
# STEP 13 : Generate bcftools stats for Normal and Tumor
################################################################################

echo "=== Generating stats for Normal ==="
bcftools stats Normal_filtered_variants.vcf.gz > normal_mutation_stats.txt

echo "=== Generating stats for Tumor ==="
bcftools stats Tumor_filtered_variants.vcf.gz > tumor_mutation_stats.txt

echo "bcftools stats completed. Summaries in normal_mutation_stats.txt and tumor_mutation_stats.txt."

echo "Pipeline complete."
