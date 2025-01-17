import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import numpy as np

# File paths for VCF files
normal_vcf_path = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/normal_filtered_variants.vcf.gz"
tumor_vcf_path = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/Tumor_variants_filtered.vcf.gz"

# Step 1: Extract allele frequencies for both normal and tumor samples
normal_allele_frequencies = []
tumor_allele_frequencies = []

# Extract allele frequencies from the normal VCF file
with pysam.VariantFile(normal_vcf_path, "r") as vcf:
    for record in vcf:
        sample = record.samples.get("Normal_Sample")
        if sample and "AF" in sample.keys():
            af_value = sample["AF"]
            if af_value and isinstance(af_value[0], float):
                normal_allele_frequencies.append(af_value[0])

# Extract allele frequencies from the tumor VCF file
with pysam.VariantFile(tumor_vcf_path, "r") as vcf:
    for record in vcf:
        sample = record.samples.get("Tumor_Sample")
        if sample and "AF" in sample.keys():
            af_value = sample["AF"]
            if af_value and isinstance(af_value[0], float):
                tumor_allele_frequencies.append(af_value[0])

# Step 2: Create DataFrame for plotting
import pandas as pd

data = {
    "Sample Type": ["Normal"] * len(normal_allele_frequencies) + ["Tumor"] * len(tumor_allele_frequencies),
    "Allele Frequency": normal_allele_frequencies + tumor_allele_frequencies
}

df = pd.DataFrame(data)

# Step 3: Generate box plot using Seaborn
plt.figure(figsize=(8, 6))
sns.boxplot(x="Sample Type", y="Allele Frequency", data=df)
plt.title("Comparison of Allele Frequency: Normal vs Tumor")
plt.show()

