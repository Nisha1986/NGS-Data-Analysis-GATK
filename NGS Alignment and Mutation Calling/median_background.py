import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import numpy as np
import pandas as pd

# File paths for VCF files
normal_vcf_path = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/normal_filtered_variants.vcf.gz"
tumor_vcf_path = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/Tumor_variants_filtered.vcf.gz"

# File paths for BAM files
normal_bam_path = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/Normal_Sample.sorted.bam"
tumor_bam_path = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/Tumor_Sample.sorted.bam"

# Output file to save results
output_file = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/output_file.txt"

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
try:
    with pysam.VariantFile(tumor_vcf_path, "r") as vcf:
        for record in vcf:
            # Check both record.info and record.samples for AF
            af_value = record.info.get("AF", None)
            if af_value is None and "Tumor_Sample" in record.samples:
                sample = record.samples["Tumor_Sample"]
                af_value = sample.get("AF", None)

            if af_value and isinstance(af_value[0], float):  # Ensure AF is a float
                tumor_allele_frequencies.append(af_value[0])  # Fix: Append to tumor_allele_frequencies
except Exception as e:
    print(f"Error reading VCF file: {e}")
    exit()

# Check if allele frequencies are empty
if not normal_allele_frequencies or not tumor_allele_frequencies:
    print("No valid allele frequencies found in the VCF files. Please check the input.")
    exit()

# Step 2: Calculate the median background mutation levels for both normal and tumor
normal_median_background = np.median(normal_allele_frequencies)
tumor_median_background = np.median(tumor_allele_frequencies)

# Function to calculate Reads Per Million (RPM)
def calculate_rpm(bam_path, median_background):
    bam_file = pysam.AlignmentFile(bam_path, "rb")
    total_bases = sum(read.query_alignment_length for read in bam_file.fetch() if not read.is_unmapped and read.query_alignment_length is not None)
    bam_file.close()
    rpm = (median_background * 1_000_000) / total_bases
    return rpm

# Calculate RPM for normal and tumor samples
normal_rpm = calculate_rpm(normal_bam_path, normal_median_background)
tumor_rpm = calculate_rpm(tumor_bam_path, tumor_median_background)

# Output the results
print(f"Normal - Median Background Mutation Level (AF): {normal_median_background}")
print(f"Normal - Reads Per Million Threshold (RPM): {normal_rpm}")
print(f"Tumor - Median Background Mutation Level (AF): {tumor_median_background}")
print(f"Tumor - Reads Per Million Threshold (RPM): {tumor_rpm}")

# Save results to a file
with open(output_file, "w") as out:
    out.write(f"Normal - Median Background Mutation Level (AF): {normal_median_background}\n")
    out.write(f"Normal - Reads Per Million Threshold (RPM): {normal_rpm}\n")
    out.write(f"Tumor - Median Background Mutation Level (AF): {tumor_median_background}\n")
    out.write(f"Tumor - Reads Per Million Threshold (RPM): {tumor_rpm}\n")

# Step 3: Create DataFrame for plotting
data = {
    "Sample Type": ["Normal"] * len(normal_allele_frequencies) + ["Tumor"] * len(tumor_allele_frequencies),
    "Allele Frequency": normal_allele_frequencies + tumor_allele_frequencies
}

df = pd.DataFrame(data)

# Step 4: Generate box plot using Seaborn
plt.figure(figsize=(8, 6))
sns.boxplot(x="Sample Type", y="Allele Frequency", data=df)
plt.title("Comparison of Allele Frequency: Normal vs Tumor")
plt.show()

# Step 5: Creating a plot for the comparison

# Data for normal and tumor samples
normal_af = normal_median_background
normal_rpm = normal_rpm
tumor_af = tumor_median_background
tumor_rpm = tumor_rpm

# Labels for the plot
labels = ['Normal', 'Tumor']
af_values = [normal_af, tumor_af]
rpm_values = [normal_rpm, tumor_rpm]

# Creating figure and axis for plotting
fig, ax1 = plt.subplots(figsize=(8, 6))

# Bar plot for Allele Frequency (AF)
bar_width = 0.35
index = np.arange(len(labels))
bar1 = ax1.bar(index, af_values, bar_width, label='Median AF', color='b', alpha=0.6)

# Adding second axis for RPM values
ax2 = ax1.twinx()
bar2 = ax2.bar(index + bar_width, rpm_values, bar_width, label='RPM Threshold', color='g', alpha=0.6)

# Set labels and title
ax1.set_xlabel('Sample Type')
ax1.set_ylabel('Median AF', color='b')
ax2.set_ylabel('Reads Per Million (RPM)', color='g')
ax1.set_title('Comparison of Normal and Tumor Mutation Levels')

# Adding legend
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

# X-ticks
ax1.set_xticks(index + bar_width / 2)
ax1.set_xticklabels(labels)

# Show plot
plt.tight_layout()
plt.show()

