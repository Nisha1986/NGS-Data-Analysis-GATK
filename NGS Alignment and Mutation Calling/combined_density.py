import pysam
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# File paths for VCF and BAM files
normal_vcf_path = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/normal_variants.vcf.gz"
tumor_vcf_path = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/Tumor_variants_filtered.vcf.gz"
normal_bam_path = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/Normal_Sample.sorted.bam"
tumor_bam_path = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/Tumor_Sample.sorted.bam"

# Output file to save results
output_file = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/combined_mutation_density_output.txt"

# Step 1: Count number of mutations in the VCF (normal and tumor)
def count_mutations(vcf_path):
    return int(subprocess.check_output(f"bcftools view -H {vcf_path} | wc -l", shell=True).strip())

num_normal_mutations = count_mutations(normal_vcf_path)
num_tumor_mutations = count_mutations(tumor_vcf_path)

# Step 2: Open the BAM files and calculate total bases sequenced (only mapped bases)
def calculate_bases_and_reads(bam_path):
    bam_file = pysam.AlignmentFile(bam_path, "rb")
    total_bases = 0
    total_reads = 0
    for ref in bam_file.references:
        for read in bam_file.fetch(ref):
            if not read.is_unmapped:  # Only count mapped reads
                total_reads += 1
                total_bases += read.query_length  # Count the bases in the read
    bam_file.close()
    return total_bases, total_reads

normal_total_bases, normal_total_reads = calculate_bases_and_reads(normal_bam_path)
tumor_total_bases, tumor_total_reads = calculate_bases_and_reads(tumor_bam_path)

# Step 3: Calculate mutation density
def calculate_mutation_density(num_mutations, total_bases):
    return num_mutations / total_bases if total_bases != 0 else 0

normal_mutation_density = calculate_mutation_density(num_normal_mutations, normal_total_bases)
tumor_mutation_density = calculate_mutation_density(num_tumor_mutations, tumor_total_bases)

# Step 4: Normalize to mutations per million bases
def normalize_mutation_density(mutation_density):
    return mutation_density * 1_000_000

normal_mutations_per_million = normalize_mutation_density(normal_mutation_density)
tumor_mutations_per_million = normalize_mutation_density(tumor_mutation_density)

# Step 5: Output the results
with open(output_file, 'w') as out:
    out.write(f"Normal - Number of mutations: {num_normal_mutations}\n")
    out.write(f"Normal - Total bases sequenced: {normal_total_bases}\n")
    out.write(f"Normal - Total reads: {normal_total_reads}\n")
    out.write(f"Normal - Mutation density (mutations per base): {normal_mutation_density}\n")
    out.write(f"Normal - Normalized mutation density (mutations per million bases): {normal_mutations_per_million}\n")
    out.write(f"Tumor - Number of mutations: {num_tumor_mutations}\n")
    out.write(f"Tumor - Total bases sequenced: {tumor_total_bases}\n")
    out.write(f"Tumor - Total reads: {tumor_total_reads}\n")
    out.write(f"Tumor - Mutation density (mutations per base): {tumor_mutation_density}\n")
    out.write(f"Tumor - Normalized mutation density (mutations per million bases): {tumor_mutations_per_million}\n")

# Print results to console
print(f"Normal - Number of mutations: {num_normal_mutations}")
print(f"Normal - Total bases sequenced: {normal_total_bases}")
print(f"Normal - Total reads: {normal_total_reads}")
print(f"Normal - Mutation density (mutations per base): {normal_mutation_density}")
print(f"Normal - Normalized mutation density (mutations per million bases): {normal_mutations_per_million}")

print(f"Tumor - Number of mutations: {num_tumor_mutations}")
print(f"Tumor - Total bases sequenced: {tumor_total_bases}")
print(f"Tumor - Total reads: {tumor_total_reads}")
print(f"Tumor - Mutation density (mutations per base): {tumor_mutation_density}")
print(f"Tumor - Normalized mutation density (mutations per million bases): {tumor_mutations_per_million}")

# Step 6: Create graphs
mutation_types = ['Normal', 'Tumor']
mutation_densities = [normal_mutations_per_million, tumor_mutations_per_million]

# Create a bar plot for mutation density comparison
plt.figure(figsize=(8, 6))
sns.barplot(x=mutation_types, y=mutation_densities, palette='viridis')
plt.title("Mutation Density Comparison: Normal vs Tumor")
plt.ylabel("Mutations per Million Bases")
plt.xlabel("Sample Type")
plt.tight_layout()

# Save the plot
mutation_density_plot_path = "/home/hp/Desktop/NIsha/Pupil_Bio/pupil_bio_ngs_dataset/bwa/mutation_density_comparison.png"
plt.savefig(mutation_density_plot_path)

# Show the plot
plt.show()

