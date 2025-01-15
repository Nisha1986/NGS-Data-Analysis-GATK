# Pupil_Bio_PMP_Analysis
This task focuses on the statistical analysis and data handling of Phased Methylation Patterns (PMPs), a novel approach to analyzing CpG methylation data for tissue differentiation. The goal is to process complex sequencing data, evaluate variability in methylation patterns, and statistically assess their potential as reliable tissue-specific biomarkers. The task involves calculating coverage, identifying patterns of interest, and testing the hypothesis that PMPs provide higher specificity compared to individual CpG sites. Through data visualization, statistical tests, and machine learning approaches, the task aims to explore the relationship between sequencing depth, confidence in PMP specificity, and their utility in distinguishing tissue types.

CpG methylation is an epigenetic marker widely studied due to its role in gene regulation. However, the reliability of single CpG site methylation as a biomarker is low because of:

#Errors in bisulfite sequencing: Bisulfite sequencing can introduce false positives or negatives in methylation status determination.
#Sampling techniques: Inconsistent methods of sample collection and preparation may introduce variability.
#Biological variability: Methylation status can vary naturally across tissues or within cells of the same tissue.

A PMP is a structured representation of:

DNA strand orientation: Either forward (‘f’) or reverse (‘r’).
Relative positions of CpG sites: A triad of CpG sites (e.g., x:y:z, where x, y, and z are integers indicating positions).
Methylation statuses: A binary string indicating the methylation state at the three positions (e.g., ‘000’ for all unmethylated, ‘111’ for all methylated, or combinations like ‘101’).
PMPs act as composite biomarkers that can potentially differentiate tissue types with higher specificity compared to analyzing single CpG sites.

Hypothesis
#The assignment assumes that PMPs, by capturing combined methylation signatures across multiple CpG sites, provide:
#Higher specificity.
#Improved reliability as biomarkers for distinguishing between tissue types.

Binary Representation of Methylation Status
Each digit corresponds to the methylation state of a specific CpG site:

0: The site is unmethylated.
1: The site is methylated.
For example:

000: All three CpG sites are unmethylated.
111: All three CpG sites are methylated.
101: The first and third CpG sites are methylated, while the second is unmethylated.

Breakdown of the Methylation Statuses
000: Unmethylated at all three CpG sites.
100: Methylated at the first CpG site, unmethylated at the second and third.
111: Methylated at all three CpG sites.
101: Methylated at the first and third CpG sites, unmethylated at the second.
011: Unmethylated at the first CpG site, methylated at the second and third.
001: Unmethylated at the first and second CpG sites, methylated at the third.
010: Unmethylated at the first and third CpG sites, methylated at the second.
110: Methylated at the first and second CpG sites, unmethylated at the third.

Key Observations
Diversity of Patterns:

The eight patterns represent all possible combinations of methylation and unmethylation states for three CpG sites. This is consistent with 
=8 possible combinations.
Biological Relevance:
Each pattern might correspond to a unique phased methylation signature. These signatures can be specific to tissue types or experimental conditions.
Analysis Potential:
These patterns can be analyzed for their frequency across tissue types to identify tissue-specific PMPs.
Statistical analysis can help determine if certain patterns (e.g., 111 or 100) are enriched in specific tissue types, supporting their use as biomarkers.

The dataset now includes:

Two Tissues: cfDNA and Islet
The goal would be to see how methylation patterns differ between these tissues.
Two Strands: f (forward) and r (reverse)
Each strand might exhibit unique methylation patterns that should be analyzed separately.
Two Replicates: Rep1 and Rep2
These technical replicates help assess experimental consistency.

Analysis Steps

#Data Cleaning:
Ensure there are no missing values in the dataset.
Verify consistency in CpG_Coordinates and strand across replicates.

#Frequency Analysis:
Sum the counts across all patterns (000 to 111) for each set of CpG_Coordinates.
Identify the most frequent PMP for each set of coordinates.

#Tissue-Specific Patterns:
Compare the distribution of PMP frequencies across tissue types (cfDNA in this dataset).
Look for patterns that are uniquely or predominantly associated with cfDNA.

#Statistical Tests:
Hypothesis Testing
#Null Hypothesis: The PMP methylation pattern is not associated with tissue type.
#Alternative Hypothesis: The PMP methylation pattern is associated with tissue type.
Use a chi-square test or similar to determine if the distribution of PMPs is significantly different across replicates or tissue types.
If you have multiple tissue types, perform pairwise comparisons to identify tissue-specific PMPs.

Visualization:
#Heatmaps: Show the frequency distribution of PMPs (000 to 111) for each tissue type.
#Bar Graphs: Compare the relative abundance of PMPs across samples or replicates.

#Machine Learning Approach:
Identifies complex patterns across PMPs.
Suitable for larger datasets or when predicting tissue types is the primary goal.
#####LightGBM (Light Gradient Boosting Machine) algorithm########
LightGBM is particularly well-suited for handling imbalanced datasets, where class 0 is underrepresented compared to class 1. Its ability to handle large datasets efficiently 
and support advanced techniques like class weighting makes it a great choice.
Key Terms
Precision, Recall (Sensitivity),F1-score 
Support: The number of actual instances of each class in the dataset.

#Correlation with Tissue Type:
Use machine learning models (e.g., logistic regression, decision trees) to predict tissue types based on PMP frequencies.
Assess the specificity and sensitivity of PMPs in tissue differentiation.
Classification Report
Check precision, recall, and F1-score for each tissue class (0: Islet, 1: cfDNA).
Focus on recall for cfDNA (to minimize false negatives) and precision for Islet (to minimize false positives).
Feature Importance
Use the ranked PMPs to identify those most associated with tissue differentiation.
PMPs with higher importance scores are potential biomarkers.
AUC-ROC
A high AUC-ROC score indicates good differentiation between tissues.
