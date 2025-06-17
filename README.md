# SetDesign
SetDesign was designed to help researchers understand the impact of data processing choices on the stability and reproducibility of their set-based analyses in genome-wide association studies (GWAS). Specifically, SetDesign considers model misspecification caused by data processing steps. The package can calculate the power of the Sequece Kernel Association Test (SKAT) under misspecified models, and it can also calculate the bias in regression coefficient estimates under misspecified models. Users can specify parameters such as the number of subjects, SNP minor allele frequencies, significance levels, and effect sizes. The package supports both linear and binary outcomes and accounts for correlated SNPs by specifying a correlation coefficient for more precise results.

## Reporting Checklist for Set-Based Genetic Association Tests

To improve clarity and reproducibility of set-based test results, we recommend reporting at least the following data processing and analysis choices in the methods section:

* **Variant Inclusion Criteria:**
    * What Minor Allele Frequency (MAF) or Minor Allele Count (MAC) thresholds were applied for including variants in the analysis sets? (e.g., MAF < 0.01, MAC > 5)
    * Were different thresholds used for different types of analyses (e.g., single-variant vs. set-based)?

* **Variant Type Inclusion:**
    * What types of genetic variants were included in the sets? (e.g., Single Nucleotide Variants (SNVs) only, SNVs + multi-nucleotide variants (MNVs), SNVs + insertions/deletions (Indels), all variant types)

* **Handling of Multi-Allelic Variants:**
    * How were variants with more than two alleles handled? Choose one or describe method:
        * Distinct non-reference alleles treated as separate variants?
        * All non-reference alleles collapsed into a single "alternate" allele category?
        * Multi-allelic sites split into multiple bi-allelic records?
        * Excluded entirely?
        * Other method? (Please specify)

* **Set Definition Criteria:**
    * How were the boundaries of the variant sets defined? (e.g., specific gene coordinates from a reference genome build, gene coordinates +/- X kb upstream/downstream, predefined genomic regions)
    * If functional annotations were used to select variants for sets (e.g., "coding variants", "LoF variants", "regulatory variants"):
        * Which annotation source and version was used?
        * What specific annotation categories/terms were included in the definition for each set type?

*Reporting these specific choices clearly will help others understand exactly how the analysis was performed and aid in reproducing or comparing results across studies.*
