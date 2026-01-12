## Downstream SNP Merging and Window Definition

After per-sample SNP generation, SNPs are aggregated across individuals and
summarized in genomic windows. This step defines a reference-based windowing
scheme and quantifies SNP density across the genome.

### Script

- `merge_snps_miss60_density.sbatch`

### Purpose

This step is used to:

- Merge per-sample SNP VCFs into a multi-sample VCF
- Filter SNP sites by overall missingness
- Define fixed, reference-based genomic windows
- Quantify SNP density per window

This step is **not** intended to produce a final genotype matrix for inference.

---

### Processing steps

1. Merge per-sample SNP VCFs using `bcftools merge`
2. Filter SNPs by site-level missingness (`F_MISSING ≤ 0.60`)
3. Retain SNPs only (indels excluded)
4. Define non-overlapping genomic windows based on the reference genome
5. Count SNPs per window to estimate genome-wide SNP density

---

### Treatment of multiallelic sites

Multiallelic SNPs are intentionally retained at this stage.

The merged VCF is used exclusively to:
- Define genomic windows
- Estimate SNP density

SNPs are treated as **positional markers**, not as genotypes for phylogenetic or
population-genetic inference. Multiallelic filtering can be applied later
without affecting window definitions.

### Workflow diagram

![SNP merging and window-based density workflow](merge_snps_missingness_density_workflow.png)

---

## SNP-based Window Definition and FASTA Extraction

This step converts SNP-defined windows into reference-coordinate BED files and
extracts the corresponding sequences from per-sample FASTA files.

### Script

- `03.GetAlign_make_bed_windows_and_fastas.sbatch`

---

### Inputs

This step requires:

1. A **merged multi-sample SNP VCF**
2. A **reference genome FASTA**
3. **Per-sample FASTA files**

The reference FASTA is used to validate chromosome boundaries and ensure that
all window coordinates are consistent with the reference genome.

---

### Processing overview

1. Define genomic windows based on a fixed number of consecutive SNPs
   (e.g. 50 SNPs per window), resetting at chromosome boundaries
2. Validate window coordinates against reference chromosome lengths
3. Generate a BED file describing all SNP-based windows
4. Apply the BED windows uniformly to each per-sample FASTA file

This produces per-sample FASTA files in which each sequence corresponds to a
single SNP-defined genomic window expressed in reference coordinates.

---

### Conceptual note

Although window boundaries are derived from SNP positions in the merged VCF,
all windows are expressed in reference genome coordinates. This ensures that:

- Window definitions are independent of sample-specific FASTA content
- FASTA extraction is purely coordinate-based
- No circularity is introduced between SNP selection and sequence analysis

---

## Window-based Alignment Assembly

After window-specific FASTA files have been generated for each sample, sequences
are combined across samples to produce one FASTA file per genomic window.

### Script

- `04.windows_to_alignments.sbatch`

---

### Processing logic

1. Extract window-specific sequences for each sample
2. Store sequences temporarily on a per-sample basis
3. Merge sequences across samples for each genomic window
4. Generate one FASTA file per window containing all samples

---

### Outputs

