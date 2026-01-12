# Modern BAM → Majority-rule Consensus FASTA

This repository contains a SLURM-based pipeline to process **modern resequencing BAM files**
and generate **majority-rule variant calls** as an intermediate step toward
**sample-specific consensus FASTA sequences**.

The workflow is designed to minimize reference bias by relying on observed
read data and explicit quality thresholds.

---

## Script

- `01.modern_bam_to_fasta_majority.sbatch`  
  SLURM batch script to generate majority-rule VCFs from BAM files using `bcftools`.

---

## Workflow Overview

<p align="center">
  <img src="bam_to_consensus_workflow.png" width="850">
</p>

<p align="center">
  <strong>Figure 1.</strong> Workflow used to convert mapped BAM files into
  sample-specific majority-rule consensus sequences.
  BAM files are summarized into VCFs using <code>bcftools mpileup</code> and
  <code>bcftools call -c</code>, followed by filtering, masking of low-confidence
  sites, and consensus sequence generation.
</p>

---

## Input Files

### Required inputs

- **BAM list file**

## Downstream SNP Merging and Window Definition

After per-sample SNP generation, SNPs are aggregated across individuals and
summarized in fixed genomic windows. This step is used to define a
reference-based windowing scheme and to quantify SNP density across the genome.

### Script

- `merge_snps_miss60_density.sbatch`

### Purpose

The purpose of this step is to:
- Merge per-sample SNP VCFs into a multi-sample VCF
- Filter sites by overall missingness
- Define fixed, reference-based genomic windows
- Quantify SNP density per window

This step is **not** intended to produce a final analysis-ready genotype matrix.

---

### Processing steps

1. Merge per-sample SNP VCFs using `bcftools merge`
2. Filter SNPs by site-level missingness (`F_MISSING ≤ 0.60`)
3. Retain SNPs only (indels excluded)
4. Define non-overlapping 10 kb windows based on the reference genome
5. Count SNPs per window to estimate genome-wide SNP density

---

### Treatment of multiallelic sites

Multiallelic SNPs are intentionally retained at this stage.

The merged VCF produced here is used exclusively for:
- Defining genomic windows
- Estimating SNP density
### Workflow diagram

![SNP merging and window-based density workflow](merge_snps_missingness_density_workflow.png)

SNPs are treated as **positional markers**, not as genotypes for phylogenetic or
population-genetic inference. Multiallelic filtering can be applied later,
during FASTA extraction, MSA generation, or tree inference, without affecting
window definitions.

---
## Window-based FASTA Extraction and Alignment Assembly

After genomic windows have been defined from the merged SNP VCF (either using
fixed-size windows or fixed-SNP windows), these windows are applied uniformly
across all samples to generate window-specific FASTA sequences and
multi-sample alignments.

This stage links SNP-based window definition to downstream sequence-based
analyses such as multiple sequence alignment (MSA) and phylogenetic inference.

---

### Input dependencies

This step assumes the following inputs have already been generated:

1. **Window definitions (BED format)**  
   Derived from the merged SNP VCF, using either:
   - Fixed physical windows (e.g. 10 kb), or
   - Fixed SNP-count windows (e.g. 50 SNPs per window)

2. **Per-sample FASTA files**  
   Sample-specific FASTA sequences generated in earlier steps.

All window definitions are reference-based and shared across samples.

---

### Script: window-based alignment construction

- `windows_to_alignments.sbatch`

This script constructs window-specific multi-sample FASTA alignments from
per-sample windowed FASTA files.

---

### Processing logic

1. **Per-sample window extraction**  
   For each sample, window coordinates are used to extract sequence fragments
   corresponding to each genomic window. Window identifiers are parsed from
   FASTA headers and preserved.

2. **Temporary per-sample storage**  
   Extracted window sequences are written to per-sample temporary directories,
   ensuring that each window and sample combination is tracked explicitly.

3. **Window-wise merging across samples**  
   For each genomic window, sequences from all samples are concatenated into a
   single FASTA file, producing one alignment-ready FASTA per window.

4. **Coordinate-consistent naming**  
   Output files are named using window indices and reference coordinates,
   ensuring traceability back to the reference genome.

---
## Window-based FASTA Extraction and Alignment Assembly

After genomic windows have been defined from the merged SNP VCF (either using
fixed-size windows or fixed-SNP windows), these windows are applied uniformly
across all samples to generate window-specific FASTA sequences and
multi-sample alignments.

This stage links SNP-based window definition to downstream sequence-based
analyses such as multiple sequence alignment (MSA) and phylogenetic inference.

---

### Input dependencies

This step assumes the following inputs have already been generated:

1. **Window definitions (BED format)**  
   Derived from the merged SNP VCF, using either:
   - Fixed physical windows (e.g. 10 kb), or
   - Fixed SNP-count windows (e.g. 50 SNPs per window)

2. **Per-sample FASTA files**  
   Sample-specific FASTA sequences generated in earlier steps.

All window definitions are reference-based and shared across samples.

---

### Script: window-based alignment construction

- `windows_to_alignments.sbatch`

This script constructs window-specific multi-sample FASTA alignments from
per-sample windowed FASTA files.

---

### Processing logic

1. **Per-sample window extraction**  
   For each sample, window coordinates are used to extract sequence fragments
   corresponding to each genomic window. Window identifiers are parsed from
   FASTA headers and preserved.

2. **Temporary per-sample storage**  
   Extracted window sequences are written to per-sample temporary directories,
   ensuring that each window and sample combination is tracked explicitly.

3. **Window-wise merging across samples**  
   For each genomic window, sequences from all samples are concatenated into a
   single FASTA file, producing one alignment-ready FASTA per window.

4. **Coordinate-consistent naming**  
   Output files are named using window indices and reference coordinates,
   ensuring traceability back to the reference genome.

---

### Outputs




### Outputs

