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

SNPs are treated as **positional markers**, not as genotypes for phylogenetic or
population-genetic inference. Multiallelic filtering can be applied later,
during FASTA extraction, MSA generation, or tree inference, without affecting
window definitions.

---

### Outputs

