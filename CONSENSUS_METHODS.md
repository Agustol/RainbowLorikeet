# Consensus Methods Note: Legacy vs Modern Approaches

This document explains why this project uses a modern, VCF-driven consensus
workflow (`bcftools consensus`) rather than legacy FASTQ consensus generation
(`vcfutils.pl vcf2fq`).

---

## Comparison Figure

<p align="center">
  <img src="legacy_vs_modern_consensus_comparison.png" width="900">
</p>

---

## Legacy approach: `mpileup → vcfutils.pl vcf2fq`

**Summary**
- Generates a consensus FASTQ directly from BAM + reference
- Typically collapses heterozygous sites implicitly
- Often masks low-depth sites aggressively

**Limitations (relevant to phylogenomic windows)**
- Collapses heterozygous sites without explicit policy control
- Can introduce reference bias
- Poor / brittle handling of indels in many workflows
- Depends on deprecated or brittle scripts
- Can inflate similarity among samples, affecting window-based analyses

---

## Modern approach: `bcftools consensus`

**Summary**
- Separates variant calling from consensus generation
- Uses an explicit VCF as input plus the reference FASTA
- Allows explicit control of missing data and heterozygosity policies

**Advantages**
- Reproducible and transparent (VCF-driven)
- Better supported / maintained tooling
- Explicit policies for heterozygotes and missing data
- Robust handling of indels
- More appropriate for window-based FASTA extraction and downstream alignments

---

## Implications for this pipeline

Because this project relies on:
- SNP-defined genomic windows
- Window-based FASTA extraction across many samples
- Alignment-ready per-window FASTA outputs

the modern VCF-based consensus strategy improves interpretability and reduces
bias compared to legacy FASTQ-based approaches.

For these reasons, legacy `vcfutils.pl vcf2fq` workflows are not used in this
pipeline.

# Consensus Methods Note: Three way comparisson 

<p align="center">
  <img src="fasta_generation_lc _methods.png" width="900">
</p>

## INPUT DATA
----------
BAM (Low Coverage)
Low-coverage BAM files contain substantial uncertainty at many genomic positions.
How these uncertain positions are treated determines whether reference bias is
introduced during consensus generation.

----------------------------------------------------------------
1) VARIANT-ONLY CALLING + PATCHING (bcftools, variants only)
----------------------------------------------------------------
• Variant calling is restricted to polymorphic sites.
• Invariant sites are not explicitly represented in the VCF.
• During consensus generation, positions without variant records are filled
  with the reference base (“patching”).
• Low-coverage or missing positions are typically masked only after consensus
  generation.

Key assumption:
  Absence of a variant call implies identity with the reference genome.

Consequence:
  Reference bases may be introduced at sites with no supporting read evidence
  unless those sites are explicitly masked downstream.

----------------------------------------------------------------
2) ALL-SITES CALLING (-A) + MASKING (NON-PATCHING, bcftools)
----------------------------------------------------------------
• Variant calling records both variant and invariant sites.
• Invariant sites (0/0) are retained only when supported by sequencing reads.
• Low-coverage sites and heterozygous genotypes are explicitly converted to
  missing (./.) prior to consensus generation.
• During consensus generation, missing genotypes are emitted as N rather than
  reference bases.

Key assumption:
  Absence of evidence is treated as uncertainty, not as reference identity.

Consequence:
  Invariant sites are preserved only when justified by data, eliminating
  implicit reference bias.

----------------------------------------------------------------
3) ANGSD doFasta 2 (GENOTYPE LIKELIHOOD APPROACH)
----------------------------------------------------------------
• Uses genotype likelihoods instead of hard genotype calls.
• High-probability bases are emitted as A/C/G/T.
• Low-confidence positions are masked as N.

Key assumption:
  Bases are emitted only when posterior support exceeds a defined threshold.

Consequence:
  Conceptually equivalent to the all-sites non-patching approach when similar
  masking thresholds are applied.

----------------------------------------------------------------
EQUIVALENCE CONDITION
----------------------------------------------------------------
When identical masking rules are applied (e.g. DP ≥ 3×, heterozygous sites masked):

• Variant sites resolve to ALT in all methods.
• Supported invariant sites resolve to REF in all methods.
• Low-coverage or missing sites resolve to N in all methods.

Under these conditions, patched + masked and all-sites + masked pipelines
produce identical FASTA outputs despite differing upstream assumptions.

----------------------------------------------------------------
INTERPRETATION TABLE
----------------------------------------------------------------
Site category             Patched + Masked        -A + Masked
---------------------------------------------------------------
Variant site              ALT                     ALT
Supported invariant (0/0) REF                     REF
No data / low DP          N                       N

----------------------------------------------------------------
KEY TAKEAWAY
----------------------------------------------------------------
Patching and non-patching consensus methods differ conceptually in how missing
data are interpreted. However, once low-information sites are masked using the
same criteria, both approaches converge on the same set of informative genomic
positions and therefore produce equivalent FASTA sequences.
