# SARS-CoV-2 Variant Analysis Pipeline

This repository contains a reproducible bioinformatics pipeline for
SARS-CoV-2 variant detection and missense mutation identification
using NGS data.

The pipeline was developed as part of a **Bioinformatics & Molecular Genetics Internship**
and focuses on mutation discovery, interpretation, and lineage annotation.

---

##  Overview of the Workflow

1. Download SARS-CoV-2 reference genome (NC_045512.2)
2. Download sequencing data from NCBI SRA
3. Convert SRA to paired-end FASTQ files
4. Align reads to reference genome using **BWA-MEM**
5. Process alignments using **SAMtools**
6. Call variants using **bcftools**
7. Filter variants by quality (QUAL ≥ 30)
8. Generate a consensus genome
9. Annotate mutations using **Nextclade**
10. Extract amino-acid substitutions (missense mutations)

---

##  Missense Mutation Identification

Missense mutations are identified using **Nextclade** annotation.
These are reported in the `aaSubstitutions` column of the Nextclade output TSV.

Example mutations detected include:

- **Spike protein**: D614G, N501Y, P681R
- **Nucleocapsid**: R203K, G204R
- **ORF1a / ORF1b** substitutions

This approach avoids reliance on custom SnpEff databases and uses
curated SARS-CoV-2 annotation standards.

---

## Repository Structure

sars-cov-2-variant-analysis/
  scripts/sars_cov2_variant_pipeline.sh
  docs/METHODS.md
  .gitignore
  README.md

---

## How to Run the Pipeline

```bash
bash scripts/sars_cov2_variant_pipeline.sh
```

Required tools:
sra-tools

bwa

samtools

bcftools

vcftools

nextclade
