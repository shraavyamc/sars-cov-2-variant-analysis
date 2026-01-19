# SARS-CoV-2 Variant Analysis Pipeline

This repository contains a bioinformatics workflow developed during my
**Bioinformatics & Molecular Genetics Internship at Insilicomics (Feb–May 2025)**.

The pipeline performs NGS-based variant calling and missense mutation analysis
on SARS-CoV-2 sequencing data.

## Workflow Overview
SRA → FASTQ → QC → Alignment → BAM → Variant Calling → VCF → Annotation

## Tools Used
- BWA
- SAMtools
- BCFtools
- VCFtools
- SnpEff
- FastQC

## Key Mutations Studied
- D614G
- N501Y
- P681R

## Disclaimer
This repository is for academic and computational research purposes only.
