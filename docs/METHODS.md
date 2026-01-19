# Methodology

## Data Collection
SARS-CoV-2 sequencing data was retrieved from NCBI SRA and ENA using Bioproject IDs.
A total of 157 Illumina sequencing samples were selected.

## Preprocessing & Quality Control
FASTQ files were assessed using FastQC to evaluate base quality, GC content,
adapter contamination, and sequence length.

## Reference Genome Preparation
The SARS-CoV-2 reference genome (NC_045512.2) was downloaded from NCBI GenBank
and indexed using BWA.

## Alignment & Variant Calling
Reads were aligned using BWA-MEM, processed using SAMtools, and variants
were called using BCFtools.

## Variant Annotation
Variants were annotated using SnpEff and analyzed for missense mutations.

