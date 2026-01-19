#!/bin/bash
set -e  # Exit immediately if a command fails

############################
# VARIABLES
############################

SRA_ID="SRR29011221"

FASTA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/NC_045512.2.fasta"
FASTA="NC_045512.2.fasta"

SRA_FILE="${SRA_ID}.sra"
FASTQ1="${SRA_ID}_1.fastq"
FASTQ2="${SRA_ID}_2.fastq"

SAM_FILE="${SRA_ID}_aligned.sam"
BAM_FILE="${SRA_ID}_aligned.bam"
SORTED_BAM="${SRA_ID}_aligned_sorted.bam"

PILEUP_FILE="${SRA_ID}.pileup"
VCF_FILE="${SRA_ID}.vcf.gz"
FILTERED_VCF="${SRA_ID}_filtered.vcf"
ANNOTATED_VCF="${SRA_ID}_annotated.vcf"

############################
# DOWNLOAD REFERENCE GENOME
############################

REF_FASTA="NC_045512.2.fasta"
REF_URL="https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_045512.2&db=nuccore&report=fasta"

echo "Downloading SARS-CoV-2 reference genome..."
wget -O "$REF_FASTA" "$REF_URL" || { echo "Failed to download reference genome"; exit 1; }

############################
# FETCH SRA DATA
############################
echo "Fetching SRA data..."
prefetch "$SRA_ID" || { echo "Failed to fetch SRA data"; exit 1; }

# -----------------------------
# Convert SRA to FASTQ
# -----------------------------
echo "Converting SRA to FASTQ..."
fasterq-dump "$SRA_ID" --split-files --threads 4 || { echo "Failed to convert SRA to FASTQ"; exit 1; }

FASTQ1="${SRA_ID}_1.fastq"
FASTQ2="${SRA_ID}_2.fastq"

if [[ ! -f "$FASTQ1" || ! -f "$FASTQ2" ]]; then
    echo "FASTQ files not found after fasterq-dump"
    exit 1
fi

############################
# INDEX REFERENCE GENOME
############################

echo "Indexing reference genome..."
bwa index $FASTA || { echo "Failed to index reference genome"; exit 1; }

############################
# ALIGNMENT
############################

echo "Aligning reads to reference genome..."
bwa mem $FASTA $FASTQ1 $FASTQ2 > $SAM_FILE || { echo "Alignment failed"; exit 1; }

samtools flagstat $SAM_FILE > "${SRA_ID}_flagstat.txt"

############################
# SAM → BAM → SORT
############################

echo "Converting SAM to BAM..."
samtools view -Sb $SAM_FILE > $BAM_FILE

echo "Sorting BAM file..."
samtools sort $BAM_FILE -o $SORTED_BAM

echo "Indexing BAM file..."
samtools index $SORTED_BAM

############################
# VARIANT CALLING
############################

echo "Generating pileup..."
samtools mpileup -f $FASTA $SORTED_BAM > $PILEUP_FILE

echo "Calling variants..."
bcftools mpileup -f $FASTA $SORTED_BAM | \
bcftools call -mv -Oz -o $VCF_FILE

bcftools index $VCF_FILE

############################
# FILTER VARIANTS
############################

############################
# FILTER VARIANTS
############################

echo "Filtering variants (QUAL ≥ 30)..."

FILTERED_VCF="${SRA_ID}_filtered.recode.vcf"

vcftools --gzvcf "$VCF_FILE" \
         --minQ 30 \
         --recode \
         --recode-INFO-all \
         --out "${SRA_ID}_filtered" \
         2> vcftools_error.log

if [[ ! -f "$FILTERED_VCF" ]]; then
    echo "Filtered VCF not created. Exiting."
    cat vcftools_error.log
    exit 1
fi

echo "Filtered VCF created: $FILTERED_VCF"

############################
# VARIANT ANNOTATION (SnpEff)
############################

ANNOTATED_VCF="${SRA_ID}_annotated.vcf"

echo "Annotating variants using SnpEff..."

snpEff NC_045512.2 \
    "$FILTERED_VCF" > "$ANNOTATED_VCF" || { echo "SnpEff annotation failed"; exit 1; }

echo "Annotated VCF created: $ANNOTATED_VCF"
