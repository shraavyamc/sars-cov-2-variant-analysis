
#!/usr/bin/env bash
set -euo pipefail

#Configuration

SRA_ID="SRR29011221"
THREADS=4

REF_FASTA="NC_045512.2.fasta"
REF_NAME="NC_045512.2"

RESULTS_DIR="results"
DATA_DIR="data"
NEXTCLADE_DATASET="data/nextclade"

mkdir -p "$RESULTS_DIR" "$DATA_DIR"

# CHECK DEPENDENCIES


for tool in bwa samtools bcftools vcftools nextclade fasterq-dump prefetch bgzip tabix; do
    if ! command -v $tool &>/dev/null; then
        echo "ERROR: $tool not found in PATH"
        exit 1
    fi
done


# DOWNLOAD REFERENCE


if [[ ! -f "$REF_FASTA" ]]; then
    echo "Downloading SARS-CoV-2 reference genome..."
    wget -O "$REF_FASTA" \
      "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_045512.2&db=nuccore&report=fasta"
fi


# INDEX REFERENCE


bwa index "$REF_FASTA"
samtools faidx "$REF_FASTA"


# FETCH SRA


if [[ ! -f "${SRA_ID}.fastq" && ! -f "${SRA_ID}_1.fastq" ]]; then
    echo "Fetching SRA data..."
    prefetch "$SRA_ID"
    fasterq-dump "$SRA_ID" --split-files -O .
fi


# ALIGNMENT


echo "Aligning reads..."
bwa mem -t "$THREADS" "$REF_FASTA" \
    "${SRA_ID}_1.fastq" "${SRA_ID}_2.fastq" \
    > "${SRA_ID}.sam"

samtools view -bS "${SRA_ID}.sam" > "${SRA_ID}.bam"
samtools sort "${SRA_ID}.bam" -o "${SRA_ID}_sorted.bam"
samtools index "${SRA_ID}_sorted.bam"


# VARIANT CALLING

echo "Calling variants..."
bcftools mpileup -f "$REF_FASTA" "${SRA_ID}_sorted.bam" \
| bcftools call -mv -Oz -o "${SRA_ID}.vcf.gz"

tabix -p vcf "${SRA_ID}.vcf.gz"

# FILTER VARIANTS

FILTERED_VCF="${RESULTS_DIR}/${SRA_ID}_filtered.vcf"

echo "Filtering variants (QUAL ≥ 30)..."
vcftools --gzvcf "${SRA_ID}.vcf.gz" \
    --minQ 30 \
    --recode \
    --recode-INFO-all \
    --stdout > "$FILTERED_VCF"

if [[ ! -s "$FILTERED_VCF" ]]; then
    echo "Filtered VCF not created"
    exit 1
fi

# CONSENSUS GENOME

echo "Generating consensus genome..."
bgzip -c "$FILTERED_VCF" > "${FILTERED_VCF}.gz"
tabix -p vcf "${FILTERED_VCF}.gz"

bcftools consensus \
    -f "$REF_FASTA" \
    "${FILTERED_VCF}.gz" \
    > "${RESULTS_DIR}/${SRA_ID}_consensus.fasta"

# NEXTCLADE DATASET

if [[ ! -d "$NEXTCLADE_DATASET" ]]; then
    echo "Downloading Nextclade dataset..."
    nextclade dataset get \
        --name sars-cov-2 \
        --output-dir "$NEXTCLADE_DATASET"
fi

# NEXTCLADE ANNOTATION

echo "Running Nextclade..."
nextclade run \
    -D "$NEXTCLADE_DATASET" \
    "${RESULTS_DIR}/${SRA_ID}_consensus.fasta" \
    --output-tsv "${RESULTS_DIR}/${SRA_ID}_nextclade.tsv" \
    --output-json "${RESULTS_DIR}/${SRA_ID}_nextclade.json"

echo "Pipeline completed successfully"
