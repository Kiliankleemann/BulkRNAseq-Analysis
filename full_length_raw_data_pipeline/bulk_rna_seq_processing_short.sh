#!/bin/bash

# -------------------- VARIABLES - MANUAL ADUSTMENT! -------------------- #
# Directory for raw fastq files
RAW_DIR="./fastq_files"

# Fastq file endings - edit if necessary
R1_ENDING="_1.fq.gz"
R2_ENDING="_2.fq.gz"

# Number of threads for fastqc and cutadapt
THREADS=12

# Salmon transcriptome index (download for human or mouse! http://refgenomes.databio.org/)
SALMON_INDEX="/directory/to/salmon_transcriptome"

# Output directories
TRIMMED_DIR="$RAW_DIR/trimmed_reads"
QC_DIR="./QC"
TRANSCRIPT_QUANT_DIR="./transcript_quant"
# -------------------------------------------------------- #

# Move fastq files to fastq_files directory
mkdir -p "$RAW_DIR"
mv *.fq.gz "$RAW_DIR"

# Generate sample list
find "$RAW_DIR" -maxdepth 1 -type f -name "*.fq.gz" -exec basename {} \; | cut -d '_' -f1 | sort -u > sample_list.txt
echo "Samples found:"
cat sample_list.txt

# -------------------- FASTQC on raw reads -------------------- #
mkdir -p "$QC_DIR"
fastqc "$RAW_DIR"/*.fq.gz -t "$THREADS" -o "$QC_DIR"

# Quality overview
multiqc "$QC_DIR" -o "$QC_DIR"

# -------------------- Cutadapt for adapter trimming -------------------- #
mkdir -p "$TRIMMED_DIR"

while read sample; do
    cutadapt \
        --cores "$THREADS" \
        --minimum-length 15 \
        --poly-a \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o "$TRIMMED_DIR/${sample}.1.trimmed.fastq.gz" \
        -p "$TRIMMED_DIR/${sample}.2.trimmed.fastq.gz" \
        "$RAW_DIR/${sample}${R1_ENDING}" \
        "$RAW_DIR/${sample}${R2_ENDING}"
done < sample_list.txt

# -------------------- FASTQC on trimmed reads -------------------- #
mkdir -p "$QC_DIR/trimmed_reads"
fastqc -t "$THREADS" "$TRIMMED_DIR"/*.fastq.gz
multiqc "$TRIMMED_DIR" -o "$QC_DIR/trimmed_reads"

# -------------------- Salmon quantification -------------------- #
mkdir -p "$TRANSCRIPT_QUANT_DIR"

while read sample; do
    salmon quant \
        -i "$SALMON_INDEX" \
        -l A \
        -1 "$TRIMMED_DIR/${sample}.1.trimmed.fastq.gz" \
        -2 "$TRIMMED_DIR/${sample}.2.trimmed.fastq.gz" \
        --validateMappings \
        -o "$TRANSCRIPT_QUANT_DIR/${sample}_quant" \
        -p 8
done < sample_list.txt

echo "All samples processed successfully."
