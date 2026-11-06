#!/bin/bash

set -eo pipefail

TOTAL_CORES=$(sysctl -n hw.ncpu)
THREADS=$(( TOTAL_CORES > 1 ? TOTAL_CORES - 1 : 1 ))

echo "Using ${THREADS} threads for processing."

# Path setting
stem_path="."
sample_list="sample_list.txt"
reference_name="Paguma_larvata_cytb.fas"

# Conda
eval "$(conda shell.bash hook)"
conda activate Cytb_mapping

cd "${stem_path}"
mkdir -p cleaned_reads mapping_results Cytb_consensus

# Loop processing for each sample
while IFS= read -r sample_ID || [ -n "$sample_ID" ]; do
    [ -z "$sample_ID" ] && continue

    echo "=========================================="
    echo "Processing sample: ${sample_ID}"
    echo "=========================================="

    # 1. Trim and quality filter NGS reads
    trim_galore --paired \
        --cores "${THREADS}" \
        "/Fujita_Civet cell RNAseq/2025.10.16_CQ627_RNA-seq_6G_24sample/BMK250729-CQ627-ZX01-0101/BMK_DATA_20250926152634_1/Data/rawdata/${sample_ID}_1.fq.gz" \
        "/Fujita_Civet cell RNAseq/2025.10.16_CQ627_RNA-seq_6G_24sample/BMK250729-CQ627-ZX01-0101/BMK_DATA_20250926152634_1/Data/rawdata/${sample_ID}_2.fq.gz" \
        --output_dir cleaned_reads

    # 2 & 3. Mapping, Sorting, and Direct BAM output
    minimap2 -ax sr -t "${THREADS}" "${reference_name}" \
        "cleaned_reads/${sample_ID}_1_val_1.fq.gz" \
        "cleaned_reads/${sample_ID}_2_val_2.fq.gz" \
        | samtools sort -@"${THREADS}" -o "mapping_results/${sample_ID}.sorted.bam"

    # Index BAM file
    samtools index -@"${THREADS}" "mapping_results/${sample_ID}.sorted.bam"

    # 4. Make consensus
    out_fa="Cytb_consensus/${sample_ID}.cons.fas"
    tmp_fa="${out_fa}.tmp"

    samtools consensus -a --show-ins yes --show-del yes "mapping_results/${sample_ID}.sorted.bam" -o "${tmp_fa}"
    awk -v id="${sample_ID}" 'NR==1{print ">"id; next} {print}' "${tmp_fa}" > "${out_fa}"
    rm -f "${tmp_fa}"

    echo "Finished sample: ${sample_ID}"

done < "${sample_list}"

echo "All samples finished completely!"
