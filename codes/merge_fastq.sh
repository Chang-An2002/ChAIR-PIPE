#!/usr/bin/env bash
set -eo pipefail

# Usage: merge_fastq.sh PETSAMPLEID,ATACSAMPLEID,RNASAMPLEID,SAMPLEID
# Set INPATH to the directory containing the three input-library directories.
# Set OUTPATH to the sample's output directory, or edit the defaults below.
die() { echo "$*" >&2; exit 1; }
[[ $# -eq 1 ]] || die "Usage: bash $0 PETSAMPLEID,ATACSAMPLEID,RNASAMPLEID,SAMPLEID"
IFS=, read -r SAMPLEID1 SAMPLEID2 SAMPLEID3 SAMPLEID4 extra <<< "$1"
[[ -n "$SAMPLEID1" && -n "$SAMPLEID2" && -n "$SAMPLEID3" && -n "$SAMPLEID4" && -z "$extra" && "$1" == "$SAMPLEID1,$SAMPLEID2,$SAMPLEID3,$SAMPLEID4" ]] || die "Expected exactly four non-empty sample identifiers"

INPATH=${INPATH:-/path/to/your/input}
OUTPATH=${OUTPATH:-/path/to/your/output/${SAMPLEID4}}

# Collect R1 files in a fixed order and require the matching reads before writing.
# Accept both *_R1.fastq.gz and Illumina *_R1_001.fastq.gz names.
collect_files() {
    local sample=$1 reads=$2 file mate i
    [[ -d "$INPATH/$sample" ]] || die "Missing input directory: $INPATH/$sample"
    mapfile -d '' -t files < <(find "$INPATH/$sample" -type f \( -name '*_R1.fastq.gz' -o -name '*_R1_001.fastq.gz' \) -print0 | sort -z)
    [[ ${#files[@]} -gt 0 ]] || die "No R1 FASTQ files found in $INPATH/$sample"
    for file in "${files[@]}"; do
        for ((i=2; i<=reads; i++)); do
            if [[ "$file" == *_R1_001.fastq.gz ]]; then
                mate=${file%_R1_001.fastq.gz}_R${i}_001.fastq.gz
            else
                mate=${file%_R1.fastq.gz}_R${i}.fastq.gz
            fi
            [[ -f "$mate" ]] || die "Missing matching read: $mate"
        done
    done
}

collect_files "$SAMPLEID1" 3
dna_files=("${files[@]}")
collect_files "$SAMPLEID2" 3
dna_files+=("${files[@]}")
collect_files "$SAMPLEID3" 2
rna_files=("${files[@]}")

dna_prefix=$OUTPATH/$SAMPLEID4
rna_prefix=${dna_prefix}_GTAACATGCG-AGGTAACACT
outputs=()
for i in 1 2 3; do outputs+=("${dna_prefix}_S1_L001_R${i}_001.fastq.gz"); done
for i in 1 2; do outputs+=("${rna_prefix}_S1_L001_R${i}_001.fastq.gz"); done
for file in "${outputs[@]}"; do
    [[ ! -e "$file" ]] || die "Output already exists: $file; choose a new OUTPATH"
done
mkdir -p "$OUTPATH"

# Merge PET and ATAC DNA reads, preserving the same file order for R1/R2/R3.
for file in "${dna_files[@]}"; do
    for i in 1 2 3; do
        if [[ "$file" == *_R1_001.fastq.gz ]]; then
            mate=${file%_R1_001.fastq.gz}_R${i}_001.fastq.gz
        else
            mate=${file%_R1.fastq.gz}_R${i}.fastq.gz
        fi
        cat "$mate" >> "${dna_prefix}_S1_L001_R${i}_001.fastq.gz"
    done
done

# Merge RNA reads and truncate both sequence and quality, as in the protocol.
for file in "${rna_files[@]}"; do
    for i in 1 2; do
        if [[ "$file" == *_R1_001.fastq.gz ]]; then
            mate=${file%_R1_001.fastq.gz}_R${i}_001.fastq.gz
        else
            mate=${file%_R1.fastq.gz}_R${i}.fastq.gz
        fi
        length=28
        [[ "$i" -eq 1 ]] || length=90
        gzip -cd "$mate" |
            awk -v length_limit="$length" 'NR % 4 == 2 || NR % 4 == 0 {print substr($0, 1, length_limit); next} {print}' |
            gzip >> "${rna_prefix}_S1_L001_R${i}_001.fastq.gz"
    done
done

printf 'DNA prefix: %s\nRNA prefix: %s\n' "$dna_prefix" "$rna_prefix"
