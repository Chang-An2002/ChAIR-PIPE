#!/bin/bash

# Usage: merge_fastq.sh PETSAMPLEID,ATACSAMPLEID,RNASAMPLEID,SAMPLEID

SAMPLEID1=$(echo $1 | awk -F, '{print $1}')
SAMPLEID2=$(echo $1 | awk -F, '{print $2}')
SAMPLEID3=$(echo $1 | awk -F, '{print $3}')
SAMPLEID4=$(echo $1 | awk -F, '{print $4}')

INPATH="/path/to/your/input"
OUTPATH="/path/to/your/output/${SAMPLEID4}"

SUF=".fastq.gz"

# Merge DNA fastq files from multiple lanes
for SAMPLEID in $SAMPLEID1 $SAMPLEID2; do
    DNAFILES=$(find ${INPATH}/${SAMPLEID} -type f -name "*_R[1-3]${SUF}" | sed "s/_R[1-3]${SUF}//g" | sort -u)
    for file in $DNAFILES; do
        for i in {1..3}; do
            cat ${file}_R${i}${SUF} >>${OUTPATH}/${SAMPLEID4}_S1_L001_R${i}_001.fastq.gz
        done
    done
done

# Merge RNA fastq files from multiple lanes
RNAFILES=$(find ${INPATH}/${SAMPLEID3} -type f -name "*_R[1-2]${SUF}" | sed "s/_R[1-2]${SUF}//g" | sort -u)
for file in $RNAFILES; do
    zcat ${file}_R1${SUF} |
        awk '{if(NR % 4 == 1 || NR % 4 == 3){print $0}else{print substr($0, 0, 28)}}' |
        gzip >>${OUTPATH}/${SAMPLEID4}_GTAACATGCG-AGGTAACACT_S1_L001_R1_001.fastq.gz
    zcat ${file}_R2${SUF} |
        awk '{if(NR % 4 == 1 || NR % 4 == 3){print $0}else{print substr($0, 0, 90)}}' |
        gzip >>${OUTPATH}/${SAMPLEID4}_GTAACATGCG-AGGTAACACT_S1_L001_R2_001.fastq.gz
done
