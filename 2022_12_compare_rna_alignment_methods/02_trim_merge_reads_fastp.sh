#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-12-12
# Updated... 2022-12-12
#
# Description:
# Trim reads using fastp
# ------------------------------------------------------------------------------

# Variables
mkdir -p /workdir/tf259/hackathon_dec2022/reads/b73/trimmed_reads
PROJ_DIR=/workdir/tf259/hackathon_dec2022/reads/b73
RAW_READ_DIR=/workdir/tf259/hackathon_dec2022/reads/b73
FASTQ_TRIM=$PROJ_DIR/trimmed_reads

# Trim reads
for i in $RAW_READ_DIR/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Trimming: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

    /programs/fastp-0.23.2/fastp \
      --in1 "${RAW_READ_DIR}/${SAMPLE}_1.fq.gz" \
      --in2 "${RAW_READ_DIR}/${SAMPLE}_2.fq.gz" \
      --out1 "${FASTQ_TRIM}/${SAMPLE}_1_trim.fq.gz" \
      --out2 "${FASTQ_TRIM}/${SAMPLE}_2_trim.fq.gz" \
      --detect_adapter_for_pe \
      --json "${FASTQ_TRIM}/${SAMPLE}.fastp.json" \
      --html "${FASTQ_TRIM}/${SAMPLE}.fastp.html" \
      --report_title ${SAMPLE}
done


# Make directories
mkdir -p /workdir/mbb262/b73
mkdir /workdir/mbb262/b73/trimmed_reads
cp /workdir/tf259/hackathon_dec2022/reads/b73/trimmed_reads/*.fq.gz /workdir/mbb262/b73/trimmed_reads
cp /workdir/tf259/hackathon_dec2022/reads/NAM/trimmed_reads/*.fq.gz /workdir/mbb262/b73/trimmed_reads
cp /workdir/tf259/hackathon_dec2022/reads/simulated_cdna/*.fq.gz /workdir/mbb262/b73/trimmed_reads


##  Merge paired reads with BBmerge --------------------------------------------

# For reads all in the main folder
# Variables
mkdir /workdir/mbb262/b73/merged_reports
mkdir /workdir/mbb262/b73/merged
PROJ_DIR=/workdir/mbb262/b73
FASTQ_TRIM=$PROJ_DIR/trimmed_reads
MERGE_DIR=$PROJ_DIR/merged
OUT_DIR_INSERT_SIZE=$PROJ_DIR/reads/merged_reports

for i in $FASTQ_TRIM/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Trimming: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

    /programs/bbmap-38.96/bbmerge.sh \
        in1=$FASTQ_TRIM/${SAMPLE}_1.fq.gz  \
        in2=$FASTQ_TRIM/${SAMPLE}_2.fq.gz \
        out=$MERGE_DIR/${SAMPLE}_trimmed_bbmerge.fq.gz \
        ihist=$OUT_DIR_INSERT_SIZE/${SAMPLE}_trimmed_bbmerge.txt
done

# For other reads
for i in $FASTQ_TRIM/*_1.trim.fq.gz
do
    SAMPLE=$(basename ${i} _1.trim.fq.gz)

    echo "Trimming: " ${SAMPLE}_1.trim.fq.gz ${SAMPLE}_2.trim.fq.gz
print
    /programs/bbmap-38.96/bbmerge.sh \
        in1=$FASTQ_TRIM/${SAMPLE}_1.trim.fq.gz  \
        in2=$FASTQ_TRIM/${SAMPLE}_2.trim.fq.gz \
        out=$MERGE_DIR/${SAMPLE}_trimmed_bbmerge.fq.gz \
        ihist=$OUT_DIR_INSERT_SIZE/${SAMPLE}_trimmed_bbmerge.txt
done


# For simulated reads ------------------------------------------
FASTQ_TRIM=/workdir/tf259/hackathon_dec2022/reads/simulated_cdna_B73_v5
/programs/bbmap-38.96/bbmerge.sh \
        in1=$FASTQ_TRIM/sample_01_1.fq.gz  \
        in2=$FASTQ_TRIM/sample_01_2.fq.gz \
        out=$MERGE_DIR/sample_01_trimmed_bbmerge.fq.gz \
        ihist=$OUT_DIR_INSERT_SIZE/sample_01_trimmed_bbmerge.txt

/programs/bbmap-38.96/bbmerge.sh \
        in1=$FASTQ_TRIM/sample_02_1.fq.gz  \
        in2=$FASTQ_TRIM/sample_02_2.fq.gz \
        out=$MERGE_DIR/sample_02_trimmed_bbmerge.fq.gz \
        ihist=$OUT_DIR_INSERT_SIZE/sample_02_trimmed_bbmerge.txt

/programs/bbmap-38.96/bbmerge.sh \
        in1=$FASTQ_TRIM/sample_03_1.fq.gz  \
        in2=$FASTQ_TRIM/sample_03_2.fq.gz \
        out=$MERGE_DIR/sample_03_trimmed_bbmerge.fq.gz \
        ihist=$OUT_DIR_INSERT_SIZE/sample_03_trimmed_bbmerge.txt


# For TERRA_MEPP reads that joined the pipeline much later --------------------
PROJ_DIR=/workdir/mbb262/b73
FASTQ_TRIM=/workdir/tf259/hackathon_dec2022/reads/TERRA-MEPP
MERGE_DIR=$PROJ_DIR/merged
OUT_DIR_INSERT_SIZE=$PROJ_DIR/reads/merged_reports

for i in $FASTQ_TRIM/*_1.trim.fq.gz
do
    SAMPLE=$(basename ${i} _1.trim.fq.gz)

    echo "Trimming: " ${SAMPLE}_1.trim.fq.gz ${SAMPLE}_2.trim.fq.gz

    /programs/bbmap-38.96/bbmerge.sh \
        in1=$FASTQ_TRIM/${SAMPLE}_1.trim.fq.gz  \
        in2=$FASTQ_TRIM/${SAMPLE}_2.trim.fq.gz \
        out=$MERGE_DIR/${SAMPLE}_trimmed_bbmerge.fq.gz \
        ihist=$OUT_DIR_INSERT_SIZE/${SAMPLE}_trimmed_bbmerge.txt
done


# For Karl's 282 reads that joined the pipeline much later ---------------------

# Don't need to merge anything - they're good as is, copy over to merged read dir
cp /workdir/tf259/hackathon_dec2022/reads/maize282/trimmed_reads/*.trim.fq.gz $MERGE_DIR









