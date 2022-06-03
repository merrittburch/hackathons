#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-06-03 
# Updated... 2022-06-03

# Description 
# Align RNAseq reads to refernce.consensus transcriptome using 
# minimap2 for shortread
# NOTE: some development could be done with parameters for optimization
# ---------------------------------------------------------------

# Mapping parameters & directories ---------------
mkdir /workdir/minimap2_alignments

export PATH=/programs/minimap2-2.17:$PATH
N_THREADS=25
FASTQC_TRIM_DIR=/workdir/expression/trimmed
GENOME_DIR=/path/to/ref/transcriptome
OUT_PATH=/workdir/minimap2_alignments

# Align paired end reads with each other
for i in $FASTQC_TRIM_DIR/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Aligning: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

    minimap2 -ax sr \
        $GENOME_DIR \
        $FASTQC_TRIM_DIR/${SAMPLE}_1.fq.gz \
        $FASTQC_TRIM_DIR/${SAMPLE}_2.fq.gz \
        -t $N_THREADS \
        > /workdir/minimap2_alignments/${SAMPLE}_alm.sam
done

