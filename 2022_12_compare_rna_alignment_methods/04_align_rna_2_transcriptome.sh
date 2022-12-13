#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-12-13
# Updated... 2022-12-13
#
# Description:
# Align rnaseq reads back to their reference transcriptomes that are now
# just the primary alignments, high quality (>0.9 mapping rate), <20% mismatches
# of the canonical transcripts
# ------------------------------------------------------------------------------

## Align merged reads to reference transcriptomes -------------------------------------

# Create directories
mkdir -p /workdir/mbb262/b73/output/minimap_alignments

# Variables
N_THREADS=20
PROJ_DIR=/workdir/mbb262/b73
ALIGN_OUT=$PROJ_DIR/output/minimap_alignments
REF_TRANS=$PROJ_DIR/references/fa_transcriptomes/Zm-B73-REFERENCE-NAM-5.0_canonical_named.fa
FASTQC_TRIM=$PROJ_DIR/merged

# unit tests to b73 #1
# need to get sinle end code from other script
for i in $FASTQC_TRIM/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Aligning: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz
done
    /programs/minimap2-2.17/minimap2 -ax sr \
        -t $N_THREADS \
        --max-qlen 350 \
        $REF_TRANS \
        $FASTQC_TRIM/${SAMPLE}_1.fq.gz \
        $FASTQC_TRIM/${SAMPLE}_2.fq.gz \
        > $ALIGN_OUT/${SAMPLE}_minimap_to_${GENOME}.sam
done
