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

# Align reads to filtered down transcriptome
for i in $FASTQC_TRIM/*.fq.gz
do
    SAMPLE=$(basename ${i} .fq.gz)

    echo "Aligning: " ${SAMPLE}.fq.gz

    minimap2 -ax sr \
      -t $N_THREADS \
      $REF_TRANS \
      $FASTQC_TRIM/${SAMPLE}.fq.gz > $ALIGN_OUT/${SAMPLE}.sam

done
