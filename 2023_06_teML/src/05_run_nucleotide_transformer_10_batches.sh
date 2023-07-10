#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-06-08
# Updated... 2023-06-08
#
# Description:
# Create a loop that iterates through all TE files, runs nucleotide transformer,
# outputs results to file
# Combine all results within ridge regression or lasso model script
# ------------------------------------------------------------------------------

# Go to all files
TE_SEQ_DIR=/workdir/mbb262/te
cd $TE_SEQ_DIR

# Create a loop
for i in *_te_sequence_with_walley_expression.txt
do
    echo "Running: " ${i}

    python /workdir/mbb262/te_hugging_face.py ${i} /workdir/mbb262/te/te_nucleotideTransformer_results/${i}_result_output_chunk.txt
done