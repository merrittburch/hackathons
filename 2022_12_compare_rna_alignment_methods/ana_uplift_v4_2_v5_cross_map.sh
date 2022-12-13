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
# Uplift hapmap 321 snps from v4 to v5
# ------------------------------------------------------------------------------

# Make directories
mkdir -p /workdir/mbb262/snps/v4

# Get SNPs from cbsu
scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/genotypes/nam/imputed/nam_and_ibm_imputed_taxa/*vcf.gz /workdir/mbb262/snps/v4

# Get cross-map file from maizegdb and v5 genome
cd /workdir/mbb262/snps
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/chain_files/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain 
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz

# Use crossmap to uplift files
# Set path
export PYTHONPATH=/programs/CrossMap-0.6.1/lib64/python3.9/site-packages:/programs/CrossMap-0.6.1/lib/python3.9/site-packages
export PATH=/programs/CrossMap-0.6.1/bin:$PATH

# Create a loop that iterates through each chromosome & uplifts for INBRED data
mkdir /workdir/mbb262/snps/v5
for CHROM in {1..10}
do
  echo "I am on chr ${CHROM}"

  CrossMap.py vcf \
    ./B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain \
    /workdir/mbb262/snps/v4/nam_ibm_imputed_chr${CHROM}.vcf.gz \
    ./Zm-B73-REFERENCE-NAM-5.0.fa \
    /workdir/mbb262/snps/v5/nam_ibm_imputed_v5_chr${CHROM}.vcf

  bgzip --threads 20 /workdir/mbb262/snps/v5/nam_ibm_imputed_v5_chr${CHROM}.vcf

  echo "I just finished chr ${CHROM}"
done
