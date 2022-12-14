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
# Count reads from sam files using custom kotlin script
# ------------------------------------------------------------------------------

# Set up kotlin
wget https://github.com/JetBrains/kotlin/releases/download/v1.7.10/kotlin-compiler-1.7.10.zip
unzip kotlin-compiler-1.7.10.zip 
export PATH=/workdir/mbb262/kotlinc_1.7.10/bin:$PATH
export _JAVA_OPTIONS=-Xmx50g

# make output dir
mkdir /workdir/mbb262/b73/output/counts

# Get all transcript ids
cd /workdir/mbb262/b73/references/fa_transcriptomes
grep -e ">" Zm-B73-REFERENCE-NAM-5.0_canonical_named.fa > temp.txt
sed 's/>//g' temp.txt > all_b73_transcript_ids.txt
mv all_b73_transcript_ids.txt /workdir/mbb262/b73/output/counts


# Run Kotlin script 
/workdir/mbb262/kotlinc/bin/kotlinc -script \
    /home/mbb262/git_projects/hackathons/2022_12_compare_rna_alignment_methods/05b_count.main.kts
