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


# Run Kotlin script 
/workdir/mbb262/kotlinc/bin/kotlinc -script \
    /home/mbb262/git_projects/te_ase_nam/src/04b_count_rnaseq_reads_minimap2.main.kts
