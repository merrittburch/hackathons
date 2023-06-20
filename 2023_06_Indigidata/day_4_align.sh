#!/bin/bash

#SBATCH --job-name=mer_23_candida
#SBATCH --account=PAS2510
#SBATCH --time=6:00:00
#SBATCH --nodes=1 --ntasks-per-node=28
#SBATCH --mail-user=mbb262@cornell.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL


# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-06-15
# Updated... 2023-06-15
#
# Description:
# - Map reads to the Candidia genome
# ------------------------------------------------------------------------------

# Load correct modules
module load bowtie2

# Change into my own directory
cd /fs/scratch/PAS2510/Indigidata_2023/merritt/genome

# Align reads to the genome
bowtie2 -x candida_albicans_alignment \
    -q -1 \
    /fs/scratch/PAS2510/Indigidata_2023/merritt/Candida_raw_files/trimmed_reads/paired_trimmed_SC5314_CRISPR_MAY1376_R1.fastq.gz -2 \
    /fs/scratch/PAS2510/Indigidata_2023/merritt/Candida_raw_files/trimmed_reads/paired_trimmed_SC5314_CRISPR_MAY1376_R2.fastq.gz -S \
    /fs/scratch/PAS2510/Indigidata_2023/merritt/genome/SC5314_CRISPR_MAY1376.sam
