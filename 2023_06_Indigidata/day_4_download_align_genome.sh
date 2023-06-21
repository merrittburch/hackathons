# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-06-15
# Updated... 2023-06-15
#
# Description:
# - Download Candidia albicans genome, index it
# ------------------------------------------------------------------------------


# Make directory
mkdir /fs/scratch/PAS2510/Indigidata_2023/merritt/genome
cd /fs/scratch/PAS2510/Indigidata_2023/merritt/genome

# Download candidia albicans sequence
# Genome link: http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/
wget http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_A22_current_chromosomes.fasta.gz

# Start conda environment
conda activate indigidata

# Load samtools
module load samtools

# Index samtools
gunzip  C_albicans_SC5314_A22_current_chromosomes.fasta.gz
samtools faidx ./C_albicans_SC5314_A22_current_chromosomes.fasta

# Build index file from the reference genome
module load bowtie2
bowtie2-build -f ./C_albicans_SC5314_A22_current_chromosomes.fasta \
    ./candida_albicans_alignment

# Align reads in the other script

# Convert to bam file
cd /fs/scratch/PAS2510/Indigidata_2023/merritt/genome
samtools view -bh SC5314_CRISPR_MAY1376.sam > SC5314_CRISPR_MAY1376.bam

# Sort and index bam files
samtools sort SC5314_CRISPR_MAY1376.bam > SC5314_CRISPR_MAY1376_sorted.bam
samtools index SC5314_CRISPR_MAY1376_sorted.bam

