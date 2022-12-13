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
# Align B73 RNAseq reads from NAM and Anju's experiment using my experimental pipeline  
# ------------------------------------------------------------------------------


## Gather and fish out transcript sequences from reference genomes -------------

# Get trascripts
scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/ajs692/panand_data/gene_model_annotation/fastas/Zea_mays_v5_annotatedCDS.fa /workdir/mbb262/b73

# Format the long names
sed 's/:.*//' /workdir/mbb262/b73/Zea_mays_v5_annotatedCDS.fa > /workdir/mbb262/b73/Zea_mays_v5_annotatedCDS_shortenedString.fa

# Gather genomes
mkdir -p /workdir/mbb262/b73/references
wget -P /workdir/mbb262/b73/references https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz

# Find matching B73 transcript sequence in each genome
mkdir /workdir/mbb262/b73/alignments
mkdir -p /workdir/mbb262/b73/references/nam_aligned_transcriptomes
mkdir -p /workdir/mbb262/b73/alignments/nam_aligned_transcriptomes

N_THREADS=20
/programs/minimap2-2.17/minimap2 -ax splice --eqx \
    -t $N_THREADS \
    -I 6G \
    /workdir/mbb262/b73/references/Zm-B73-REFERENCE-NAM-5.0.fa.gz \
    /workdir/mbb262/b73/Zea_mays_v5_annotatedCDS_shortenedString.fa \
    > /workdir/mbb262/b73/references/nam_aligned_transcriptomes/Zm-B73-REFERENCE-NAM-5.0.sam


## Section picks the primary alignment thats >90% aligned -----------------------

# Gather primary alignment of each transcript that aligns >.9, export as sam file

# Set up kotlin (if not already set up)
cd /workdir/mbb262
wget https://github.com/JetBrains/kotlin/releases/download/v1.7.10/kotlin-compiler-1.7.10.zip
unzip kotlin-compiler-1.7.10.zip 
export PATH=/workdir/mbb262/kotlinc_1.7.10/bin:$PATH
export _JAVA_OPTIONS=-Xmx50g


# Call Kotlin script 
/workdir/mbb262/kotlinc/bin/kotlinc -script \
    /home/mbb262/git_projects/hackathons/2022_12_compare_rna_alignment_methods/03b_select_sam_primary.main.kts \
    /workdir/mbb262/b73/references/nam_aligned_transcriptomes/Zm-B73-REFERENCE-NAM-5.0.sam 0.9 0.05

# Count the number of things filtered
wc -l Zm-B73-REFERENCE-NAM-5.0.sam
wc -l Zm-B73-REFERENCE-NAM-5.0_filtered.sam

## Format sam files and get actual fasta sequence to align to -------------------

# Get canonical transcripts: to be replaced with better gene choices later
cd /workdir/mbb262/b73/references
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts

# Make directories
mkdir -p /workdir/mbb262/b73/references/bam_transcriptomes
mkdir -p /workdir/mbb262/b73/references/fa_transcriptomes

# Variables
SAM_TRANSCRIPTOME_DIR=/workdir/mbb262/b73/references/nam_aligned_transcriptomes
BAM_TRANSCRIPTOME_DIR=/workdir/mbb262/b73/references/bam_transcriptomes
FA_TRANSCRIPTOME_DIR=/workdir/mbb262/b73/references/fa_transcriptomes

# Loop through all transcriptomes and format them
for i in $SAM_TRANSCRIPTOME_DIR/*filtered.sam
do
    SAMPLE=$(basename ${i} _filtered.sam)

    # echo "SAM to BAM to FA to subsampled FA: " $SAM_TRANSCRIPTOME_DIR/${SAMPLE}.sam
    echo -e "SAM to BAM to FA to FA subset: ${SAMPLE}"

    # sam to bam
    /programs/samtools-1.15.1-r/bin/samtools view -bS $SAM_TRANSCRIPTOME_DIR/${SAMPLE}.sam > $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam

    # bam to fa
    /programs/samtools-1.15.1-r/bin/samtools fasta $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam > $FA_TRANSCRIPTOME_DIR/${SAMPLE}.fa

    # Subsample to specific gene list (in this case: canonical)
    /programs/samtools-1.15.1-r/bin/samtools faidx $FA_TRANSCRIPTOME_DIR/${SAMPLE}.fa \
        -r /workdir/mbb262/b73/references/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts \
        -o $FA_TRANSCRIPTOME_DIR/${SAMPLE}_canonical.fa

    # Append genome name to fasta files
    sed "/^>/s/$/_${SAMPLE}/" $FA_TRANSCRIPTOME_DIR/${SAMPLE}_canonical.fa > $FA_TRANSCRIPTOME_DIR/${SAMPLE}_canonical_named.fa
done

# Clean up fa_transcriptomes directory
mkdir /workdir/mbb262/b73/references/fai_files
mkdir /workdir/mbb262/b73/references/unnamed_transcripts

mv $FA_TRANSCRIPTOME_DIR/*fai /workdir/mbb262/b73/references/fai_filess





