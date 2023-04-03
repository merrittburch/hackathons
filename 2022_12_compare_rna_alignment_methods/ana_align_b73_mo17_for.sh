#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-12-12
# Updated... 2023-01-05
#
# Description:
# Run B73 and Mo17 through the expression alignment pipeline for Ana to test ASE
# capabilities
# ------------------------------------------------------------------------------

# Make directories
mkdir -p /workdir/mbb262/ase
mkdir /workdir/mbb262/ase/raw_reads
mkdir /workdir/mbb262/ase/trimmed_reads


## Gather data -------------------------------------------------------------------
imeta qu -d sample_title like '2021RNAHybrids_GCLBp1A05' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_0801p2E12' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_GCLTp4D7' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'


# B73 are the first two, m017 the second two
cd /workdir/mbb262/ase/raw_reads
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R087/MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R087/MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R088/MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R088/MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_2.fq.gz

# additional hybrids 
#b73xmo17
imeta qu -d sample_title like '2021RNAHybrids_GCLTp4D7' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_GCLTp4A8' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_GCLTp4C9' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'

iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R119/MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R119/MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_2.fq.gz

# rename for simplicity
mv MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_1.fq.gz b73_1.fq.gz
mv MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_2.fq.gz b73_2.fq.gz
mv MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_1.fq.gz mo17_1.fq.gz
mv MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_2.fq.gz mo17_2.fq.gz
mv MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_1.fq.gz hybrid_b73_mo17_1.fq.gz
mv MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_2.fq.gz hybrid_b73_mo17_2.fq.gz


## Trim reads ----------------------------------------------------------------------

# Make somewhere for the summary files to go
mkdir -p /workdir/mbb262/ase/summary_trimming

# Variables
PROJ_DIR=/workdir/mbb262/ase
RAW_READ_DIR=/workdir/mbb262/ase/raw_reads
FASTQ_TRIM=$PROJ_DIR/trimmed_reads
SUMMARY_TRIM=$PROJ_DIR/summary_trimming
N_THREADS=30

# B73
/programs/fastp-0.23.2/fastp \
  --thread $N_THREADS \
  --in1 "${RAW_READ_DIR}/b73_1.fq.gz" \
  --in2 "${RAW_READ_DIR}/b73_2.fq.gz" \
  --merge \
  --merged_out "${FASTQ_TRIM}/b73.trim.fq.gz" \
  --detect_adapter_for_pe \
  --html "${SUMMARY_TRIM}/b73.fastp.html" \
  --report_title b73

# Mo17
/programs/fastp-0.23.2/fastp \
  --thread $N_THREADS \
  --in1 "${RAW_READ_DIR}/mo17_1.fq.gz" \
  --in2 "${RAW_READ_DIR}/mo17_2.fq.gz" \
  --merge \
  --merged_out "${FASTQ_TRIM}/mo17.trim.fq.gz" \
  --detect_adapter_for_pe \
  --html "${SUMMARY_TRIM}/mo17.fastp.html" \
  --report_title mo17

# Hybrid
/programs/fastp-0.23.2/fastp \
  --thread $N_THREADS \
  --in1 "${RAW_READ_DIR}/hybrid_b73_mo17_1.fq.gz" \
  --in2 "${RAW_READ_DIR}/hybrid_b73_mo17_2.fq.gz" \
  --merge \
  --merged_out "${FASTQ_TRIM}/hybrid_b73_mo17.trim.fq.gz" \
  --detect_adapter_for_pe \
  --html "${SUMMARY_TRIM}/hybrid_b73_mo17.fastp.html" \
  --report_title hybrid_b73_mo17


## Gather and format reference transcriptomes --------------------------------------------------------------

# Get trascripts
scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/ajs692/panand_data/gene_model_annotation/fastas/Zea_mays_v5_annotatedCDS.fa /workdir/mbb262/ase

# Format the long names
sed 's/:.*//' /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS.fa > /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString.fa

# Subsample fasta to specific genes (in this case: canonical)
cd /workdir/mbb262/ase/references
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts
/programs/samtools-1.15.1-r/bin/samtools faidx \
    /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString.fa \
    -r /workdir/mbb262/ase/references/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts \
    -o /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString_canonical.fa


# Gather genomes
mkdir -p /workdir/mbb262/ase/references
wget -P /workdir/mbb262/ase/references https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget -P /workdir/mbb262/ase/references https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Mo17-REFERENCE-YAN-1.0/Zm-Mo17-REFERENCE-YAN-1.0.fa.gz

# Find matching transcript sequence in each genome by minimapping all B73 CDS sequence to each NAM genome
mkdir /workdir/mbb262/ase/alignments
mkdir -p /workdir/mbb262/ase/references/nam_aligned_transcriptomes
N_THREADS=20

/programs/minimap2-2.17/minimap2 -ax splice --eqx \
    -t $N_THREADS \
    -I 6G \
    /workdir/mbb262/ase/references/Zm-B73-REFERENCE-NAM-5.0.fa \
    /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString_canonical.fa \
    > /workdir/mbb262/ase/references/nam_aligned_transcriptomes/Zm-B73-REFERENCE-NAM-5.0.sam

/programs/minimap2-2.17/minimap2 -ax splice --eqx \
    -t $N_THREADS \
    -I 6G \
    /workdir/mbb262/ase/references/Zm-Mo17-REFERENCE-YAN-1.0.fa \
    /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString_canonical.fa \
    > /workdir/mbb262/ase/references/nam_aligned_transcriptomes/Zm-Mo17-REFERENCE-YAN-1.0.sam


## Format reference transcriptomes ---------------------------------------------------
# (pull out sequence in each reference that aligns the best to the b73 cds seqeunce)

# Make directories
mkdir -p /workdir/mbb262/ase/references/bam_transcriptomes
mkdir -p /workdir/mbb262/ase/references/fa_transcriptomes
mkdir /workdir/mbb262/ase/references/bam_transcriptomes/primary

# Variables
SAM_TRANSCRIPTOME_DIR=/workdir/mbb262/ase/references/nam_aligned_transcriptomes
BAM_TRANSCRIPTOME_DIR=/workdir/mbb262/ase/references/bam_transcriptomes
FA_TRANSCRIPTOME_DIR=/workdir/mbb262/ase/references/fa_transcriptomes
PRIMARY_DIR=/workdir/mbb262/ase/references/bam_transcriptomes/primary

export PATH=/programs/samtools-1.15.1/bin:$PATH
export PATH=/programs/bedtools2-2.29.2/bin:$PATH

# Loop through all transcriptomes and format them
for i in $SAM_TRANSCRIPTOME_DIR/*.sam
do
    SAMPLE=$(basename ${i} .sam)

    # echo "SAM to BAM to FA to subsampled FA: " $SAM_TRANSCRIPTOME_DIR/${SAMPLE}.sam
    echo -e "SAM to BAM to FA to FA subset: ${SAMPLE}"

    # sam to bam
    /programs/samtools-1.15.1-r/bin/samtools view -bS -@ 15 $SAM_TRANSCRIPTOME_DIR/${SAMPLE}.sam | /programs/samtools-1.15.1-r/bin/samtools sort -@ 15 -o $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam

    # Index bam
    /programs/samtools-1.15.1-r/bin/samtools index $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam

    # Get primary alignments from bam
    /programs/samtools-1.15.1-r/bin/samtools view -e 'flag <= 16' $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam -o $PRIMARY_DIR/primary_${SAMPLE}.bam

    # bam to bed
    bedtools bamtobed -i $PRIMARY_DIR/primary_${SAMPLE}.bam > $PRIMARY_DIR/primary_${SAMPLE}.bed
done

# Go from bed to fasta
cd /workdir/mbb262/ase/references
gunzip *.fa.gz
find /workdir/mbb262/ase/references -maxdepth 1 -name '*.fa' > $PRIMARY_DIR/refGenome.list
find $PRIMARY_DIR/ -maxdepth 1 -name '*.bed' > $PRIMARY_DIR/genomeBeds.list
parallel --link "bedtools getfasta -name -s -fi {1} -bed {2/} -fo {1/.}.fasta" :::: refGenome.list :::: genomeBeds.list 


# Do last bit of formatting to read names
for i in $PRIMARY_DIR/*.fasta
do
    SAMPLE=$(basename ${i} .fasta)

    # Remove where in each reference the seqeunce was pulled from
    sed 's/:.*//' $PRIMARY_DIR/${SAMPLE}.fasta > $PRIMARY_DIR/${SAMPLE}_shortNames.fasta

    # Append genome name to fasta files
    sed "/^>/s/$/_${SAMPLE}/" $PRIMARY_DIR/${SAMPLE}_shortNames.fasta > $FA_TRANSCRIPTOME_DIR/${SAMPLE}_named.fa
done


## Create hybrid genomes ------------------------------------------

# create hybrid b73 and mo17
cd /workdir/mbb262/ase/references/fa_transcriptomes
cat Zm-B73-REFERENCE-NAM-5.0_named.fa Zm-Mo17-REFERENCE-YAN-1.0_named.fa > b73_mo17_combined_genome.fa


## Align back to genomes ------------------------------------------

# Create directories
mkdir -p /workdir/mbb262/ase/output/minimap_alignments

# Mapping parameters ---------------
N_THREADS=20
PROJ_DIR=/workdir/mbb262/ase
ALIGN_OUT=$PROJ_DIR/output/minimap_alignments
REF_TRANS=$PROJ_DIR/references/fa_transcriptomes
FASTQ_TRIM_MERGE_DIR=$PROJ_DIR/trimmed_reads

# Align b73 to b73
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_named.fa \
    $FASTQ_TRIM_MERGE_DIR/b73.trim.fq.gz > $ALIGN_OUT/b73_to_b73_minimap2.sam

# Align mo17 to mo17 
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    $REF_TRANS/Zm-Mo17-REFERENCE-YAN-1.0_named.fa \
    $FASTQ_TRIM_MERGE_DIR/mo17.trim.fq.gz > $ALIGN_OUT/mo17_to_mo17_minimap2.sam


# Align hybrid mo17xb73 to hybrid genome ---> forgot this one
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    $REF_TRANS/b73_mo17_combined_genome.fa \
    $FASTQ_TRIM_MERGE_DIR/hybrid_b73_mo17.trim.fq.gz > $ALIGN_OUT/hybridB73Mo17_to_jointB73Mo17_minimap2.sam

# Align hybrid mo17xb73 to b73
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_named.fa \
    $FASTQ_TRIM_MERGE_DIR/hybrid_b73_mo17.trim.fq.gz > $ALIGN_OUT/hybridB73Mo17_to_b73_minimap2.sam

# Align hybrid mo17xb73 to mo17
minimap2 -ax sr \
    --secondary=yes \
    -t $N_THREADS \
    $REF_TRANS/Zm-Mo17-REFERENCE-YAN-1.0_named.fa \
    $FASTQ_TRIM_MERGE_DIR/hybrid_b73_mo17.trim.fq.gz > $ALIGN_OUT/hybridB73Mo17_to_mo17_minimap2.sam


## Count mapped reads ------------------------------------------------------------

# Get transcript names
cd $REF_TRANS
grep -e ">" Zm-B73-REFERENCE-NAM-5.0_named.fa > temp.txt
sed 's/>//g' temp.txt > $ALIGN_OUT/all_b73_transcript_ids.txt

grep -e ">" Zm-Mo17-REFERENCE-YAN-1.0_named.fa > temp.txt
sed 's/>//g' temp.txt > $ALIGN_OUT/all_mo17_transcript_ids.txt

grep -e ">" b73_mo17_combined_genome.fa > temp.txt
sed 's/>//g' temp.txt > $ALIGN_OUT/all_b73_and_mo17_transcript_ids.txt


# Set up kotlin
wget https://github.com/JetBrains/kotlin/releases/download/v1.7.10/kotlin-compiler-1.7.10.zip
unzip kotlin-compiler-1.7.10.zip 
export PATH=/home/mbb262/bioinformatics/kotlinc_1.7.10/bin:$PATH
export _JAVA_OPTIONS=-Xmx50g

# make output dir
mkdir /workdir/mbb262/ase/output/counts


# Run Kotlin script 
# Modify in the script the name of the transcript file, the directory to the sams, and the out dir
/home/mbb262/bioinformatics/kotlinc/bin/kotlinc -script \
    /home/mbb262/git_projects/hackathons/2022_12_compare_rna_alignment_methods/05b_ana_count_noClip_sampleData.main.kts




