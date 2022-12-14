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
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R119/MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R119/MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_2.fq.gz

# additional hybrids 
#b73xmo17
imeta qu -d sample_title like '2021RNAHybrids_GCLTp4D7' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_GCLTp4A8' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_GCLTp4C9' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'

# B73xky21
imeta qu -d sample_title like '2021RNAHybrids_GCGPp4A6' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_0716p4C6' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_0723p4G6' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'

iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R119/MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R119/MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R124/MS21R124_CKDL210018333-2a-AK690-7UDI304_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R124/MS21R124_CKDL210018333-2a-AK690-7UDI304_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R134/MS21R134_CKDL210018333-2a-AK850-AK2430_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R134/MS21R134_CKDL210018333-2a-AK850-AK2430_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R233/MS21R233_CKDL210018333-2a-AK17213-7UDI233_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R233/MS21R233_CKDL210018333-2a-AK17213-7UDI233_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R295/MS21R295_CKDL210018333-2a-AK6662-AK30512_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R295/MS21R295_CKDL210018333-2a-AK6662-AK30512_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R362/MS21R362_CKDL210018333-2a-AK30484-7UDI1611_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R362/MS21R362_CKDL210018333-2a-AK30484-7UDI1611_HH5V7DSX2_L1_2.fq.gz


# rename for simplicity
mv MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_1.fq.gz b73_1.fq.gz
mv MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_2.fq.gz b73_2.fq.gz
mv MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_1.fq.gz mo17_1.fq.gz
mv MS21R088_CKDL210018333-2a-AK1952-AK2044_HH5V7DSX2_L1_2.fq.gz mo17_2.fq.gz
mv MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_1.fq.gz hybrid_b73_mo17_1.fq.gz
mv MS21R119_CKDL210018333-2a-GE09-AK1954_HH5V7DSX2_L1_2.fq.gz hybrid_b73_mo17_2.fq.gz

## Trim reads ----------------------------------------------------------------------

# Variables
PROJ_DIR=/workdir/mbb262/ase
RAW_READ_DIR=/workdir/mbb262/ase/raw_reads
FASTQ_TRIM=$PROJ_DIR/trimmed_reads

# B73
/programs/fastp-0.23.2/fastp \
  --in1 "${RAW_READ_DIR}/b73_1.fq.gz" \
  --in2 "${RAW_READ_DIR}/b73_2.fq.gz" \
  --out1 "${FASTQ_TRIM}/b73_1.trim.fq.gz" \
  --out2 "${FASTQ_TRIM}/b73_2.trim.fq.gz" \
  --detect_adapter_for_pe \
  --json "${FASTQ_TRIM}/b73.fastp.json" \
  --html "${FASTQ_TRIM}/b73.fastp.html" \
  --report_title b73

# Mo17
/programs/fastp-0.23.2/fastp \
  --in1 "${RAW_READ_DIR}/mo17_1.fq.gz" \
  --in2 "${RAW_READ_DIR}/mo17_2.fq.gz" \
  --out1 "${FASTQ_TRIM}/mo17_1.trim.fq.gz" \
  --out2 "${FASTQ_TRIM}/mo17_2.trim.fq.gz" \
  --detect_adapter_for_pe \
  --json "${FASTQ_TRIM}/mo17.fastp.json" \
  --html "${FASTQ_TRIM}/mo17.fastp.html" \
  --report_title mo17

# Hybrid
/programs/fastp-0.23.2/fastp \
  --in1 "${RAW_READ_DIR}/hybrid_b73_mo17_1.fq.gz" \
  --in2 "${RAW_READ_DIR}/hybrid_b73_mo17_2.fq.gz" \
  --out1 "${FASTQ_TRIM}/hybrid_b73_mo17_1.trim.fq.gz" \
  --out2 "${FASTQ_TRIM}/hybrid_b73_mo17_2.trim.fq.gz" \
  --detect_adapter_for_pe \
  --json "${FASTQ_TRIM}/hybrid_b73_mo17.fastp.json" \
  --html "${FASTQ_TRIM}/hybrid_b73_mo17.fastp.html" \
  --report_title hybrid_b73_mo17


##  Merge paired reads with BBmerge -------------------------------------------------------------------

# Variables
mkdir /workdir/mbb262/ase/merged_reports
mkdir /workdir/mbb262/ase/merged
PROJ_DIR=/workdir/mbb262/ase
FASTQ_TRIM=$PROJ_DIR/trimmed_reads
MERGE_DIR=$PROJ_DIR/merged
OUT_DIR_INSERT_SIZE=$PROJ_DIR/reads/merged_reports

for i in $FASTQ_TRIM/*_1.trim.fq.gz
do
    SAMPLE=$(basename ${i} _1.trim.fq.gz)

    echo "Trimming: " ${SAMPLE}_1.trim.fq.gz ${SAMPLE}_2.trim.fq.gz

    /programs/bbmap-38.96/bbmerge.sh \
        in1=$FASTQ_TRIM/${SAMPLE}_1.trim.fq.gz  \
        in2=$FASTQ_TRIM/${SAMPLE}_2.trim.fq.gz \
        out=$MERGE_DIR/${SAMPLE}_trimmed_bbmerge.fq.gz \
        ihist=$OUT_DIR_INSERT_SIZE/${SAMPLE}_trimmed_bbmerge.txt
done

# Do for left out samples because I forgot them
/programs/bbmap-38.96/bbmerge.sh \
        in1=$FASTQ_TRIM/hybrid_b73_mo17_1.trim.fq.gz  \
        in2=$FASTQ_TRIM/hybrid_b73_mo17_2.trim.fq.gz \
        out=$MERGE_DIR/hybrid_b73_mo17_trimmed_bbmerge.fq.gz \
        ihist=$OUT_DIR_INSERT_SIZE/hybrid_b73_mo17_trimmed_bbmerge.txt


## Gather and format reference genomes --------------------------------------------------------------

# Get trascripts
scp mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/ajs692/panand_data/gene_model_annotation/fastas/Zea_mays_v5_annotatedCDS.fa /workdir/mbb262/ase

# Format the long names
sed 's/:.*//' /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS.fa > /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString.fa

# Gather genomes
mkdir -p /workdir/mbb262/ase/references
wget -P /workdir/mbb262/ase/references https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget -P /workdir/mbb262/ase/references https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Mo17-REFERENCE-YAN-1.0/Zm-Mo17-REFERENCE-YAN-1.0.fa.gz

# Find matching transcript sequence in each genome
mkdir /workdir/mbb262/ase/alignments
N_THREADS=20
/programs/minimap2-2.17/minimap2 -ax splice --eqx \
    -t $N_THREADS \
    -I 6G \
    /workdir/mbb262/ase/references/Zm-B73-REFERENCE-NAM-5.0.fa.gz \
    /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString.fa \
    > /workdir/mbb262/ase/alignments/Zm-B73-REFERENCE-NAM-5.0.sam

/programs/minimap2-2.17/minimap2 -ax splice --eqx \
    -t $N_THREADS \
    -I 6G \
    /workdir/mbb262/ase/references/Zm-Mo17-REFERENCE-YAN-1.0.fa.gz \
    /workdir/mbb262/ase/Zea_mays_v5_annotatedCDS_shortenedString.fa \
    > /workdir/mbb262/ase/alignments/Zm-Mo17-REFERENCE-YAN-1.0.sam

mkdir -p /workdir/mbb262/ase/references/nam_aligned_transcriptomes
mv /workdir/mbb262/ase/alignments/*.sam /workdir/mbb262/ase/references/nam_aligned_transcriptomes

# Get canonical transcripts: to be replaced with better gene choices later
cd /workdir/mbb262/ase/references
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts

# Make directories
mkdir -p /workdir/mbb262/ase/references/bam_transcriptomes
mkdir -p /workdir/mbb262/ase/references/fa_transcriptomes

# Variables
SAM_TRANSCRIPTOME_DIR=/workdir/mbb262/ase/references/nam_aligned_transcriptomes
BAM_TRANSCRIPTOME_DIR=/workdir/mbb262/ase/references/bam_transcriptomes
FA_TRANSCRIPTOME_DIR=/workdir/mbb262/ase/references/fa_transcriptomes
export PATH=/programs/samtools-1.15.1/bin:$PATH

# Loop through all transcriptomes and format them
for i in $SAM_TRANSCRIPTOME_DIR/*.sam
do
    SAMPLE=$(basename ${i} .sam)

    # echo "SAM to BAM to FA to subsampled FA: " $SAM_TRANSCRIPTOME_DIR/${SAMPLE}.sam
    echo -e "SAM to BAM to FA to FA subset: ${SAMPLE}"

    # sam to bam
    /programs/samtools-1.15.1-r/bin/samtools view -bS $SAM_TRANSCRIPTOME_DIR/${SAMPLE}.sam > $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam

    # bam to fa
    /programs/samtools-1.15.1-r/bin/samtools fasta $BAM_TRANSCRIPTOME_DIR/${SAMPLE}.bam > $FA_TRANSCRIPTOME_DIR/${SAMPLE}.fa
 
    # Subsample to specific gene list (in this case: canonical)
    /programs/samtools-1.15.1-r/bin/samtools faidx $FA_TRANSCRIPTOME_DIR/${SAMPLE}.fa \
        -r /workdir/mbb262/ase/references/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts \
        -o $FA_TRANSCRIPTOME_DIR/${SAMPLE}_canonical.fa

    # Append genome name to fasta files
    sed "/^>/s/$/_${SAMPLE}/" $FA_TRANSCRIPTOME_DIR/${SAMPLE}_canonical.fa > $FA_TRANSCRIPTOME_DIR/${SAMPLE}_canonical_named.fa
done

# Clean up fa_transcriptomes directory
mkdir /workdir/mbb262/ase/references/fai_files
mkdir /workdir/mbb262/ase/references/unnamed_transcripts

mv /workdir/mbb262/ase/references/fa_transcriptomes/*fai /workdir/mbb262/ase/references/fai_files
mv /workdir/mbb262/ase/references/fa_transcriptomes/*_canonical.fa /workdir/mbb262/ase/references/unnamed_transcripts
mv /workdir/mbb262/ase/references/fa_transcriptomes/*.0.fa /workdir/mbb262/ase/references/unnamed_transcripts


# Create hybrid genomes ------------------------------------------

# create hybrid b73 and mo17
cd /workdir/mbb262/ase/references/fa_transcriptomes
cat Zm-B73-REFERENCE-NAM-5.0_canonical_named.fa Zm-Mo17-REFERENCE-YAN-1.0_canonical_named.fa > b73_mo17_combined_genome.fa


# Align back to genomes ------------------------------------------

# Create directories
mkdir -p /workdir/mbb262/ase/output/minimap_alignments

# Mapping parameters ---------------
N_THREADS=20
PROJ_DIR=/workdir/mbb262/ase
ALIGN_OUT=$PROJ_DIR/output/minimap_alignments
REF_TRANS=$PROJ_DIR/references/fa_transcriptomes
FASTQ_TRIM_MERGE_DIR=$PROJ_DIR/merged

# Align b73 to b73
minimap2 -ax sr \
    -t $N_THREADS \
    $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_canonical_named.fa \
    $FASTQ_TRIM_MERGE_DIR/b73_trimmed_bbmerge.fq.gz > $ALIGN_OUT/b73_to_b73_minimap2.sam

# Align mo17 to mo17 
minimap2 -ax sr \
    -t $N_THREADS \
    $REF_TRANS/Zm-Mo17-REFERENCE-YAN-1.0_canonical_named.fa \
    $FASTQ_TRIM_MERGE_DIR/mo17_trimmed_bbmerge.fq.gz > $ALIGN_OUT/mo17_to_mo17_minimap2.sam


# Align hybrid mo17xb73 to hybrid genome ---> forgot this one
minimap2 -ax sr \
    -t $N_THREADS \
    $REF_TRANS/b73_mo17_combined_genome.fa \
    $FASTQ_TRIM_MERGE_DIR/hybrid_b73_mo17_trimmed_bbmerge.fq.gz > $ALIGN_OUT/hybridB73Mo17_to_jointB73Mo17_minimap2.sam

# Align hybrid mo17xb73 to b73
minimap2 -ax sr \
    -t $N_THREADS \
    $REF_TRANS/Zm-B73-REFERENCE-NAM-5.0_canonical_named.fa \
    $FASTQ_TRIM_MERGE_DIR/hybrid_b73_mo17_trimmed_bbmerge.fq.gz > $ALIGN_OUT/hybridB73Mo17_to_b73_minimap2.sam

# Align hybrid mo17xb73 to mo17
minimap2 -ax sr \
    -t $N_THREADS \
    $REF_TRANS/Zm-Mo17-REFERENCE-YAN-1.0_canonical_named.fa \
    $FASTQ_TRIM_MERGE_DIR/hybrid_b73_mo17_trimmed_bbmerge.fq.gz > $ALIGN_OUT/hybridB73Mo17_to_mo17_minimap2.sam


# Get transcript names
cd $REF_TRANS
grep -e ">" Zm-B73-REFERENCE-NAM-5.0_canonical_named.fa > temp.txt
sed 's/>//g' temp.txt > $ALIGN_OUT/all_b73_transcript_ids.txt

grep -e ">" Zm-Mo17-REFERENCE-YAN-1.0_canonical_named.fa > temp.txt
sed 's/>//g' temp.txt > $ALIGN_OUT/all_mo17_transcript_ids.txt

grep -e ">" b73_mo17_combined_genome.fa > temp.txt
sed 's/>//g' temp.txt > $ALIGN_OUT/all_b73_and_mo17_transcript_ids.txt

