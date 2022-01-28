#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-01- 24
#
# Description 
# Download RNAseq expression data for Andropogoeneae and process it
#   - 
# ---------------------------------------------------------------

# -------------------
# Make directories
# -------------------

mkdir /workdir/hack
mkdir /workdir/hack/reads
mkdir /workdir/hack/reads/raw
mkdir /workdir/hack/reads/trimmed
mkdir /workdir/hack/output


# -------------------
# Get data from blfs1
# -------------------

# RNAseq data
# test gathering and parsing paths
# Cinta said: AN20RN* 
# 73 samples total, 5-10 per species
imeta qu -d sample_title like '%AN20RN%' organism like '%Andropogon gerardii%'

# download data
iget /ibl/home/RawSeqData/RNASeq/andropogoneae/Novogene_10_23_2020/raw_data/AN20RN27/AN20RN27_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/andropogoneae/Novogene_10_23_2020/raw_data/AN20RN27/AN20RN27_2.fq.gz

# Search for latest assemblies in key genome
imeta qu -d filetype = 'fasta' and version_status like '%latest%' and isPanAnd = 'yes' and organism like '%Andropogon gerardii%'

# Search + get latest assemblies all 20 genomes
imeta qu -d filetype = 'fasta' and version_status like '%latest%' and isPanAnd = 'yes'
parallel -j 15 "cd /workdir/panand/arun/; iget -T {}" :::: <(imeta qu -d filetype = 'fasta' and version_status like '%latest%' and isPanAnd = 'yes' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}')

# Gather all of Aimee's genomes
scp -r mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/ajs692/panand_data/megahit_shortread_assemblies/final_contigs/ /workdir/panand/aimee


# Selecting different one as travis for RNAseq alignment
# collection: /ibl/home/assemblies/andropogoneae/private/Corteva_pop
# dataObj: Ag-CAM1351-Draft-PanAnd-1.0
# Andropogon gerardii
iget /ibl/home/assemblies/andropogoneae/private/Corteva_pop/Ag-CAM1351-Draft-PanAnd-1.0.fasta.gz


# ------------------------------
# Trim reads using trim galore
# ------------------------------

trim_galore \
  --cores 8 \
  --paired /workdir/hack/AN20RN27_1.fq.gz /workdir/hack/AN20RN27_2.fq.gz \
  --output_dir /workdir/hack \
  --quality 20 \
  --fastqc \
  --fastqc_args "-o /workdir/hack/" \
  --basename ag


# ------------------------------
# Map reads using hisat2
# ------------------------------

# Gather paths
source /programs/HISAT2/hisat2.sh

# index genome
hisat2-build \
  -p 20 \
  /workdir/hack/Ag-CAM1351-Draft-PanAnd-1.0.fasta \
  /workdir/hack

# map reads
hisat2 \
  -p 30 \
  -x /workdir/hack/hisat2_index/hack \
  -1 /workdir/hack/ag_val_1.fq.gz \
  -2 /workdir/hack/ag_val_2.fq.gz \
  -S /workdir/hack/AN20RN27_output.sam --new-summary \
  2> /workdir/hack/AN20RN27_summary.txt


# ------------------------------
# Map reads using minimap2
# ------------------------------

# Minimap splice aware
/programs/minimap2-2.17/minimap2 \
  -t 15 -ax splice Ag-CAM1351-Draft-PanAnd-1.0.fasta \
  ag_val_1.fq.gz ag_val_1.fq.gz > ag_aln.sam 


# Minimap short read, not splice aware
/programs/minimap2-2.17/minimap2 \
  -t 15 -ax sr Ag-CAM1351-Draft-PanAnd-1.0.fasta \
  ag_val_1.fq.gz ag_val_1.fq.gz > ag_aln_sr.sam 

