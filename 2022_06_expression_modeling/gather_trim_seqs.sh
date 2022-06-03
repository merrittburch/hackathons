# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-06-02 
# Updated... 2022-06-02

# Description 
# Collect gff files of PanAnd lines(10 high quality assemblies from Arun)
# Collect expression data from multiple tissues + some NAM
# trim sequences with trimgalore
# ---------------------------------------------------------------

# -------------------
# Make directories
# -------------------

mkdir /workdir/gffs
mkdir /workdir/expression
mkdir /workdir/expression/raw
mkdir /workdir/expression/trimmed
mkdir /workdir/expression/output_summary


# -------------------
# Gather gff files
# -------------------

# find gff files
ils /ibl/home/assembly_annotations/andropogoneae

# Download gff files (can't figure out how to do this recursively)
cd /workdir/gffs
iget /ibl/home/assembly_annotations/andropogoneae/Andropogon-gerardii_BIND.v1.gff3.gz
iget /ibl/home/assembly_annotations/andropogoneae/Andropogon-virginicus_BIND.v1.gff3.gz
iget /ibl/home/assembly_annotations/andropogoneae/Bothriochloa-laguroides_BIND.v1.gff3.gz
iget /ibl/home/assembly_annotations/andropogoneae/Chrysopogon-serrulatus_BIND.v1.gff3.gz
iget /ibl/home/assembly_annotations/andropogoneae/Cymbopogon-refractus_BIND.v1.gff3.gz
iget /ibl/home/assembly_annotations/andropogoneae/Elionorus-tripsacoides_BIND.v1.gff3.gz
iget /ibl/home/assembly_annotations/andropogoneae/Heteropogon-contortus_BIND.v1.gff3.gz
iget /ibl/home/assembly_annotations/andropogoneae/Poa-pratensis_BIND.v1.gff3.gz
iget /ibl/home/assembly_annotations/andropogoneae/Schizachyrium-scoparium_BIND.v1.gff3.gz
iget /ibl/home/assembly_annotations/andropogoneae/Sorghastrum-nutans_BIND.v1.gff3.gz
iget /ibl/home/assembly_annotations/andropogoneae/Themeda-triandra_BIND.v1.gff3.gz
iget /ibl/home/assembly_annotations/andropogoneae/Tripsacum-dactyloides-northern_BIND.v1.gff3.gz
iget /ibl/home/assembly_annotations/andropogoneae/Tripsacum-dactyloides-southern_BIND.v1.gff3.gz


# -------------------
# Gather RNAseq data
# -------------------

# Full metadata pulled from: https://docs.google.com/spreadsheets/d/1S-8Xb-Q92m6yM7FYCeU4-xHY58JP285FVAXwlaXZ90I/edit#gid=177587468
# Available here: /workdir/expression/AN20_expression_metadata.csv

# RNAseq data, Cinta said samples starting with AN20RN* 
imeta qu -d library_strategy like '%RNA-Seq%' and sample_name like '%AN20RN%' and isPanAnd = 'yes'

# Parse imeta results and then download everything it finds
cd /workdir/expression/raw
parallel -j 15 "cd /workdir/expression/raw; iget -T {}" :::: <(imeta qu -d library_strategy like '%RNA-Seq%' and sample_name like '%AN20RN%' and isPanAnd = 'yes' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}')

# gather a few NAM inbred and hybrid reads
# All growing point tissues
# B73/ky21: MS21R233; B73: MS21R001; Ky21: MS21R021

# B73/ky21: MS21R233
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R233/MS21R233_CKDL210018333-2a-AK17213-7UDI233_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R233/MS21R233_CKDL210018333-2a-AK17213-7UDI233_HH5V7DSX2_L1_2.fq.gz

# B73: MS21R020
# Note B73: MS21R001 --> failed library I think
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R020/MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R020/MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_2.fq.gz

# Ky21: MS21R021
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R021/MS21R021_CKDL210018333-2a-AK2147-7UDI223_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R021/MS21R021_CKDL210018333-2a-AK2147-7UDI223_HH5V7DSX2_L1_2.fq.gz


# ------------------------------
# Trim reads using trim galore
# ------------------------------

N_THREADS=38
FASTQC_RAW_DIR=/workdir/expression/raw
FASTQC_OUT=/workdir/expression/output_summary
FASTQC_TRIM=/workdir/expression/trimmed


## Trim Galore! for paired end reads ----
# would benefit from parallel
for i in $FASTQC_RAW_DIR/*_1.fq.gz
do
    SAMPLE=$(basename ${i} _1.fq.gz)

    echo "Trimming: " ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

    trim_galore \
    --cores $N_THREADS \
    --paired $FASTQC_RAW_DIR/${SAMPLE}_1.fq.gz $FASTQC_RAW_DIR/${SAMPLE}_2.fq.gz \
    --output_dir $FASTQC_TRIM \
    --quality 20 \
    --fastqc \
    --fastqc_args "-o $FASTQC_OUT" \
    --basename $SAMPLE
done

# move summary outputs to a differnt directory
mv /workdir/expression/trimmed/*_trimming_report.txt /workdir/expression/output_summary


# ------------------------------
# Run multiqc on summary outputs
# ------------------------------

export LC_ALL=en_US.UTF-8
export PYTHONPATH=/programs/multiqc-1.10.1/lib64/python3.6/site-packages:/programs/multiqc-1.10.1/lib/python3.6/site-packages
export PATH=/programs/multiqc-1.10.1/bin:$PATH
multiqc /workdir/expression/output_summary/


