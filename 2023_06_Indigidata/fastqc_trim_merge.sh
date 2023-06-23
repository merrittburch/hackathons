#!/usr/bin/env bash

set \
  -o errexit \
  -o nounset \
  -o pipefail

# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-06-13
# Updated... 2023-06-13
#
# Description:
# - Run fast qc, run trimmomatic, and merge paired end reads for Indigidata data
# ------------------------------------------------------------------------------

# Start conda environment
conda activate indigidata

# load fastqc
module load fastqc/0.11.7 

# Change directory to my own folder
mkdir /fs/scratch/PAS2510/Indigidata_2023/merritt
cd /fs/scratch/PAS2510/Indigidata_2023/merritt

# copy all files over
cp -r /fs/scratch/PAS2510/Indigidata_2023/Candida_raw_files /fs/scratch/PAS2510/Indigidata_2023/merritt

# Move to my directory
cd /fs/scratch/PAS2510/Indigidata_2023/merritt

# Analyze just our files
fastqc -t 20 SC5314_CRISPR_MAY1376_R1.fastq.gz -o /fs/scratch/PAS2510/Indigidata_2023/merritt/Candida_raw_files/fastqc_results
fastqc -t 20 SC5314_CRISPR_MAY1376_R2.fastq.gz -o /fs/scratch/PAS2510/Indigidata_2023/merritt/Candida_raw_files/fastqc_results

# Analyze all files 
# RUN_PATH=/fs/scratch/PAS2510/Indigidata_2023/merritt/Candida_raw_files
# cd $RUN_PATH
# mkdir /fs/scratch/PAS2510/Indigidata_2023/merritt/Candida_raw_files/fastqc_results
# for file in $(ls $RUN_PATH)
# do
#     SAMPLE=`basename $file`
#     fastqc -t 20 ${SAMPLE} -o /fs/scratch/PAS2510/Indigidata_2023/merritt/Candida_raw_files/fastqc_results
# done

# Load trimmomatic
module load trimmomatic/0.36

# Run trimmomatic
cd /fs/scratch/PAS2510/Indigidata_2023/merritt/Candida_raw_files
java -jar $TRIMMOMATIC PE -phred33 -threads 28 \
  ./SC5314_CRISPR_MAY1376_R1.fastq.gz \
  ./SC5314_CRISPR_MAY1376_R2.fastq.gz \
  ./trimmed_reads/paired_trimmed_SC5314_CRISPR_MAY1376_R1.fastq.gz \
  ./trimmed_reads/unpaired_trimmed_SC5314_CRISPR_MAY1376_R1.fastq.gz \
  ./trimmed_reads/paired_trimmed_SC5314_CRISPR_MAY1376_R2.fastq.gz \
  ./trimmed_reads/unpaired_trimmed_SC5314_CRISPR_MAY1376_R2.fastq.gz \
  ILLUMINACLIP:/path/to/trim/file.fa:2:30:10 HEADCROP:10 MINLEN:50 SLIDINGWINDOW:4:20 LEADING:5  TRAILING:5 

# Input Read Pairs: 5218953 Both Surviving: 5150401 (98.69%) Forward Only Surviving: 34711 (0.67%) Reverse Only Surviving: 4418 (0.08%) Dropped: 29423 (0.56%)
# TrimmomaticPE: Completed successfully

# Analyze the trimmed files
cd /fs/scratch/PAS2510/Indigidata_2023/merritt/Candida_raw_files
fastqc -t 20 ./trimmed_reads/paired_trimmed_SC5314_CRISPR_MAY1376_R1.fastq.gz \
 -o /fs/scratch/PAS2510/Indigidata_2023/merritt/Candida_raw_files/fastqc_results

fastqc -t 20 ./trimmed_reads/paired_trimmed_SC5314_CRISPR_MAY1376_R2.fastq.gz \
 -o /fs/scratch/PAS2510/Indigidata_2023/merritt/Candida_raw_files/fastqc_results

# Load more packages into conda environment
conda install -c bioconda pandaseq

# Use pandaseq to merge reads
cd /fs/scratch/PAS2510/Indigidata_2023/merritt/Candida_raw_files/trimmed_reads
pandaseq \
  -F \
  -f ./paired_trimmed_SC5314_CRISPR_MAY1376_R1.fastq.gz \
  -r ./paired_trimmed_SC5314_CRISPR_MAY1376_R2.fastq.gz \
  -u unmerged_pandaseq.fa 2> pandastat.txt 1> merged_trimmed_SC5314_CRISPR_MAY1376.fastq




pandaseq -F -f /path/to/your/file/paired_trimmed_SC5314_CRISPR_MAY1376_R1.fastq.gz -r /path/to/your/file/paired_trimmed_SC5314_CRISPR_MAY1376_R2.fastq.gz -u /path/to/your/outputfile/unmerged_pandaseq.fa 2> /path/to/your/outputfile/pandastat.txt 1> /path/to/your/outputfile/merged_trimmed_SC5314_CRISPR_MAY1376.fastq



