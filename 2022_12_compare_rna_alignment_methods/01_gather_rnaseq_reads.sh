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
# Gather B73 reads from Anju's dataset
# ------------------------------------------------------------------------------

# Find individual checks of B73
imeta qu -d sample_title like '2021RNAHybrids_GCGPp1G09' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_GCGPp3B08' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_0716p1D01' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_0716p2H10' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_GCLTp1G09' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_GCLTp3B08' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_GCLBp1A05' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'
imeta qu -d sample_title like '2021RNAHybrids_GCLBp3H10' and instrument_model like '%Illumina NovaSeq 6000%' | grep -v "^-" | awk '{print $2}' | awk '{tmp = $1; getline; print tmp "/" $1}'

# Download individual checks of b73
cd /workdir/tf259/hackathon_dec2022/reads/b73
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R001/MS21R001_CKDL210018333-2a-AK705-GG04_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R001/MS21R001_CKDL210018333-2a-AK705-GG04_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R020/MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R020/MS21R020_CKDL210018333-2a-GG04-AK705_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R041/MS21R041_CKDL210018333-2a-7UDI2717-AK23823_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R041/MS21R041_CKDL210018333-2a-7UDI2717-AK23823_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R046/MS21R046_CKDL210018333-2a-AK30561-AK30560_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R046/MS21R046_CKDL210018333-2a-AK30561-AK30560_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R055/MS21R055_CKDL210018333-2a-AK30571-AK30570_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R055/MS21R055_CKDL210018333-2a-AK30571-AK30570_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R070/MS21R070_CKDL210018333-2a-AK4415-AK5231_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R070/MS21R070_CKDL210018333-2a-AK4415-AK5231_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R076/MS21R076_CKDL210018333-2a-AK30573-AK2486_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R076/MS21R076_CKDL210018333-2a-AK30573-AK2486_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R087/MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R087/MS21R087_CKDL210018333-2a-AK30527-AK30018_HH5V7DSX2_L1_2.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R096/MS21R096_CKDL210018333-2a-GF11-GH12_HH5V7DSX2_L1_1.fq.gz
iget /ibl/home/RawSeqData/RNASeq/Zea/X202SC21060387-Z01-F003/raw_data/MS21R096/MS21R096_CKDL210018333-2a-GF11-GH12_HH5V7DSX2_L1_2.fq.gz





