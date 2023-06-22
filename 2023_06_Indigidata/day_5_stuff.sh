# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-06-16
# Updated... 2023-06-16
#
# Description:
# - Day 5
# ------------------------------------------------------------------------------

# Change directories
cd /fs/scratch/PAS2510/Indigidata_2023/Candida-genome-stuff

# Load module
module load gatk

# Copy files into my own directory
cd /fs/scratch/PAS2510/Indigidata_2023/merritt/day5
cp /fs/scratch/PAS2510/Indigidata_2023/Candida-genome-stuff/C_albicans_SC5314_A21_current_chromosomes.* ./ 

# Add or replace groups
gatk AddOrReplaceReadGroups \
    -I /fs/scratch/PAS2510/Indigidata_2023/merritt/genome/SC5314_CRISPR_MAY1376_sorted.bam \
    -O /fs/scratch/PAS2510/Indigidata_2023/merritt/day5/SC5314_CRISPR_MAY1376_GATK_sorted.bam \
    --RGLB SC5314_CRISPR_MAY1376 \
    --RGPL ILLUMINA \
    --RGPU unit1 \
    --RGSM SC5314_CRISPR_MAY1376

# Add or replace groups with thier files
gatk AddOrReplaceReadGroups \
    -I /fs/scratch/PAS2510/Indigidata_2023/MAY1717output.sorted.bam \
    -O /fs/scratch/PAS2510/Indigidata_2023/merritt/day5/MAY1717output \
    --RGLB MAY1717output \
    --RGPL ILLUMINA \
    --RGPU unit1 \
    --RGSM MAY1717output


# Mark duplicates
gatk MarkDuplicates \
    -I SC5314_CRISPR_MAY1376_GATK_sorted.bam \
    -O SC5314_CRISPR_MAY1376_markdupl.bam \
    -M SC5314_CRISPR_MAY1376_markdupl.txt \
    -VALIDATION_STRINGENCY SILENT

gatk MarkDuplicates \
    -I MAY1717output.bam \
    -O MAY1717_markdupl.bam \
    -M MAY1717_markdupl.txt \
    -VALIDATION_STRINGENCY SILENT

# Call haplotypes
gatk HaplotypeCaller

 gatk HaplotypeCaller  \
    -R /fs/scratch/PAS2510/Indigidata_2023/merritt/genome/C_albicans_SC5314_A22_current_chromosomes.fasta \
    -I SC5314_CRISPR_MAY1376_markdupl.bam \
    -O SC5314_CRISPR_MAY1376.vcf.gz \
    -ERC GVCF

 gatk HaplotypeCaller  \
    -R /fs/scratch/PAS2510/Indigidata_2023/merritt/genome/C_albicans_SC5314_A22_current_chromosomes.fasta \
    -I MAY1717_markdupl.bam \
    -O SC5314_CRISPR_MAY1376.vcf.gz \
    -ERC GVCF







