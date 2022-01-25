# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-01-25 
#
# Description 
# Gather genomes, find sequence within maize v5, align to 
# other refernce genomes and count prevalance
# ---------------------------------------------------------------


# Get genomes
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz


# -----------------------------------
# Run aimee script #1
# ExtractGFFSequences_filterGenes.py
# -----------------------------------

# Order of function export_annotation_string(reference_path, annotation_path, out_path, keep_id_path = None)
# need to create an empty keep_id_path file with correct name (query_v5_annotatedCDS.fa)
# query_genes.csv = v5 gene_transcript names

# modules needed
pip install Bio
pip install gffutils
pip install pyfaidx

# run script
python3 ExtractGFFSequences_filterGenes.py \
    ./Zm-B73-REFERENCE-NAM-5.0.fa \
    ./Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 \
    ./query_genes.csv \
    ./query_v5_annotatedCDS.fa


# -----------------------------------
# Run aimee script #2
# minimapAnnotatedTranscripts.sh
# -----------------------------------


