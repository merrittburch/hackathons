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

# What the genome alignment file looks like
# /workdir/panand/arun/Ag-CAM1351-Draft-PanAnd-1.0.fasta.gz
# /workdir/panand/arun/Td-FL_9056069_6-Draft-PanAnd-1.0.fasta.gz
# /workdir/panand/arun/Av-Kellogg1287-8-Draft-PanAnd-1.0.fasta.gz
# /workdir/panand/arun/Td-KS_B6_1-Draft-PanAnd-1.0.fasta.gz

# add these two files together to create `genome_2_align_2.list`
find /workdir/panand/arun -type f > arun.list
find /workdir/panand/aimee/final_contigs -type f > aimee.list
cat *.list > genome_2_align_2.list

# Run minimap in parallel
cd /workdir/hack/aimee_pipeline/sam_alignments
/programs/parallel/bin/parallel -j 5 \
    "/programs/minimap2-2.17/minimap2 -ax splice --eqx -t 20 -I 6G {} /workdir/hack/aimee_pipeline/query_v5_annotatedCDS.fa > {/.}_v5_aligns.sam" :::: /workdir/hack/aimee_pipeline/genome_2_align_2.list


# -----------------------------------
# Run aimee script #3
# CountCDSSamBps.py (convered from ipynb)
# -----------------------------------

# install the kotlin kernel and jupyter lab on cbsu machine
conda create -c conda-forge -c jetbrains -n merrittEnv jupyterlab kotlin-jupyter-kernel

# See what conda environments are present
conda env list

# Activate my conda environment
conda activate merrittEnv

#Start Jupyter notebook server
jupyter lab

# open the website with the link given
# Choose the Kotlin kernel











