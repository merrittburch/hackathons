# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-12-15
# Updated... 2022-12-15
#
# Description:
# Normalize hackathon rnaseq reads
# ------------------------------------------------------------------------------

# Load helpful packages
library(DESeq2)
library(data.table)
library(dplyr)
library(Hmisc)
library(ggrepel)


# ------------------------------
# Format Data
# ------------------------------

# Load in metadata
metadata <- read.csv("~/git_projects/hackathons/2022_12_compare_rna_alignment_methods/data/hackathon_expression_metadata.csv")
metadata$id_tissue <- paste0(metadata$short_description, "-", metadata$X.tissue)
metadata <- metadata %>% 
  select("file_name", "id_tissue", "age", "X.tissue")
metadata$file_name <- gsub("_R1.fatq.gz", "", metadata$file_name) # remove fastq file extension from file names

# Remove row with bad sample
metadata <- metadata %>% 
  filter(!file_name == "MS21R001_CKDL210018333-2a-AK705-GG04_HH5V7DSX2_L1")

# Format to deseq2 style
rownames(metadata) <- metadata$file_name
metadata <- metadata %>% select(-file_name)
rownames(metadata) <- gsub("_R1.fastq.gz", "", rownames(metadata)) # remove extra name

# See ordering of metadata file
temp <- rownames(metadata)

# Load in read counts
tpm <- data.table::fread("~/git_projects/hackathons/2022_12_compare_rna_alignment_methods/data/b73_minimap_count.txt")

# Remove a row with all missing data
tpm <- tpm %>% 
  filter(!V1 == "MS21R001_CKDL210018333-2a-AK705-GG04_HH5V7DSX2_L1_trimmed_bbmerge")

# Remove pipeline names
tpm$V1 <- gsub("_trimmed_bbmerge", "", tpm$V1)

# Look at dimensions
dim(tpm)

# Add a pseudocount to the data 
tpm <- cbind(tpm[,1], (tpm[,-1] + 0.01))

# transpose
colnames_tpm <- t(tpm[,1])
trans_tpm <- t(tpm[,-1])
dim(trans_tpm)
head(trans_tpm)
colnames(trans_tpm) <- colnames_tpm[1,]

# look at dimensions again
dim(trans_tpm)
nrow(trans_tpm)
trans_tpm[1:5,1:5]

# Remove .trim from colnames
colnames(trans_tpm) <- gsub(".trim", "", colnames(trans_tpm) )

# Order the columns in the same way as the metadata and check
trans_tpm <- trans_tpm[ ,order(match(colnames(trans_tpm), rownames(metadata))) ]
lala <- cbind(data.frame(rownames(metadata)), data.frame(colnames(trans_tpm)))

# Check the dimnesitons going into deseq2
ncol(trans_tpm)
nrow(metadata)

# Turn to numeric
trans_tpm <- as.numeric(trans_tpm)

# Turn into a deseq object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = trans_tpm,
                              colData = metadata,
                              design = ~ 1)
# Add a pseudocount to data

# Normalize
temp <- DESeq2::DESeq(dds)

