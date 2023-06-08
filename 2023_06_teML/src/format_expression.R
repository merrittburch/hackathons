# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-06-07
# Updated... 2023-06-07
#
# Description:
# Format B73 gene expression data for input into nucleotide transformer
# ------------------------------------------------------------------------------

# Load packages
library(dplyr)

# Load in expression dataset
# wget mbb262@cbsublfs1.biohpc.cornell.edu:/data1/users/mbb262/pleiotropy/interval_data/expression/v4_Transcripts_meanLogExpression.csv
exp_df <- read.csv("/Users/mbb262-admin/Library/CloudStorage/Box-Box/git_projects/hackathons/2023_06_teML/data/v4_Transcripts_meanLogExpression.csv") %>% 
  select("v4_geneIDs", "Mature_Leaf_8")

# Load in xref file, format
xref <- read.delim("/Users/mbb262-admin/Library/CloudStorage/Box-Box/git_projects/hackathons/2023_06_teML/data/B73v4_B73v5_liftoff_genemodel_CDS_xref_shortened.txt")
colnames(xref) <- c("v4", "v5")
xref <- xref[!grepl("chr", xref$v5),]
xref <- xref[!grepl("scaf", xref$v5),]

# Go from v4 to v5
exp_together <- merge(exp_df, xref, by.x = "v4_geneIDs", by.y = "v4") %>% 
  select("v5", "Mature_Leaf_8") %>% 
  distinct(v5, .keep_all = TRUE)

# Load in sequence data
te_seq <- Biostrings::readDNAStringSet("te_sequence.fasta")
te_seq <- data.frame(ID=names(te_seq),sequences=as.character(te_seq))
rownames(te_seq) <- NULL
te_seq$cleanTranscripts <- gsub("_.*", "", te_seq$ID)

# combine expression and sequence data
te_seq_exp <- merge(exp_together, te_seq, by.x = "v5", by.y = "cleanTranscripts") %>% 
  select("ID", "Mature_Leaf_8", "sequences")

# rename columns
colnames(te_seq_exp) <- c("gene", "exp", "teSeq")

# Export as txt file
write.table(te_seq_exp, file = "te_sequence_with_walley_expression", 
            quote = F, row.names = F, sep = "\t")
