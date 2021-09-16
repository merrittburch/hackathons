# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2021-09-10 
#
# Description 
#   - Converting Anju's HARE expression values from v3 gene ids
#   - to v4 gene ids
# ---------------------------------------------------------------

# load packages
library(data.table)
library(dplyr)
library(stringr)

# load in and format translation table v3 --> v4
genes <- data.table::fread("data/maize.v3TOv4.geneIDhistory.txt", 
                           header = FALSE, 
                           drop = c(3:5))
colnames(genes) <- c("v3", "v4")

# remove scaffolds and unnamed genes (have :)
genes <- genes %>% filter(!str_detect(genes$v4, ":")) 


# load in Anju expression dataset and format
hare <- data.table::fread("data/exp_282_max_cis_trans_model.csv")
hare <- hare[,c(18416, 1:18415)] # make taxa names come first

# transpose hare for merging
hare_colnames <- hare$taxa # save taxa names 
hare <- t(hare[,-1]) %>% as.data.frame() # transpose
colnames(hare) <- hare_colnames # change column names in transposed matrix
hare$v3_genes <- rownames(hare) # add rownames


# load in tf dataset
tf <- data.table::fread("data/Zma_TF_list.txt")
tf$tf_name <- gsub(".*_", "", tf$TF_ID) #parse out tf name

# load in Washburn gene family training dataset
wash <- data.table::fread("data/unexpressed_v_expressed.csv",
                          select = c(1:4))


# ------------------------------
# Do gene convsersions
# ------------------------------

# HARE
hare_v4 <- merge(x = hare, y = genes, by.x = "v3_genes", by.y = "v3", all.y = TRUE)
hare_v4 <- hare_v4 %>% select("v4", "B73") # rearragnge
colnames(hare_v4) <- c("v4_gene", "max_cis_rna_value")
hare_v4[is.na(hare_v4)] <- 0
data.table::fwrite(hare_v4, "data/exp_b73_max_cis_trans_model.csv")


# tf binding
tf_v4 <- merge(x = tf, y = genes, by.x = "Gene_ID", by.y =  "v3")
tf_v4$v4_tf_id <- paste0(tf_v4$v4, "_", tf_v4$tf_name) # append on tf names
data.table::fwrite(tf_v4, "data/Zma_TF_list_v3_v4.txt")

# washburn
wash_v4 <- merge(x = wash, y = genes, by.x = "Geneid", by.y = "v3")
wash_v4 <- wash_v4 %>% select("v4", "category", "family_index")
data.table::fwrite(wash_v4, "data/washburn_training_test_v4_unexpressed_v_expressed.csv")


# -------------------------------------------
# Subset out B73 expression from Karl leaf
# -------------------------------------------

# read in file
temp <- data.table::fread("data/L3Base_kremling_formatted_v4_hapmapids.csv")
temp <- data.frame(t(temp[37,]))
temp$gene <- rownames(temp)
colnames(temp) <- c("value", "gene")
temp<- temp[-1,]
temp$value <- as.numeric(as.character(temp$value))
data.table::fwrite(temp, "data/b73_l3base_kremling.csv")



