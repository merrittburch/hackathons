# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-11-13
# Updated... 2022-11-13
#
# Description:
# Compare read counts between different flavors of b73
# ------------------------------------------------------------------------------

# Load packages
library(data.table)
library(dplyr)
library(corrplot)
library(Hmisc)
library(ggrepel)

# Load in metadata
metadata <- read.csv("~/git_projects/hackathons/2022_12_compare_rna_alignment_methods/data/hackathon_expression_metadata.csv")
metadata$id_tissue <- paste0(metadata$short_description, "-", metadata$X.tissue)
metadata <- metadata %>% 
  select("file_name", "id_tissue")

# Load in read counts
tpm <- data.table::fread("~/git_projects/hackathons/2022_12_compare_rna_alignment_methods/data/b73_minimap_count.txt")

# Remove a row with all missing data
tpm <- tpm %>% 
  filter(!V1 == "MS21R001_CKDL210018333-2a-AK705-GG04_HH5V7DSX2_L1_trimmed_bbmerge")

# Remove pipeline names
tpm$V1 <- gsub("_trimmed_bbmerge", "", tpm$V1)

# Merge with metadata names
tpm <- merge(metadata, tpm, by.x = "file_name", by.y = "V1") %>% 
  select(-file_name)

# Look at dimensions
dim(tpm)

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


# Make a correlaiton plot with ggcorplt
corr <- round(cor(trans_tpm), 1)
a <- ggcorrplot::ggcorrplot(corr, hc.order = FALSE, 
                       type = "lower",
                       lab = TRUE,
                       lab_size = 4.5,
                       ggtheme = ggplot2::theme_bw,
                       show.diag = TRUE)
ggsave(plot = a, filename = "~/git_projects/hackathons/2022_12_compare_rna_alignment_methods/images/correlation_plot.png",
        height = 15, width = 25, units = "in")

# Create a PCA
tpm_for_pca <- tpm
rownames(tpm_for_pca) <- tpm$id_tissue
tpm_for_pca <- tpm_for_pca[,-1]
# tpm_for_pca <- tpm_for_pca + 0.01 # add a small pseudocount
pca_b73 <- prcomp(tpm_for_pca, scale = FALSE)
summary(pca_b73)

# Format output
df_pca <- data.frame(pca_b73$x)
df_pca$ids <- rownames(df_pca)
df_pca$tissue <- gsub(".*-", "", df_pca$ids)
df_pca$name <- gsub("-.*", "", df_pca$ids)
df_pca$name <- gsub("_PE75_B73", "", df_pca$name)
df_pca$name <- gsub("PE150_b73_", "", df_pca$name)
df_pca$name <- gsub("PE100_B73_", "", df_pca$name)

# Plot PCA 1 and 2
b <- ggplot(df_pca, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point()+
  xlab("PC1 (34.5%)")+
  ylab("PC2 (25.6%)") +
  ggtitle("PC1 & 2 Ed's Method. Values: NAM, Anju, Simulated") +
  geom_label_repel(aes(label = name), max.overlaps = Inf) 
c <- ggplot(df_pca, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point()+
  xlab("PC1 (34.5%)")+
  ylab("PC2 (25.6%)") +
  ggtitle("PC1 & 2 Ed's Method. Values: NAM, Anju, Simulated") 
ggsave("~/git_projects/hackathons/2022_12_compare_rna_alignment_methods/images/pc1_pc2_label.png",
       b, 
       height = 5, width = 7, units = "in")
ggsave("~/git_projects/hackathons/2022_12_compare_rna_alignment_methods/images/pc1_pc2.png",
       c, 
       height = 5, width = 7, units = "in")


# Plot 2 and 3
d <- ggplot(df_pca, aes(x = PC2, y = PC3, color = tissue)) +
  geom_point() +
  xlab("PC2 (25.6%)")+
  ylab("PC3 (11.6%)")+
  ggtitle("PC1 & 2 Ed's Method. Values: NAM, Anju, Simulated")+
  geom_label_repel(aes(label = name), max.overlaps = Inf) 
e <- ggplot(df_pca, aes(x = PC2, y = PC3, color = tissue)) +
  geom_point() +
  xlab("PC2 (25.6%)")+
  ylab("PC3 (11.6%)")+
  ggtitle("PC1 & 2 Ed's Method. Values: NAM, Anju, Simulated") 
ggsave("~/git_projects/hackathons/2022_12_compare_rna_alignment_methods/images/pc2_pc3_label.png",
       d, 
       height = 5, width = 7, units = "in")
ggsave("~/git_projects/hackathons/2022_12_compare_rna_alignment_methods/images/pc2_pc3.png",
       e, 
       height = 5, width = 7, units = "in")


# ---------------------------------------------
# correlate simulated reads
# ---------------------------------------------

# Correlate simulated reads and plot
sim_sub <- trans_tpm %>% data.frame() %>% 
  select("Simulated_PE100_B73_v01.", "Simulated_PE100_B73_v02.", "Simulated_PE100_B73_v03.")
sim_sub$ids <- rownames(sim_sub)
sim_sub$ids <- gsub("_Zm-B73-REFERENCE-NAM-5.0", "", sim_sub$ids)

# Load in actual count data
obs_sim <- read.delim("~/git_projects/hackathons/2022_12_compare_rna_alignment_methods/data/simulated_counts.txt")

# Merge the two together
together <- merge(sim_sub, obs_sim, by.y = "transcriptID", by.x = "ids")

# Make plots
plot(together$Simulated_PE100_B73_v01., together$sample_01)
# plot(together$Simulated_PE100_B73_v02., together$sample_02)
cor(together$sample_01, together$Simulated_PE100_B73_v01.)

ggplot(together, aes(y = "Simulated_PE100_B73_v01.", x = "sample_01")) +
  geom_point() +
  geom_line()


