# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-10-21
# Script to filter 282 mQTLs to mono-genic
# or Mendeliean QTLs
# ---------------------------------------

# Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon")

# Packages
library(ggplot2)

# Load in data
tip_leaf <- read.csv("leaf_tip_gwas_results.csv", header = TRUE)
leaf_base <- read.csv("leaf_bases_gwas_results.csv", header = TRUE)

# Plot all R^2 values
# Tip leaf
ggplot(tip_leaf, aes(x = r2)) + 
  geom_density() +
  ggtitle("Distribution of R^2 values in leaf tip mQTLs")

# All traits and SNPs
ggplot(tip_leaf, aes(Pos, r2)) +
  geom_point() +
  facet_grid(cols = vars(Chr), scales = "free") +
  ggtitle("Leaf tip mQTLs with r^2 > 0.6 by chromosome") +
  theme(axis.text.x = element_text(angle = 90))

# Base leaf
ggplot(leaf_base, aes(x = r2)) + 
  geom_density() +
  ggtitle("Distribution of R^2 values in leaf base mQTLs")


# --------------
# Leaf base
# --------------

# Make a dummy variable
tissue_type = leaf_base

# Subset data
tissue_temp <- tissue_type[which(tissue_type$r2 > 0.6),]
ggplot(tissue_temp, aes(x = r2)) + 
  geom_density() +
  ggtitle("Leaf base mQTLs with r^2 > 0.6")

# Break out by facet
ggplot(tissue_temp, aes(Pos, r2)) +
  geom_point() +
  facet_grid(cols = vars(Chr), scales = "free") +
  ggtitle("Leaf base mQTLs with r^2 > 0.6 by chromosome") +
  theme(axis.text.x = element_text(angle = 90))

# Count the number of features
length(unique(tissue_temp$Trait))

# Count the number of SNPs
length(unique(tissue_temp$Marker))

# --------------
# Leaf tip
# --------------

# Make a dummy variable
tissue_type = tip_leaf

# Subset data
tissue_temp <- tissue_type[which(tissue_type$r2 > 0.6),]
ggplot(tissue_temp, aes(x = r2)) + 
  geom_density() +
  ggtitle("Leaf tip mQTLs with r^2 > 0.6")

# Break out by facet
ggplot(tissue_temp, aes(Pos, r2)) +
  geom_point() +
  facet_grid(cols = vars(Chr), scales = "free") +
  ggtitle("Leaf tip mQTLs with r^2 > 0.6 by chromosome") +
  theme(axis.text.x = element_text(angle = 90))

# Count the number of features
length(unique(tissue_temp$Trait))

# Count the number of SNPs
length(unique(tissue_temp$Marker))

# -------------------
# Look at chr4 peaks
# -------------------

# Investigate chromosome 4 peak in leaf tip data
tissue_type_chr4 <- tissue_temp[which(tissue_temp$Chr == 4 & tissue_temp$Pos > 2.325e08 & tissue_temp$Pos < 2.3437e08),]
ggplot(tissue_type_chr4, aes(Pos, r2)) +
  geom_point() +
  ggtitle("Leaf tip mQTLs with r^2 > 0.6 in chr 4L peak") +
  theme(axis.text.x = element_text(angle = 90))

# Number of unique SNPs here
length(unique(tissue_type_chr4$Marker))

# Number of traits mapping here
length(unique(tissue_type_chr4$Trait))


# Whole dataset - leaf tip
# Investigate chromosome 4 peak in leaf tip data
tissue_type_chr4 <- tip_leaf[which(tip_leaf$Chr == 4 & tip_leaf$Pos > 2.3434e08 & tip_leaf$Pos < 2.3437e08),]
ggplot(tissue_type_chr4, aes(Pos, r2)) +
  geom_point() +
  ggtitle("Leaf tip mQTLs with r^2 > 0.6 in chr 4L peak") +
  theme(axis.text.x = element_text(angle = 90))

# Number of unique SNPs here
length(unique(tissue_type_chr4$Marker))

# Number of traits mapping here
length(unique(tissue_type_chr4$Trait))

# Length of region
234370000-234340000


# Whole dataset - leaf base
# Investigate chromosome 4 peak in leaf tip data
tissue_type_chr4 <- leaf_base[which(leaf_base$Chr == 4 & leaf_base$Pos > 2.3434e08 & leaf_base$Pos < 2.3437e08),]
ggplot(tissue_type_chr4, aes(Pos, r2)) +
  geom_point() +
  ggtitle("Leaf Base mQTLs with r^2 > 0.6 in chr 4L peak") +
  theme(axis.text.x = element_text(angle = 90))

# Number of unique SNPs here
length(unique(tissue_type_chr4$Marker))

# Number of traits mapping here
length(unique(tissue_type_chr4$Trait))



