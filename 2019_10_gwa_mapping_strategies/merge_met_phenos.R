# ---------------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-10-22
# Script to merge NAM and 282 phenotypes
# for metabolism and related traits
# ---------------------------------------

# Set directory
setwd("~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon")


# ----------
# 282 phenos
# -----------

# Main dataframe (empty right now)
association_blups <- data.frame()

# Function to do merging (remove columns you don't want before running this)
# Assumes you want to merge by pop, entry, and geno code
combine_me_wGeno <- function(dataFrame, nameAppend, byWhat){
  # Add author names
  colnames(dataFrame) <- paste(colnames(dataFrame), nameAppend, sep = "_")
  # Merge with main file
  association_blups <- merge(x=association_blups, y=dataFrame, 
                             by.x = c("Entry_ID_noSpaces", "entry", "Geno_Code"), 
                             by.y = byWhat,
                             all.x = TRUE)
}

# Assumes you want to merge by pop and entry only
combine_me_woGeno <- function(dataFrame, nameAppend, byWhat){
  # Add author names
  colnames(dataFrame) <- paste(colnames(dataFrame), nameAppend, sep = "_")
  # Merge with main file
  association_blups <- merge(x=association_blups, y=dataFrame,
                             by.x = c("Entry_ID_noSpaces", "entry"),
                             by.y = byWhat,
                             all.x = TRUE)
}

# Assumes you want to merge by entry only
combine_me_entry <- function(dataFrame, nameAppend, byWhat){
  # Add author names
  colnames(dataFrame) <- paste(colnames(dataFrame), nameAppend, sep = "_")
  # Merge with main file
  association_blups <- merge(x=association_blups, y=dataFrame,
                             by.x = c("entry"),
                             by.y = byWhat,
                             all.x = TRUE)
}


# -------------------
# Peiffer et al 2014
# -------------------

# Read data
peiffer2014 <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Peiffer2014Genetics_blupPhenos20150325.csv",
                        header = TRUE, na.strings = ".")

# Keep only 282 stuff
peiffer2014 <- peiffer2014[which(peiffer2014$Panel == "ASSO"),]

# # Take mean over duplicated b73 row
# cols <- 6:length(peiffer2014)
# peiffer2014 <- data.table::setDT(peiffer2014)[, lapply(.SD, mean), by=c(names(peiffer2014)[1:5]), .SDcols=cols]

# Add last name
temp <- peiffer2014[,6:16]
colnames(temp) <- paste(colnames(temp), "BLUP_Peiffer2014", sep = "_")

# Make consistent and simple names for panel, pop, entry, and geno code
# Consistent with Buckler 2006
names <- peiffer2014[,1:5]
colnames(names) <- c("Panel", "pop", "Entry_ID", "entry", "Geno_Code")


# Add to main data frame
association_blups <- cbind(names, temp)

# Modify columns to not have dashes, spaces, or underscores in names
association_blups$Entry_ID_noSpaces <- association_blups$Entry_ID
association_blups$Entry_ID_noSpaces <- gsub("-", "", association_blups$Entry_ID_noSpaces)
association_blups$Entry_ID_noSpaces <- gsub(" ", "", association_blups$Entry_ID_noSpaces)
association_blups$Entry_ID_noSpaces <- gsub("_", "", association_blups$Entry_ID_noSpaces)
association_blups$Entry_ID_noSpaces <- gsub("[.]", "", association_blups$Entry_ID_noSpaces)
association_blups$Entry_ID_noSpaces <- tolower(association_blups$Entry_ID_noSpaces)

# Reorganize
association_blups <- association_blups[,c(1:3,17,4:16)]

# Testing by removing some stuff
str(association_blups)
association_blups <- association_blups[,-c(1,2,3,5)]
str(association_blups)
# Take mean over duplicated b73 row
cols <- 3:length(association_blups)
tmp <- data.table::setDT(association_blups)[, lapply(.SD, mean), by=c(names(association_blups)[1:2]), .SDcols=cols]


# Remove columns I don't care about for the hackathon
association_blups <- association_blups[,c(1:3)]

# -----------
# Cook 2012
# -----------

cook2012 <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Cook_etal_2012_kernel_comp_pheno_data-111202/Cook_etal_2012_KernelComposition_BLUPs.csv",
                     header = TRUE, na.strings = ".")
cook2012 <- cook2012[which(cook2012$pop == 27),-c(2,3)]
colnames(cook2012) <- paste(colnames(cook2012), "BLUP_Cook2012", sep = "_")
association_blups <- merge(x=association_blups, y=cook2012, 
                           by.x = "Geno_Code", 
                           by.y = "Geno_Code_BLUP_Cook2012",
                           all.x = TRUE)


# ---------------
# Ionmics data
# ---------------


# ----------------------------
# IBM.AssPan.09PUdataset.csv; 
# ---------------------------

ionomics1 <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Ionomics Hub (iHUB) Pulled 2019-06-14/IBM.AssPan.09PUdataset.csv",
                      header = TRUE, na.strings = ".")

# Subset out 282 based on pop
ionomics1 <- ionomics1[which(ionomics1$pop == 27),-c(1:11, 14, 15,16)]

# Take average over rows with same name (average over 3 taxa replicates)
fun1 <- function(x) aggregate(x~pedigree, ionomics1, mean)
ionomics1 <- data.frame(apply(ionomics1[3:22], 2, fun1))

# Remove pedigree name (my last function wasn't very clean)
ionomics1 <- ionomics1[,-c(seq(3,40,2))]

# Fix column names
names(ionomics1) <- gsub(x = names(ionomics1), pattern = ".x", replacement = "_mean_raw_IBM.AssPan.09PUdataset")

# Merge with main
ionomics1$B11.pedigree <- gsub("_", "", ionomics1$B11.pedigree)
ionomics1$B11.pedigree <- gsub(" ", "", ionomics1$B11.pedigree)
ionomics1$B11.pedigree <- gsub("-", "", ionomics1$B11.pedigree)
ionomics1$B11.pedigree <- gsub("[.]", "", ionomics1$B11.pedigree)
ionomics1$B11.pedigree <- tolower(ionomics1$B11.pedigree)

association_blups <- merge(x=association_blups, y=ionomics1, 
                           by.x = c("Entry_ID_noSpaces"), 
                           by.y = "B11.pedigree",
                           all.x = TRUE)


# --------------------------------
# IBM.AssPan.Ky21.10PUdataset.csv
# --------------------------------

ionomics2 <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Ionomics Hub (iHUB) Pulled 2019-06-14/IBM.AssPan.Ky21.10PUdataset.csv",
                      header = TRUE, na.strings = ".")

# Subset out 282 based on pop
ionomics2 <- ionomics2[which(ionomics2$pop == 27),-c(1:12, 14, 15,16)]

# Take average over rows with same name (average over 3 taxa replicates)
fun1 <- function(x) aggregate(x~znum, ionomics2, mean)
ionomics2 <- data.frame(apply(ionomics2[2:21], 2, fun1))

# Remove pedigree name (my last function wasn't very clean)
ionomics2 <- ionomics2[,-c(seq(3,40,2))]

# Fix column names
names(ionomics2) <- gsub(x = names(ionomics2), pattern = ".x", replacement = "_mean_raw_IBM.AssPan.Ky21.10PUdataset")

# Merge with main
association_blups <- merge(x=association_blups, y=ionomics2, 
                           by.x = "Geno_Code", 
                           by.y = "B11.znum",
                           all.x = TRUE)

# -----------
# Owens 2014
# -----------

owens2014 <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Owens_etal_2014_Genetics/BLUP_24carotenoidTraits.csv",
                      header = TRUE, na.strings = "-")
# Manually merge
colnames(owens2014) <- paste(colnames(owens2014), "BLUP_Owens2014", sep = "_")
owens2014$Maize.Inbred.Line_BLUP_Owens2014 <- gsub("-", "", owens2014$Maize.Inbred.Line_BLUP_Owens2014)
owens2014$Maize.Inbred.Line_BLUP_Owens2014 <- gsub(" ", "", owens2014$Maize.Inbred.Line_BLUP_Owens2014)
owens2014$Maize.Inbred.Line_BLUP_Owens2014 <- gsub("_", "", owens2014$Maize.Inbred.Line_BLUP_Owens2014)
owens2014$Maize.Inbred.Line_BLUP_Owens2014 <- gsub("[.]", "", owens2014$Maize.Inbred.Line_BLUP_Owens2014)
owens2014$Maize.Inbred.Line_BLUP_Owens2014 <- tolower(owens2014$Maize.Inbred.Line_BLUP_Owens2014)

# Merge with main file
association_blups <- merge(x=association_blups, y=owens2014, 
                           by.x = c("Entry_ID_noSpaces"), 
                           by.y = "Maize.Inbred.Line_BLUP_Owens2014",
                           all.x = TRUE)

# ------------
# Lipka 2013
# ------------

lipka2013 <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Lipka_etal_2013_G3/BLUPs_20tocochromanolGrainTraits.csv",
                      header = TRUE, na.strings = ".")
# Manually merge
lipka2013$Sample.ID <- gsub("_", "", lipka2013$Sample.ID)
lipka2013$Sample.ID <- gsub(" ", "", lipka2013$Sample.ID)
lipka2013$Sample.ID <- gsub("-", "", lipka2013$Sample.ID)
lipka2013$Sample.ID <- gsub("[.]", "", lipka2013$Sample.ID)
lipka2013$Sample.ID <- tolower(lipka2013$Sample.ID)
colnames(lipka2013) <- paste(colnames(lipka2013), "BLUP_Lipka2014", sep = "_")
# Merge with main file
association_blups <- merge(x=association_blups, y=lipka2013, 
                           by.x = c("Entry_ID_noSpaces"), 
                           by.y = "Sample.ID_BLUP_Lipka2014",
                           all.x = TRUE)

# -----------
# Harjes 2008
# -----------

harjes2008 <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Harjes_etal_Science_2008/Harjes_etal_Science_2008.csv",
                       header = TRUE, na.strings = ".")
# Format names
harjes2008$Taxa <- gsub("_", "", harjes2008$Taxa)
harjes2008$Taxa <- gsub(" ", "", harjes2008$Taxa)
harjes2008$Taxa <- gsub("-", "", harjes2008$Taxa)
harjes2008$Taxa <- gsub("[.]", "", harjes2008$Taxa)
harjes2008$Taxa <- tolower(harjes2008$Taxa)

colnames(harjes2008) <- paste(colnames(harjes2008), "raw_Harjes2008", sep = "_")
association_blups <- merge(x=association_blups, y=harjes2008, 
                           by.x = c("Entry_ID_noSpaces"), 
                           by.y = "Taxa_raw_Harjes2008",
                           all.x = TRUE)
# 
# # ---------
# # Panzea
# # --------
# 
# panzea <- read.csv("~/Box Sync/Cornell_PhD/phenotype_mapping_gwas/Panzea_traitMatrix/traitMatrix_maize282NAM_v15-130212.csv",
#                    header = TRUE, na.strings = "NaN")
# temp <- panzea[1,] # isolate the environment info row
# temp[] <- paste(col(temp, TRUE), as.matrix(temp), sep = "_") # add environment names after traits
# colnames(panzea) <- temp[1,] # make modified names the column names
# colnames(panzea) <- paste(colnames(panzea), "raw_Panzea", sep = "_") # Add data type
# panzea <- panzea[-1,] # remove environment row
# 
# # Turn data into numeric
# panzea[,2:length(panzea)] <- lapply(panzea[2:length(panzea)], function(x) as.numeric(as.character(x)))
# 
# # Format names
# panzea$`Trait_Header name=env_raw_Panzea` <- gsub("_", "", panzea$`Trait_Header name=env_raw_Panzea`)
# panzea$`Trait_Header name=env_raw_Panzea` <- gsub(" ", "", panzea$`Trait_Header name=env_raw_Panzea`)
# panzea$`Trait_Header name=env_raw_Panzea` <- gsub("-", "", panzea$`Trait_Header name=env_raw_Panzea`)
# panzea$`Trait_Header name=env_raw_Panzea` <- gsub("[.]", "", panzea$`Trait_Header name=env_raw_Panzea`)
# panzea$`Trait_Header name=env_raw_Panzea` <- tolower(panzea$`Trait_Header name=env_raw_Panzea`)
# 
# # Subset data to only protein stuff
# panzea <- panzea[,c(1,183:202)]
# 
# # Merge
# association_blups <- merge(x=association_blups, y=panzea, 
#                            by.x = c("Entry_ID_noSpaces"), 
#                            by.y = "Trait_Header name=env_raw_Panzea",
#                            all.x = TRUE)
# # Remove empty columns (there should be 23) 
# # Hmisc::describe(association_blups) # Check if there is empty data
# association_blups <- Filter(function(x)!all(is.na(x)), association_blups)
# 
# 


# Remove rows with any mossing data
temp <- na.omit(association_blups)

# Add in tassel row
tassel_covar_and_data_names <- c("taxa", rep("data", ncol(association_blups)-1))

# Add in covariate, data, taxa IDs
merged_pcs_phenos <- rbind(tassel_covar_and_data_names, merged_pcs_phenos)


# Write csv
write.csv(association_blups, "association_metabolite_phenos_282_withMissing.txt", quote = F, row.names = F)


# ---------------------
# NAM phenotypes
# ---------------------

# read in nam phenos
all_NAM_phenos <- read.csv("~/git_projects/haplotype_v_snp_gwas/data/all_NAM_phenos.txt", sep="")

# Subset to only the metabolites with lots of data
nam_meta <- all_NAM_phenos[,c(2,9, 36:51,106:108,117:139,461:481)]

# Count missingness
table(colSums(is.na(nam_meta))) # Most traits missing ~500 taxa

# Left with 4k taxa
nam_meta <- na.omit(nam_meta)

# Add in family term
nam_meta$family <- gsub("E[0-9]{4}", "", nam_meta$Geno_Code)

# Rearrange
nam_meta <- nam_meta[,c(1,66,2:65)]

# Add in GBS names
# Load in taxa list with genotypes
NAMFamilies_unionMarkerJoin0.35_main200_taxaNames <- read.delim("~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon/NAMFamilies_unionMarkerJoin0.35_main200_taxaNames.txt")
names <- data.frame(NAMFamilies_unionMarkerJoin0.35_main200_taxaNames$Taxa)
names$short <- gsub(":.*", "", names$NAMFamilies_unionMarkerJoin0.35_main200_taxaNames.Taxa)

# Merge two datasets
nam_meta <- merge(x = names, y = nam_meta, by.x = "short", by.y = "Geno_Code")

# Remove short names
nam_meta <- nam_meta[,-1]

# rename first column
colnames(nam_meta)[1] <- "Taxa"

# Turn into character to allow for rows
nam_meta$Taxa <- as.character(nam_meta$Taxa)
nam_meta$family <- as.character(nam_meta$family)

# Add in tassel headers
tasselNames <- c("taxa", "factor", rep("data", ncol(nam_meta)-2))

# Make colnames the first row
nam_meta <- rbind(tasselNames, colnames(nam_meta), nam_meta)

# Write file
setwd("~/Box Sync/Cornell_PhD/labProjects/hackathons/2019-10-21_hackathon")
write.table(nam_meta, "nam_meta_no_missing.txt", quote = F, row.names = F, col.names = F)


# ---------------------------------
# Look at fast association resilts
# For nam metabolites
# ---------------------------------

# Read in results
fa_results <- read.table("fast_association_nam_allTraits_25pcs.txt", header = T )

# Load the library
library(qqman)

# DTS
dts <- fa_results[which(fa_results$Trait == "Days_To_Silk_BLUP_Sum0607_Buckler2009"),]
dts <- dts[which(dts$p != 0),]
manhattan(dts, chr="Chr", bp="Pos", snp="Marker", p="p",
          genomewideline = -log10(1e-05/4010),
          main = "Days_To_Silk_BLUP_Sum0607_Buckler2009 -> y ~ SNPS + Q + 25PCs")

# Oil results
oil <- fa_results[which(fa_results$Trait == "Oil_Cook2012"),]
oil <- oil[which(oil$p != 0),]
manhattan(oil, chr="Chr", bp="Pos", snp="Marker", p="p", 
          genomewideline = -log10(1e-05/4010),
          main = "Oil_Cook2012 -> y ~ SNPS + Q + 25PCs")

# Fumarate
fumarate <- fa_results[which(fa_results$Trait == "Fumarate_Blup_Wallace2014"),]
fumarate <- fumarate[which(fumarate$p != 0),]
manhattan(fumarate, chr="Chr", bp="Pos", snp="Marker", p="p", 
          genomewideline = -log10(1e-05/4010),
          main = "Fumarate_Blup_Wallace2014 -> y ~ SNPS + Q + 25PCs")

# Nitrogen 


# --------------------------
# Calculate PCs from a GRM
# --------------------------

# Load in kinship/GRM matrix
grm <- read.table("kinship_matrix_nam.txt", header = F)

# Calculate global PCs from kinship matriz
# Centers and decomposes matrix
decomp <- function(G, k=2, center=TRUE) {
  
  # stopifnot(isSymmetric(G))
  
  n <- nrow(G)
  
  # Centering (in case it is not already)
  if (center) {
    H <- diag(n) - matrix(1, nrow=n, ncol=n)/n
    G <- H %*% G %*% H
  }
  
  # Decomposition
  ED <- eigen(G)
  
  # Output PCs
  return(list(PC=ED$vectors[, 1:k] %*% diag(sqrt(ED$values[1:k])),
              variance=ED$values[1:k],
              variance_ratio=ED$values[1:k]/sum(ED$values)))
  
}

# Calculate PCs this way --> centered only
pca_global <- decomp(as.matrix(grm[,-1]), k=25, center = TRUE)

# Get variance
kinship_pcs_var <- data.frame(pca_global$variance_ratio*100)

# Save all pcs
write.table(pca_global$PC, "global_pcs_kinship_nam.txt", quote=F, row.names = F)

# Extraxct out just PCs
pca_global_matrix <- data.frame(pca_global$PC)

# Add back in names
pca_global_matrix <- cbind(grm$V1, pca_global_matrix)

# Fix colname
colnames(pca_global_matrix)[1] <- "taxa"

# Make a character
pca_global_matrix$taxa <- as.character(pca_global_matrix$taxa)

# colnames a row
pca_global_matrix <- rbind(colnames(pca_global_matrix), pca_global_matrix)

# Exrtact out just PCs and add tassel headers
tasselNames <- c("taxa", rep("covariate", ncol(pca_global_matrix)-1))
pca_global_matrix <- rbind(tasselNames, pca_global_matrix)

# Write this file
write.table(pca_global_matrix, "pca_nam_25_global_fromGRM.txt", row.names = F, quote = F, col.names = F)
write.table(kinship_pcs_var, "kinship_pcs25_var.table",row.names = F, quote = F,)


# ------------------------------
# Method #2
# PCA on an adjusted nxp matrix
# for family 
# ------------------------------

# Load packages
library(SNPRelate)
library(gdsfmt)

# Load in Beagle imputed SNPs that were filtered with the 'temp' (aka ames_taxa_gbs_deduplicated.txt, code@bottom)
vcf.fn <- "./nam_snps_main200_maf035.vcf"

# Load in the vcf file
nam_vcf <- snpgdsVCF2GDS(vcf.fn, "nam.gds", method="copy.num.of.ref", 
                            snpfirstdim=TRUE, ignore.chr.prefix = "S")

# Get summary output
snpgdsSummary("nam.gds")

# Check if it was imported correctly
genofile <- snpgdsOpen(nam_vcf)


# -----------------
# Start PCA process
# -----------------

# Function to adjust a genotype matrix with covariates
regress_out <- function(X, Q, include_intercept=TRUE) {
  
  n <- nrow(X)
  
  stopifnot(n == nrow(Q))
  
  # Include intercept to Q
  if (include_intercept & sum(apply(X, 2, var) ==0) == 0) {
    Q <- cbind(1, Q)
  }
  
  # Matrix of projection
  H <- diag(n) - Q %*% solve(crossprod(Q)) %*% t(Q)
  
  # Projecting X out of Q
  H %*% X
}

# Use function
X_list <- snpgdsGetGeno(genofile, snpfirstdim=FALSE, with.id=TRUE)

# Format X matrix
X <- X_list$genotype
rownames(X) <- X_list$sample.id
colnames(X) <- X_list$snp.id

# to get alternate counts instead reference counts
X <- 2-X

# Get family terms
familyTerm <- X_list$sample.id # might need to format
familyTerm <- gsub("E[0-9]{4}.*", "", X_list$sample.id)
familyTerm <- gsub("Z0", "", familyTerm)
familyTerm <- as.matrix(as.numeric(familyTerm))
                        
# Check matrix
str(familyTerm)
str(X)
nrow(familyTerm)
nrow(X)

# Adjust subset of snps using family term
# No longer has taxa ids
adj_snps <- regress_out(X, familyTerm, include_intercept=F) #not working

# Manually try adjustment
Q = familyTerm
H <- diag(nrow(X)) - Q %*% solve(crossprod(Q)) %*% t(Q)

# Projecting X out of Q
adj_snps <- H %*% X

# See if results are different
head(adj_snps)
head(X)

# Calculate and get local PCs
# prcomp is ok with nxp matrix
global_PCs <- prcomp(adj_snps)$x

# Save only a subset of the PCs
subset_local_pcs <- global_PCs[,50]

# Get variance explained by each PC
tmp <- (local_PCs$sdev^2/sum(local_PCs$sdev^2))[1:numLocalPcsOut]


