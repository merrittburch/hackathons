# ------------------------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2023-06-07
# Updated... 2023-06-07
#
# Description:
# Extract the upstream sequences within 10kb that contain TEs
# (i.e. extract all TE sequences within 10 kb of each gene)
# ------------------------------------------------------------------------------

# Install packages
library(rtracklayer)
library(dplyr)
library(Biostrings)
library(data.table)

# Set directory to b73 gff file
setwd('/Users/mbb262-admin/Library/CloudStorage/Box-Box/git_projects/hackathons/2023_06_teML/data/')

# Import gff file and TE annotation
GFF <- import('Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz')
GFF <- subset(GFF, GFF$type == 'CDS')
GFF <- as.data.frame(GFF)

# Create a function that gets the largest CDS from a group of candidates
maxRange <- function(inputDF){
  res <- range(c(inputDF$start, inputDF$end))
  res <- data.frame(chr = unique(inputDF$seqnames),
                    start = min(res),
                    end = max(res),
                    strand = unique(inputDF$strand))
  res
}

# Use function
GFF$ID <- factor(GFF$ID, levels = unique(GFF$ID))
CDSRange <- by(data = GFF, INDICES = GFF[,"ID"], FUN = maxRange)
CDSRange <- do.call(what = rbind, args = CDSRange)

# extract upstream 5kb using genomic ranges
# (i.e. add on a promoter/5kb to the start of each CDS)
CDSGRange <- GenomicRanges::GRanges(CDSRange)
upstream <- promoters(x = CDSGRange, upstream = 10000, downstream = 0)
upstream_df <- GenomicRanges::as.data.frame(upstream)

# Turn into a datatable object
upstream_df$transcript <- rownames(upstream_df)
upstream_df <- data.table::as.data.table(upstream_df)
data.table::setkey(upstream_df, seqnames, start, end)

# Load in TE gff, format and remove extra stuff
TE_GFF <- import('B73.PLATINUM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.split.gff3.gz')
TE_GFF <- as.data.frame(TE_GFF) %>% filter(type != "knob")
TE_GFF <- as.data.frame(TE_GFF) %>% filter(type != "subtelomere")
TE_GFF <- as.data.frame(TE_GFF) %>% filter(type != "low_complexity")
TE_GFF <- as.data.frame(TE_GFF) %>% filter(type != "centromeric_repeat")
TE_GFF <- as.data.frame(TE_GFF) %>% filter(type != "rDNA_intergenic_spacer_element") %>% 
  select(-"score",-"phase", -"source", -"Sequence_ontology", -"Identity", -"Method", -"motif", -"tsd", -"TSD", -"TIR")

# remove scaffolds 
TE_GFF <- TE_GFF[!grepl("B73_s", TE_GFF$seqnames),]
TE_GFF$seqnames <- gsub("B73_", "", TE_GFF$seqnames)

# turn into datatable object for overlap function
TE_GFF <- data.table::as.data.table(TE_GFF)
data.table::setkey(TE_GFF, seqnames, start, end)

# Find the TEs that overlap with the 10 KB of a gene
# Remove TEs that don't overlap with any genes
together <- data.table::foverlaps(TE_GFF, upstream_df, type = "within") %>% 
  tidyr::drop_na()

# Count the average number of TEs per gene
TEs_per_Gene_10Kb <- together %>% count(transcript)
hist(TEs_per_Gene_10Kb$n)
mean(TEs_per_Gene_10Kb$n)
max(TEs_per_Gene_10Kb$n)

# turn back into a dataframe
together$transcript <- paste(together$transcript, 
                                 together$type,
                                 together$ID,
                                 together$Classification, 
                              sep = "-")
together_df <- data.frame(together) %>% 
  select("seqnames", "i.start", "i.end", "i.width", "i.strand", "transcript")
colnames(together_df) <- c("seqnames", "start", "end", "width", "strand", "transcript")
together_df$seqnames <- as.factor(together_df$seqnames)
rownames(together_df) <- together_df$transcript

# Turn into a genomic range object
gr_together <- GenomicRanges::makeGRangesFromDataFrame(together_df)

# extract sequences
genomePath <- 'Zm-B73-REFERENCE-NAM-5.0.fa.gz'
if(file.exists(paste0(genomePath, ".2bit"))){
  genome <- TwoBitFile(paste0(genomePath, ".2bit"))
}else{
  dna <- readDNAStringSet(genomePath)
  dna <- replaceAmbiguities(dna, new = "N")
  export(dna, paste0(genomePath, ".2bit"))
  rm(dna)
  genome <- TwoBitFile(paste0(genomePath, ".2bit"))
}

# Use either the upstream df or the new together_df to get upstream or TE sequences
chrLength <- seqlengths(genome)
# prs.df <- data.frame(upstream, stringsAsFactors = F)
prs.df <- data.frame(gr_together, stringsAsFactors = F)
prs.df$chrLen <- chrLength[match(prs.df$seqnames, names(chrLength))]

# remove genes that are out of chr range
idx <- which(prs.df$start < 0 | prs.df$end >= prs.df$chrLen)
# if(length(idx) != 0){
#   upstream <- upstream[-c(idx)]
# }
if(length(idx) != 0){
  gr_together <- gr_together[-c(idx)]
}

# Format and write to file
# Seq <- DNAStringSet(getSeq(genome, upstream))
Seq <- DNAStringSet(getSeq(genome, gr_together))
writeXStringSet(x = Seq, filepath = "te_sequence.fasta", format = 'fasta')


