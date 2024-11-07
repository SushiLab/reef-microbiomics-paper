###############
## R script for analyzing BGC expression data in S. corallicola
###############

# set working directory
setwd("/gram/sunagawa/OMICS/msperfeld/ELN/03_ACIDO-RNA-SEQ/01_RNA-SEQ/24MSP05_RNA-SEQ/R/")

# load libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(svglite)
library(stringr)

###############
## Open S. corallicola feature annotations (provided by Plasmidsaurus sequencing company)
###############

# Specify the column names
column_names <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Read the file line by line
lines <- readLines("R_input/Sperfeld_V5Y_2_polished.gff3")

# Skip the first 8 lines
lines <- lines[-(1:8)]

# Find the line with "##FASTA"
fasta_start <- grep("##FASTA", lines)

# If such a line exists, remove it and all following lines
if (length(fasta_start) > 0) {
  lines <- lines[1:(fasta_start - 1)]
}

# Write the filtered lines to a temporary file
temp_file <- tempfile()
writeLines(lines, temp_file)

# Read the filtered file using read.table and set the column names
ScoralliGFF <- read.table(temp_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
colnames(ScoralliGFF) <- column_names

# remove temporary files
rm(column_names)
rm(temp_file)
rm(fasta_start)
rm(lines)

###############
## Generate S. corallicola feature annotation table
###############

# Extract the ID value and create a new column
# NOTE: the ID is identical with locus_tag, however, the types "regulatory_region" and "CRISPR" dont have a locus tag
ScoralliGFF <- ScoralliGFF %>% mutate(ID = str_extract(attributes, "(?<=ID=)[^;]+"))

# Extract the product value and create a new column
# NOTE: The values in product and Name are identical
ScoralliGFF <- ScoralliGFF %>% mutate(product = str_extract(attributes, "(?<=product=)[^;]+"))

# Extract the gene value and create a new column
ScoralliGFF <- ScoralliGFF %>% mutate(gene = str_extract(attributes, "(?<=gene=)[^;]+"))

# Extract the KEGG value and create a new column (flagged as PS for Plasmidsaurus)
# NOTE: There are never more than one KEGG cluster ID in the attributes column
ScoralliGFF <- ScoralliGFF %>% mutate(KEGG_PS = str_extract(attributes, "(?<=KEGG:)[^,]+"))

# Extract the RefSeq value and create a new column
ScoralliGFF <- ScoralliGFF %>% mutate(RefSeq = str_extract(attributes, "(?<=RefSeq:)[^,]+"))

# Extract the UniParc value and create a new column
ScoralliGFF <- ScoralliGFF %>% mutate(UniParc = str_extract(attributes, "(?<=UniParc:)[^,]+"))

# Rename some columns
ScoralliGFF <- ScoralliGFF %>% rename(product_PS = product)
ScoralliGFF <- ScoralliGFF %>% rename(gene_PS = gene)
ScoralliGFF <- ScoralliGFF %>% rename(KEGG_ko_PS = KEGG_PS)
ScoralliGFF <- ScoralliGFF %>% rename(RefSeq_PS = RefSeq)
ScoralliGFF <- ScoralliGFF %>% rename(UniParc_PS = UniParc)

# Calculate gene length
ScoralliGFF$length <- ScoralliGFF$end - ScoralliGFF$start + 1

# Select columns of GFF dataframe and make a feature table
ScoralliFeatures <- ScoralliGFF %>% select(ID, seqid, source, type, start, end, length, strand, RefSeq_PS, UniParc_PS, gene_PS, product_PS, KEGG_ko_PS)

###############
## add antiSMASH results (BGCs)
###############

BGC <- read.delim("R_input/Scoralli_summary_antiSMASH.txt", 
                          header = TRUE, 
                          sep = "\t", 
                          na.strings = c("", "NA"), 
                          stringsAsFactors = FALSE)

# Fix problem with comma
BGC$from_antiSMASH <- as.integer(gsub(",", "", BGC$from_antiSMASH))
BGC$to_antiSMASH <- as.integer(gsub(",", "", BGC$to_antiSMASH))

# Add new columns to ScoralliFeatures
ScoralliFeatures <- ScoralliFeatures %>%
  mutate(region_antiSMASH = NA,
         type_antiSMASH = NA,
         similar_cluster_compound_antiSMASH = NA,
         similar_cluster_type_antiSMASH = NA,
         similarity_antiSMASH = NA)

# Iterate over each row of BGC
for (i in 1:nrow(BGC)) {
  from <- BGC$from_antiSMASH[i]
  to <- BGC$to_antiSMASH[i]
  region <- BGC$region_antiSMASH[i]
  type2 <- BGC$type_antiSMASH[i]
  compound <- BGC$similar_cluster_compound_antiSMASH[i]
  cluster_type <- BGC$similar_cluster_type_antiSMASH[i]
  similarity <- BGC$similarity_antiSMASH[i]
  
  # Update ScoralliFeatures based on the condition
  ScoralliFeatures <- ScoralliFeatures %>%
    mutate(region_antiSMASH = ifelse(start >= from & end <= to, region, region_antiSMASH),
           type_antiSMASH = ifelse(start >= from & end <= to, type2, type_antiSMASH),
           similar_cluster_compound_antiSMASH = ifelse(start >= from & end <= to, compound, similar_cluster_compound_antiSMASH),
           similar_cluster_type_antiSMASH = ifelse(start >= from & end <= to, cluster_type, similar_cluster_type_antiSMASH),
           similarity_antiSMASH = ifelse(start >= from & end <= to, similarity, similarity_antiSMASH))
}

###############
## open raw feature counts
###############

counts <- read.delim("R_input/21-08-24_htseqcounts_raw.csv", header = TRUE, sep = ",", na.strings = c("", "NA"), stringsAsFactors = FALSE)

###############
## open meta data table
###############

meta <- read.delim("R_input/meta.tsv", header = TRUE, sep = "\t", na.strings = c("", "NA"), stringsAsFactors = FALSE)

###############
## Keep only samples that were analyzed for publication
###############

# Filter meta to include only the desired conditions
meta <- meta[meta$condition %in% c("reference_stationary", "coralPOM_stationary"), ]

# Get the sample names that meet the condition
matching_samples <- meta$sample_BMK

# Subset counts to keep only the columns that match the sample names
counts <- counts[, colnames(counts) %in% c("gene_id", matching_samples)]

###############
## calculate TPM-normalized feature counts and add results to the feature table
###############

# rename the column names of counts DF
name_mapping <- setNames(meta$sample_openBIS, meta$sample_BMK)
new_colnames <- colnames(counts)
new_colnames[new_colnames != "gene_id"] <- name_mapping[new_colnames[new_colnames != "gene_id"]]
colnames(counts) <- new_colnames

# sort the columns of the counts DF numerically
other_columns <- setdiff(colnames(counts), "gene_id")
sorted_columns <- sort(other_columns)
final_order <- c("gene_id", sorted_columns)
counts <- counts[, final_order]
counts <- counts %>% dplyr::rename(ID = gene_id)

# Add gene length
length_data <- ScoralliFeatures[, c("ID", "length")]
counts <- merge(counts, length_data, by = "ID", all.x = TRUE)

# Make vector with gene length
length <- counts$length

# Prepare for TPM script
rownames(counts) <- counts[, 1]
counts <- counts[, -1]
counts <- counts[, !names(counts) %in% 'length']

# TPM function
TPMfun <- function(counts,length)
{
  counts1 <- sweep(counts,MARGIN=1,(length/10^4),`/`)
  scf <- colSums(counts1)/(10^6)
  return(sweep(counts1,2,scf,`/`))
}

# apply TPM function
countsTPM <- TPMfun(counts,length)
colnames(countsTPM) <- paste0(colnames(countsTPM), "_TPM")

# Join
countsTPM <- cbind(rownames(countsTPM), data.frame(countsTPM, row.names=NULL, check.names = FALSE))
countsTPM <- countsTPM %>% dplyr::rename(ID = "rownames(countsTPM)")
ScoralliFeatures <- ScoralliFeatures %>% left_join(countsTPM, by = "ID")

###############
## Open relevant DESeq2 differential gene expression analysis results and add to the feature table
###############

# Open 
DGE <- read.csv("R_input/2024-08-21_24MS01_coralPOM_stationary_vs_reference_stationary_l0a0.01_results.csv")
DGE <- DGE %>% rename(ID = gene_id)

# Join log2FoldChange and padj columns
ScoralliFeatures <- ScoralliFeatures %>% left_join(DGE %>% select(ID, log2FoldChange, padj), by = "ID")

###############
## Save S. corallicola feature tables with TPM-normalized counts and DESeq2 results
###############

write.table(ScoralliFeatures, file = "R_output/ScoralliFeatures_TPM.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

###############
## Plot BGC expression
###############

# Filter for relevant BGC expression data
BGC_expression <- ScoralliFeatures[!is.na(ScoralliFeatures$region_antiSMASH), ]
BGC_expression <- BGC_expression[, c("ID", "region_antiSMASH", "type_antiSMASH", "log2FoldChange", "padj")]

# Keep only biosynthetic core genes
core <- read.delim("R_input/biosynth-core.list", header = FALSE, sep = "\t", na.strings = c("", "NA"), stringsAsFactors = FALSE)
BGC_expression <- BGC_expression %>% filter(ID %in% core$V1)

# Add statistics and labels
BGC_expression_long <- BGC_expression %>%
  mutate(RegionTypeLabel = paste(region_antiSMASH, type_antiSMASH, sep = " | ")) %>%
  group_by(RegionTypeLabel) %>%
  mutate(median_Log2FoldChange = median(log2FoldChange, na.rm = TRUE),
         count = n(),
         Significant = !is.na(padj) &
           !is.na(log2FoldChange) &
           padj < 0.05 & 
           (log2FoldChange < -0.585 | 
              log2FoldChange > 0.585)) %>%
  ungroup()

# Plot
p <- ggplot(BGC_expression_long, aes(y = reorder(RegionTypeLabel, median_Log2FoldChange), 
                                     x = log2FoldChange)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.5) + 
  geom_point(aes(), alpha=0.5, stroke = 0, show.legend = TRUE, size = 2.5) +
  geom_vline(xintercept = 0) +
  geom_text(aes(label= paste0("n=", count), x = -3.3),
            size = 3, position = position_dodge(width = 1.5), check_overlap = TRUE) +
  coord_cartesian(xlim = c(-3.4, 3.4)) +
  labs(title = "Transcriptional regulation of BGC core genes",
       y = "BGCs",
       x = "Log2 Fold Change",
       color = "Significance") +
  theme_bw() +
  theme(axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        plot.title = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black", linewidth = 0.5),
        legend.position = "none")
p

# Save the plot as an SVG file
ggsave(filename = "R_output/BGC_transcription.svg", plot = p, width = 8, height = 7)
ggsave(filename = "R_output/BGC_transcription.png", plot = p, width = 8, height = 7)

#############
# Calculate mean gene feature counts per sample
###############

#Get total counts per sample
samples_count <- colSums(counts)
samples_count <- data.frame(sample = names(counts), counts = samples_count, row.names = NULL)

# Calculate mean and standard deviation
mean_counts <- mean(samples_count$counts, na.rm = TRUE)
sd_counts <- sd(samples_count$counts, na.rm = TRUE)

# Print the results
cat("Mean of counts:", mean_counts, "\n")
cat("Standard deviation of counts:", sd_counts, "\n")
