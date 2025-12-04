# This script analyses and plots the BGC expression data in S. corallicola

# Load packages ----
require(tidyr)
require(dplyr)
require(ggplot2)
require(readr)
require(svglite)
require(stringr)

# Load data ----

gff_raw = readLines("R_input/Sperfeld_V5Y_2_polished.gff3") #S. corallicola feature annotations

BGC = read.delim("R_input/Scoralli_summary_antiSMASH.txt", header = TRUE, sep = "\t", na.strings = c("", "NA"), stringsAsFactors = FALSE) %>% #antiSMASH results (BGCs) 
  mutate(from_antiSMASH = as.integer(gsub(",", "", from_antiSMASH)), #remove commas and convert to integer
         to_antiSMASH = as.integer(gsub(",", "", to_antiSMASH)))

meta = read.delim("R_input/meta.tsv", header = TRUE, sep = "\t", na.strings = c("", "NA"), stringsAsFactors = FALSE) %>% #filter metadata
  filter(condition %in% c("reference_stationary", "coralPOM_stationary"))

counts = read.delim("R_input/21-08-24_htseqcounts_raw.csv", header = TRUE, sep = ",", na.strings = c("", "NA"), stringsAsFactors = FALSE) %>% #raw feature counts
  select(any_of(c("gene_id", meta$sample_BMK)))

core = read.delim("R_input/biosynth-core.list", header = FALSE, sep = "\t", na.strings = c("", "NA"), stringsAsFactors = FALSE)

DGE = read.csv("R_input/2024-08-21_24MS01_coralPOM_stationary_vs_reference_stationary_l0a0.01_results.csv") %>% rename(ID = gene_id) #relevant DESeq2 differential gene expression analysis results

# Clean feature annotations ----

# Identify FASTA start; keep only GFF annotation section
fasta_line = grep("^##FASTA", gff_raw)
if (length(fasta_line) > 0) {
  gff_raw = gff_raw[1:(fasta_line - 1)]
}

# Drop comment lines, including the first 8 header lines
gff_raw = gff_raw[!grepl("^#", gff_raw)]

# Read into a data.frame
ScoralliGFF = read.table(text = gff_raw, sep = "\t", header = FALSE, quote = "", stringsAsFactors = FALSE, comment.char = "")

# Assign column names
colnames(ScoralliGFF) = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
rm(gff_raw, fasta_line)

# Generate S. corallicola feature annotation table ----

ScoralliFeatures = ScoralliGFF %>%
  mutate(ID = str_extract(attributes, "(?<=ID=)[^;]+"), #extract relevant fields from the attributes column
         product_PS = str_extract(attributes, "(?<=product=)[^;]+"),
         gene_PS = str_extract(attributes, "(?<=gene=)[^;]+"),
         KEGG_ko_PS = str_extract(attributes, "(?<=KEGG:)[^,]+"),
         RefSeq_PS = str_extract(attributes, "(?<=RefSeq:)[^,]+"),
         UniParc_PS = str_extract(attributes, "(?<=UniParc:)[^,]+"),
         length = end - start + 1) %>% #calculate gene length
  select(ID, seqid, source, type, start, end, length, strand, RefSeq_PS, UniParc_PS, gene_PS, product_PS, KEGG_ko_PS) #select relevant columns for feature table

# Initialise antiSMASH columns in feature table
ScoralliFeatures[, c("region_antiSMASH", "type_antiSMASH", "similar_cluster_compound_antiSMASH", "similar_cluster_type_antiSMASH", "similarity_antiSMASH")] = NA

# Assign BGC info to genes ----

# For each BGC region, annotate genes that fall within start-to-end
for (i in seq_len(nrow(BGC))) {
  region_idx = which(ScoralliFeatures$start >= BGC$from_antiSMASH[i] &
                     ScoralliFeatures$end <= BGC$to_antiSMASH[i])
  if (length(region_idx) > 0) {
    ScoralliFeatures[region_idx, "region_antiSMASH"] = BGC$region_antiSMASH[i]
    ScoralliFeatures[region_idx, "type_antiSMASH"] = BGC$type_antiSMASH[i]
    ScoralliFeatures[region_idx, "similar_cluster_compound_antiSMASH"] = BGC$similar_cluster_compound_antiSMASH[i]
    ScoralliFeatures[region_idx, "similar_cluster_type_antiSMASH"] = BGC$similar_cluster_type_antiSMASH[i]
    ScoralliFeatures[region_idx, "similarity_antiSMASH"] = BGC$similarity_antiSMASH[i]
  }
}

# Calculate TPM-normalised feature counts and add results to the feature table ----

# Rename columns according to sample mapping
name_mapping = setNames(meta$sample_openBIS, meta$sample_BMK)
counts = counts %>%
  rename_with(~ name_mapping[.x], .cols = -gene_id) %>% #only rename non-gene_id columns
  select(gene_id, sort(setdiff(names(.), "gene_id"))) %>% #sort columns (gene_id first, then alphabetical)
  rename(ID = gene_id) %>% #rename gene_id to ID
  left_join(ScoralliFeatures %>% select(ID, length), by = "ID") #add gene length

# Prepare counts matrix and gene length vector
length_vec = counts$length
counts_mat = counts %>%
  select(-ID, -length) %>%
  as.matrix()
rownames(counts_mat) = counts$ID

# TPM function
TPMfun = function(counts, length) {
  counts_scaled = sweep(counts, 1, length / 1e4, `/`) #divide counts by length in kb
  scaling_factors = colSums(counts_scaled) / 1e6 #scale per million
  sweep(counts_scaled, 2, scaling_factors, `/`) #divide by scaling factor
}

# Apply TPM normalisation
countsTPM = TPMfun(counts_mat, length_vec)
colnames(countsTPM) = paste0(colnames(countsTPM), "_TPM")

# Convert to data frame and join back to feature table
countsTPM_df = countsTPM %>% as.data.frame() %>%
  tibble::rownames_to_column("ID")   # add ID column from rownames

# Explore differential abundance ----

ScoralliFeatures = ScoralliFeatures %>%
  left_join(countsTPM_df, by = "ID") %>%
  left_join(select(DGE, ID, log2FoldChange, padj), by = "ID")

# Plot BGC expression ----

# Filter relevant BGC expression data
BGC_expression <- ScoralliFeatures %>%
  filter(!is.na(region_antiSMASH)) %>% #keep only genes in BGCs
  select(ID, region_antiSMASH, type_antiSMASH, log2FoldChange, padj) %>% 
  filter(ID %in% core$V1) #keep only biosynthetic core genes

# Add statistics and labels
BGC_expression_long <- BGC_expression %>%
  mutate(RegionTypeLabel = paste(region_antiSMASH, type_antiSMASH, sep = " | "), #create a combined label for plotting
         Significant = !is.na(padj) & !is.na(log2FoldChange) & padj < 0.05 & (log2FoldChange < -0.585 | log2FoldChange > 0.585)) %>%
  group_by(RegionTypeLabel) %>%
  mutate(median_Log2FoldChange = median(log2FoldChange, na.rm = TRUE), #calculate median
         count = n()) %>% #calculate number of genes per BGC
  ungroup()

# Plot
ggplot(BGC_expression_long, aes(y = reorder(RegionTypeLabel, median_Log2FoldChange), x = log2FoldChange)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.5) + 
  geom_point(aes(), alpha=0.5, stroke = 0, show.legend = TRUE, size = 2.5) +
  geom_vline(xintercept = 0) +
  geom_text(aes(label= paste0("n=", count), x = -3.3), size = 3, position = position_dodge(width = 1.5), check_overlap = TRUE) +
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

# Save the plot
ggsave(filename = "R_output/BGC_transcription.png", plot = p, width = 8, height = 7)

# Calculate mean gene feature counts per sample ----

# Ensure only numeric columns are used for counting
numeric_counts = counts %>% dplyr::select(where(is.numeric))

# Compute total counts per sample
samples_count = data.frame(sample = colnames(numeric_counts),
                           counts = colSums(numeric_counts, na.rm = TRUE),
                           row.names = NULL)

# Calculate mean and standard deviation
mean_counts = mean(samples_count$counts, na.rm = TRUE)
sd_counts = sd(samples_count$counts, na.rm = TRUE)

# Print the results
cat("Mean of counts:", mean_counts, "\n")
cat("Standard deviation of counts:", sd_counts, "\n")
