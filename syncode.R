### Install R packages if not already installed ###


install.packages("tidyverse")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install("tximport")
BiocManager::install("DESeq2")

library(DESeq2)
library(tidyverse)

# Load count matrix
countData <- read.csv("C:/Users/epige/OneDrive/Documents/Synovial-rnaseq/Syn_ gene_count.csv", header=T, row.names=1)
# Load required libraries
library(DESeq2)
library(tidyverse)

# Load count matrix
countData <- read.csv("C:/Users/epige/OneDrive/Documents/Synovial-rnaseq/Syn_ gene_count.csv", header=T, row.names=1)

# Create coldata directly in R
coldata <- data.frame(
  SampleID = c("F83", "F85", "F87", "F101", "F105", "F109", "F121", "F125", "F127", "F145", "F143", "F147"),
  Condition = c("HU_No_Ex_WT_No_recovery", "HU_No_Ex_WT_No_recovery", "HU_No_Ex_WT_No_recovery",
                "NL_No_Ex_WT_No_recovery", "NL_No_Ex_WT_No_recovery", "NL_No_Ex_WT_No_recovery",
                "HU_Ex_WT_No_recovery", "HU_Ex_WT_No_recovery", "HU_Ex_WT_No_recovery",
                "NL_Ex_WT_No_recovery", "NL_Ex_WT_No_recovery", "NL_Ex_WT_No_recovery"),
  stringsAsFactors = FALSE
)
rownames(coldata) <- coldata$SampleID

# Ensure the order of samples in countData matches coldata
countData <- countData[, rownames(coldata)]

# Check if there are any non-integer values in countData
if(any(countData %% 1 != 0)) {
  # Round to nearest integer if there are decimal values
  countData <- round(countData)
}

# Convert to matrix
countData <- as.matrix(countData)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = coldata,
                              design = ~ Condition)

# Remove lowly expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds)

# Function to get results and format them
get_results <- function(dds, contrast, name) {
  res <- results(dds, contrast=contrast, alpha=0.05)
  res_df <- as.data.frame(res) %>%
    dplyr::select(log2FoldChange, pvalue, padj) %>%
    dplyr::rename_with(~paste0(name, "_", .))
  return(res_df)
}

# List of all pairwise comparisons
comparisons <- list(
  c("Condition", "HU_No_Ex_WT_No_recovery", "NL_No_Ex_WT_No_recovery"),
  c("Condition", "HU_Ex_WT_No_recovery", "NL_Ex_WT_No_recovery"),
  c("Condition", "HU_No_Ex_WT_No_recovery", "HU_Ex_WT_No_recovery"),
  c("Condition", "NL_No_Ex_WT_No_recovery", "NL_Ex_WT_No_recovery")
)

comparison_names <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

# Get results for all comparisons
all_results <- purrr::map2(comparisons, comparison_names, ~get_results(dds, .x, .y)) %>%
  purrr::reduce(cbind)

# Add gene names as the first column
all_results <- all_results %>%
  dplyr::mutate(gene = rownames(.), .before = 1)

# Write results to CSV
write.csv(all_results, "DESeq2_all_comparisons_results_wide.csv", row.names = FALSE)

# Print summary of results
cat("Summary of differentially expressed genes:\n")
summary_results <- map2_dfr(comparison_names, comparison_names, function(name, full_name) {
  data.frame(
    comparison = full_name,
    total_DE = sum(all_results[[paste0(name, "_padj")]] < 0.05, na.rm = TRUE),
    up_regulated = sum(all_results[[paste0(name, "_padj")]] < 0.05 & all_results[[paste0(name, "_log2FoldChange")]] > 0, na.rm = TRUE),
    down_regulated = sum(all_results[[paste0(name, "_padj")]] < 0.05 & all_results[[paste0(name, "_log2FoldChange")]] < 0, na.rm = TRUE)
  )
})

print(summary_results, n = Inf)



# Save summary results to CSV
write.csv(summary_results, "DESeq2_summary_results.csv", row.names = FALSE)

# Add gene IDs as a column for mapping
all_results$Gene_id <- rownames(all_results)

# Use biomaRt to get gene names from Ensembl IDs (adjust species as needed)
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl") # Change dataset for other species

gene_info <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                   filters='ensembl_gene_id',
                   values=all_results$Gene_id,
                   mart=ensembl)

# Merge gene names into results
all_results <- left_join(all_results, gene_info, by=c("Gene_id" = "ensembl_gene_id"))

# Reorder columns to have gene name first (if desired)
all_results <- all_results %>%
  dplyr::select(external_gene_name, everything())

# Write all results to CSV with gene names included
write.csv(all_results, "DESeq2_all_comparisons_results_with_gene_names.csv", row.names = FALSE)

# Print summary of results
cat("Summary of differentially expressed genes:\n")
summary_results <- purrr::map2_dfr(comparison_names, comparison_names, function(name, full_name) {
  data.frame(
    comparison = full_name,
    total_DE = sum(all_results[[paste0(name, "_padj")]] < 0.05, na.rm = TRUE),
    up_regulated = sum(all_results[[paste0(name, "_padj")]] < 0.05 & all_results[[paste0(name, "_log2FoldChange")]] > 0, na.rm = TRUE),
    down_regulated = sum(all_results[[paste0(name, "_padj")]] < 0.05 & all_results[[paste0(name, "_log2FoldChange")]] < 0, na.rm = TRUE)
  )
})

# Save summary results to CSV
write.csv(summary_results, "DESeq2_summary_results.csv", row.names = FALSE)

print(summary_results, n = Inf)

# pca plot modified 
rld <- rlog(dds, blind = TRUE)
pca_plot <- plotPCA(rld, intgroup="Condition", ntop=500)+ theme_minimal(base_size=18) +
  labs(title = "PCA plot for synovial tissue", x="PC1", Y="PC2")

pdf(file = "pca_plot_synovial.pdf", width = 10, height = 8)
print(pca_plot)
dev.off()


#  HEATMAP

# Load required libraries
library(pheatmap)
library(tidyverse)
library(RColorBrewer)

# Read in your data (replace 'your_data.csv' with your actual file name)

duplicate_genes <- all_results$external_gene_name[duplicated(all_results$external_gene_name)]
print(paste("Number of duplicate gene names:", length(duplicate_genes)))
print("Duplicate gene names:")
print(duplicate_genes)

library(dplyr)

all_results <- all_results %>%
  group_by(external_gene_name) %>%
  mutate(unique_name = if(n() > 1) paste0(external_gene_name, "_", row_number()) else external_gene_name) %>%
  ungroup()

# Now use unique_name for row names
heatmap_data <- all_results %>%
  dplyr::select(unique_name, ends_with("log2FoldChange")) %>%
  column_to_rownames("unique_name")

# Rename columns to make them shorter
colnames(heatmap_data) <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

# Step 3: Select top 50 most variable genes
top_50_genes <- heatmap_data %>%
  mutate(variance = apply(., 1, var)) %>%
  arrange(desc(variance)) %>%
  head(50) %>%
  dplyr::select(-variance)

# Step 4: Prepare annotation data
# Create a dataframe for column annotation (conditions)
column_annotation <- data.frame(
  Comparison = factor(colnames(top_50_genes)),
  row.names = colnames(top_50_genes)
)

# Create a color palette for the comparisons
comparison_colors <- brewer.pal(4, "Set2")
names(comparison_colors) <- levels(column_annotation$Comparison)

# Step 5: Prepare heatmap color palette
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

# Step 6: Create the heatmap
heatmap_plot <- pheatmap(top_50_genes,
                         color = heatmap_colors,
                         cluster_rows = TRUE,
                         cluster_cols = FALSE,
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         annotation_col = column_annotation,
                         annotation_colors = list(Comparison = comparison_colors),
                         main = "Top 50 Most Variable Genes - Log2 Fold Changes",
                         fontsize_row = 8,
                         fontsize_col = 10,
                         angle_col = 45,
                         cellwidth = 20,
                         cellheight = 12,
                         border_color = NA,
                         filename = "heatmap_top50_variable_genes.png",
                         width = 12,
                         height = 14)

# Step 7: Display the plot in R (optional)
print(heatmap_plot)

# Step 8: Save gene list (optional)
write.csv(rownames(top_50_genes), "top_50_variable_genes.csv", row.names = FALSE)

# Step 9: Print summary statistics
cat("\nSummary Statistics:\n")
cat("Number of genes in heatmap:", nrow(top_50_genes), "\n")
cat("Range of log2 fold changes:", round(range(top_50_genes), 2), "\n")
cat("Median log2 fold change:", round(median(unlist(top_50_genes)), 2), "\n")


# volcano plots
# Function to create volcano plot with white background
create_volcano_plot <- function(data, x, y, title) {
  ggplot(data, aes(x = !!sym(x), y = -log10(!!sym(y)))) +
    geom_point(aes(color = abs(!!sym(x)) > 1 & !!sym(y) < 0.05), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.5) +
    geom_text_repel(
      data = . %>% filter(abs(!!sym(x)) > 1 & !!sym(y) < 0.05) %>% top_n(20, wt = abs(!!sym(x))),
      aes(label = external_gene_name),
      max.overlaps = 20,
      box.padding = 0.5,
      segment.color = "grey50",
      size = 3
    ) +
    theme_bw() +  # White background with a border
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-Log10 P-value")
}

# Create volcano plots
volcano_plots <- list(
  HU_No_Ex_vs_NL_No_Ex = create_volcano_plot(all_results, "HU_No_Ex_vs_NL_No_Ex_log2FoldChange", "HU_No_Ex_vs_NL_No_Ex_pvalue", "HU No Ex vs NL No Ex"),
  HU_Ex_vs_NL_Ex = create_volcano_plot(all_results, "HU_Ex_vs_NL_Ex_log2FoldChange", "HU_Ex_vs_NL_Ex_pvalue", "HU Ex vs NL Ex"),
  HU_No_Ex_vs_HU_Ex = create_volcano_plot(all_results, "HU_No_Ex_vs_HU_Ex_log2FoldChange", "HU_No_Ex_vs_HU_Ex_pvalue", "HU No Ex vs HU Ex"),
  NL_No_Ex_vs_NL_Ex = create_volcano_plot(all_results, "NL_No_Ex_vs_NL_Ex_log2FoldChange", "NL_No_Ex_vs_NL_Ex_pvalue", "NL No Ex vs NL Ex")
)

# Save volcano plots
for (name in names(volcano_plots)) {
  ggsave(paste0("volcano_plot_", name, ".png"), volcano_plots[[name]], width = 10, height = 8, dpi = 300)
}


# Step 8: Print summary statistics for each comparison
for (comparison in c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")) {
  cat("\nSummary Statistics for", comparison, ":\n")
  
  total_de <- sum(all_results[[paste0(comparison, "_padj")]] < 0.05 & abs(all_results[[paste0(comparison, "_log2FoldChange")]]) > 1, na.rm = TRUE)
  up_regulated <- sum(all_results[[paste0(comparison, "_padj")]] < 0.05 & all_results[[paste0(comparison, "_log2FoldChange")]] > 1, na.rm = TRUE)
  down_regulated <- sum(all_results[[paste0(comparison, "_padj")]] < 0.05 & all_results[[paste0(comparison, "_log2FoldChange")]] < -1, na.rm = TRUE)
  
  cat("Total DE genes:", total_de, "\n")
  cat("Up-regulated:", up_regulated, "\n")
  cat("Down-regulated:", down_regulated, "\n")
  
  top_genes <- all_results %>%
    filter(!!sym(paste0(comparison, "_padj")) < 0.05) %>%
    arrange(desc(abs(!!sym(paste0(comparison, "_log2FoldChange"))))) %>%
    head(10)
  
  cat("Top 10 DE genes:\n")
  print(top_genes$external_gene_name)
}


##############
# GSEA ANALYSIS

# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # For mouse genes. Use org.Hs.eg.db for human genes
library(dplyr)
library(ggplot2)

# Load count matrix
all_results <- read.csv("C:/Users/epige/OneDrive/Documents/Synovial-rnaseq/DESeq2_all_comparisons_results_with_gene_names.csv", header=T, row.names=1)
# Load required libraries
all_results <- read.csv("C:/Users/epige/OneDrive/Documents/Synovial-rnaseq/DESeq2_all_comparisons_results_with_gene_names.csv", header=TRUE)

library(dplyr)

all_results <- all_results %>%
  mutate(row_name = if_else(is.na(.[,1]) | .[,1] == "", paste0("Unknown_", row_number()), .[,1])) %>%
  distinct(row_name, .keep_all = TRUE)

rownames(all_results) <- all_results$row_name
all_results$row_name <- NULL

# Check the structure of your data
str(all_results)

# Assuming 'gene' is the column with Entrez IDs and 'HU_No_Ex_vs_NL_No_Ex_log2FoldChange' is one of your comparisons
gene_list <- all_results$HU_No_Ex_vs_NL_No_Ex_log2FoldChange
names(gene_list) <- all_results$gene
gene_list <- sort(gene_list, decreasing = TRUE)

# Load necessary libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # For mouse genes. Use org.Hs.eg.db for human genes

library(biomaRt)
library(dplyr)

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

# Get the mapping between Ensembl IDs and Entrez IDs
gene_ids <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id'),
                  filters = 'ensembl_gene_id',
                  values = all_results$gene,
                  mart = ensembl)

# Merge this information with our results
all_results <- all_results %>%
  left_join(gene_ids, by = c("gene" = "ensembl_gene_id"))

# Remove genes without Entrez IDs
all_results <- all_results %>% filter(!is.na(entrezgene_id))


# Function to perform GSEA for a single comparison
perform_gsea <- function(data, comparison) {
  gene_list <- data[[paste0(comparison, "_log2FoldChange")]]
  names(gene_list) <- data$entrezgene_id
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea_result <- gseGO(geneList = gene_list,
                       ont = "BP",
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       verbose = TRUE)
  
  return(gsea_result)
}

# List of comparisons
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

# Perform GSEA for each comparison
gsea_results <- lapply(comparisons, function(comp) perform_gsea(all_results, comp))
names(gsea_results) <- comparisons

# Print summary of GSEA results
for (comp in names(gsea_results)) {
  cat("\nGSEA Results for", comp, ":\n")
  print(head(gsea_results[[comp]]))
}

# Save GSEA results
for (comp in names(gsea_results)) {
  write.csv(as.data.frame(gsea_results[[comp]]), file = paste0("GSEA_results_", comp, ".csv"))
}



##########dot plots 
# Load required libraries
library(biomaRt)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(UpSetR)

# Read in the data
all_results <- read.csv("C:/Users/epige/OneDrive/Documents/Synovial-rnaseq/DESeq2_all_comparisons_results_with_gene_names.csv", header=TRUE)

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

# Get the mapping between Ensembl IDs and Entrez IDs
gene_ids <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id'),
                  filters = 'ensembl_gene_id',
                  values = all_results$gene,
                  mart = ensembl)

# Merge this information with our results
all_results <- all_results %>%
  left_join(gene_ids, by = c("gene" = "ensembl_gene_id"), relationship = "many-to-many")

# Remove genes without Entrez IDs and clean up the data frame
all_results <- all_results %>% 
  filter(!is.na(entrezgene_id.y)) %>%
  select(external_gene_name, gene, ends_with("log2FoldChange"), entrezgene_id = entrezgene_id.y)

# Function to perform GSEA for a single comparison
perform_gsea <- function(data, comparison) {
  gene_list <- data[[paste0(comparison, "_log2FoldChange")]]
  names(gene_list) <- data$entrezgene_id
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea_result <- gseGO(geneList = gene_list,
                       ont = "BP",
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       verbose = TRUE)
  
  return(gsea_result)
}

# List of comparisons
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

# Perform GSEA for each comparison
gsea_results <- lapply(comparisons, function(comp) perform_gsea(all_results, comp))
names(gsea_results) <- comparisons

# Function to create and save dot plot with white background
create_dotplot <- function(gsea_result, comparison) {
  dotplot <- dotplot(gsea_result, showCategory=20, title=paste("GSEA -", comparison)) +
    theme_bw() +  # White background
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size=12),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(paste0("GSEA_dotplot_", comparison, ".png"), dotplot, width=12, height=10, dpi=300)
  
  return(dotplot)
}

# Create dot plots for all comparisons
dot_plots <- lapply(names(gsea_results), function(comp) create_dotplot(gsea_results[[comp]], comp))

# Print summary of GSEA results and save to CSV
for (comp in names(gsea_results)) {
  cat("\nGSEA Results for", comp, ":\n")
  print(head(gsea_results[[comp]]))
  
  write.csv(as.data.frame(gsea_results[[comp]]), file = paste0("GSEA_results_", comp, ".csv"), row.names = FALSE)
}

# Create combined plot of top pathways for all comparisons with white background
create_combined_dotplot <- function(gsea_results_list, top_n = 5) {
  top_pathways <- lapply(names(gsea_results_list), function(comp) {
    data.frame(
      Comparison = comp,
      Description = head(gsea_results_list[[comp]]@result$Description, top_n),
      NES = head(gsea_results_list[[comp]]@result$NES, top_n),
      p.adjust = head(gsea_results_list[[comp]]@result$p.adjust, top_n)
    )
  })
  
  combined_data <- do.call(rbind, top_pathways)
  
  combined_plot <- ggplot(combined_data, aes(x = Comparison, y = Description)) +
    geom_point(aes(size = abs(NES), color = NES)) +
    scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    theme_bw() +  # White background
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = "Top Enriched Pathways Across Comparisons",
         size = "Absolute NES",
         color = "NES")
  
  ggsave("Combined_top_pathways_dotplot.png", combined_plot, width = 14, height = 12, dpi = 300)
  
  return(combined_plot)
}

# Create and save the combined plot
combined_plot <- create_combined_dotplot(gsea_results)
print(combined_plot)

################### upset plot

library(ComplexHeatmap)


# Function to create UpSet plot using UpSetR
create_upset_plot <- function(gsea_results_list, top_n = 20) {
  # Extract top pathways for each comparison
  top_pathways <- lapply(gsea_results_list, function(x) {
    head(x@result$ID, top_n)
  })
  
  # Create a binary matrix for UpSet plot
  pathway_matrix <- fromList(top_pathways)
  
  # Create the UpSet plot
  upset(pathway_matrix,
        sets.bar.color = "darkgrey",
        main.bar.color = "black",
        matrix.color = "darkred",
        point.size = 3,
        line.size = 1,
        order.by = "freq")
}

# Create and save the UpSet plot
png("GSEA_UpSet_plot.png", width = 12, height = 8, units = "in", res = 300)
create_upset_plot(gsea_results)
dev.off()

############################ kegg analysis

# 1. Prepare the gene list
# We'll use log2FoldChange for ranking and adjusted p-value for filtering
gene_list <- all_results$HU_No_Ex_vs_NL_No_Ex_log2FoldChange
names(gene_list) <- all_results$entrezgene_id
gene_list <- sort(gene_list, decreasing = TRUE)

# 2. Perform KEGG pathway analysis
kegg_result <- enrichKEGG(gene = names(gene_list),
                          organism = 'mmu',  # 'mmu' for mouse, 'hsa' for human
                          keyType = 'ncbi-geneid',
                          pvalueCutoff = 0.05,
                          pAdjustMethod = 'BH',
                          minGSSize = 10,
                          maxGSSize = 500,
                          use_internal_data = FALSE)

# 3. View results
head(kegg_result)

# 4. Save results to a CSV file
write.csv(as.data.frame(kegg_result), "KEGG_enrichment_results.csv", row.names = FALSE)

# 5. Visualize results

# Dotplot
dotplot <- dotplot(kegg_result, showCategory = 20) +
  ggtitle("Top 20 Enriched KEGG Pathways") +
  theme(axis.text.y = element_text(size = 8))

ggsave("KEGG_dotplot.png", dotplot, width = 12, height = 10, dpi = 300)

# Barplot
barplot <- barplot(kegg_result, showCategory = 20) +
  ggtitle("Top 20 Enriched KEGG Pathways") +
  theme(axis.text.y = element_text(size = 8))

ggsave("KEGG_barplot.png", barplot, width = 12, height = 10, dpi = 300)

# 6. Perform GSEA with KEGG pathways
kegg_gsea <- gseKEGG(geneList = gene_list,
                     organism = 'mmu',
                     keyType = 'ncbi-geneid',
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = 'BH',
                     verbose = FALSE)

# 7. Visualize GSEA results
gsea_plot <- gseaplot2(kegg_gsea, geneSetID = 1:5, pvalue_table = TRUE)
ggsave("KEGG_GSEA_plot.png", gsea_plot, width = 12, height = 10, dpi = 300)

# 8. Save GSEA results
write.csv(as.data.frame(kegg_gsea), "KEGG_GSEA_HU_No_Ex_vs_NL_No_Ex_results.csv", row.names = FALSE)

########################################################################
# Assuming 'all_results' is your dataframe with DESeq2 results

# Define the comparisons
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

# Function to perform KEGG analysis for a single comparison
perform_kegg_analysis <- function(data, comparison) {
  # Prepare the gene list
  gene_list <- data[[paste0(comparison, "_log2FoldChange")]]
  names(gene_list) <- data$entrezgene_id
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Perform KEGG pathway analysis
  kegg_result <- enrichKEGG(gene = names(gene_list),
                            organism = 'mmu',  # 'mmu' for mouse, 'hsa' for human
                            keyType = 'ncbi-geneid',
                            pvalueCutoff = 0.05,
                            pAdjustMethod = 'BH',
                            minGSSize = 10,
                            maxGSSize = 500,
                            use_internal_data = FALSE)
  
  # Save results to a CSV file
  write.csv(as.data.frame(kegg_result), paste0("KEGG_enrichment_results_", comparison, ".csv"), row.names = FALSE)
  
  # Create and save dotplot
  dotplot <- dotplot(kegg_result, showCategory = 20) +
    ggtitle(paste("Top 20 Enriched KEGG Pathways -", comparison)) +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(paste0("KEGG_dotplot_", comparison, ".png"), dotplot, width = 12, height = 10, dpi = 300)
  
  # Perform GSEA with KEGG pathways
  kegg_gsea <- gseKEGG(geneList = gene_list,
                       organism = 'mmu',
                       keyType = 'ncbi-geneid',
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = 'BH',
                       verbose = FALSE)
  
  # Save GSEA results
  write.csv(as.data.frame(kegg_gsea), paste0("KEGG_GSEA_results_", comparison, ".csv"), row.names = FALSE)
  
  # Create and save GSEA plot
  gsea_plot <- gseaplot2(kegg_gsea, geneSetID = 1:5, pvalue_table = TRUE)
  ggsave(paste0("KEGG_GSEA_plot_", comparison, ".png"), gsea_plot, width = 12, height = 10, dpi = 300)
  
  return(list(enrichment = kegg_result, gsea = kegg_gsea))
}

# Perform KEGG analysis for all comparisons
kegg_results <- map(comparisons, ~perform_kegg_analysis(all_results, .x))
names(kegg_results) <- comparisons

# Create a combined dotplot for all comparisons
create_combined_dotplot <- function(kegg_results_list, top_n = 5) {
  top_pathways <- map_dfr(names(kegg_results_list), function(comp) {
    result <- kegg_results_list[[comp]]$enrichment
    if (nrow(result) == 0) return(NULL)
    data.frame(
      Comparison = comp,
      Description = head(result$Description, top_n),
      p.adjust = head(result$p.adjust, top_n),
      Count = head(result$Count, top_n)
    )
  })
  
  combined_plot <- ggplot(top_pathways, aes(x = Comparison, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = "Top Enriched KEGG Pathways Across Comparisons",
         size = "Gene Count",
         color = "Adjusted p-value")
  
  ggsave("Combined_KEGG_dotplot.png", combined_plot, width = 14, height = 12, dpi = 300)
  
  return(combined_plot)
}

# Create and save the combined dotplot
combined_plot <- create_combined_dotplot(kegg_results)
print(combined_plot)

# Create a summary table of top pathways for all comparisons
create_summary_table <- function(kegg_results_list, top_n = 5) {
  summary_table <- map_dfr(names(kegg_results_list), function(comp) {
    result <- kegg_results_list[[comp]]$enrichment
    if (nrow(result) == 0) return(NULL)
    data.frame(
      Comparison = comp,
      Pathway = head(result$Description, top_n),
      p.adjust = head(result$p.adjust, top_n),
      Count = head(result$Count, top_n),
      GeneRatio = head(result$GeneRatio, top_n)
    )
  })
  
  write.csv(summary_table, "KEGG_summary_table.csv", row.names = FALSE)
  return(summary_table)
}

# Create and save the summary table
summary_table <- create_summary_table(kegg_results)
print(head(summary_table))

########################################### venn diagram practice



# Install and load required packages
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) install.packages("ggVennDiagram")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(ggVennDiagram)
library(ggplot2)

# Assuming 'all_results' is your dataframe with DESeq2 results

# Function to get genes with significant log2fold change for a comparison
get_sig_genes <- function(data, comparison, lfc_cutoff = 1) {
  column_name <- paste0(comparison, "_log2FoldChange")
  significant_genes <- data[abs(data[[column_name]]) > lfc_cutoff, "external_gene_name"]
  return(significant_genes)
}

# List of comparisons
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

# Get genes with significant log2fold change for each comparison
sig_genes_list <- lapply(comparisons, function(comp) get_sig_genes(all_results, comp))
names(sig_genes_list) <- comparisons

# Create the Venn diagram
venn_plot <- ggVennDiagram(
  sig_genes_list, 
  label_alpha = 0,
  category.names = comparisons
) +
  scale_fill_gradient(low = "white", high = "skyblue") +
  labs(title = "Overlap of Genes with Significant Log2 Fold Change",
       subtitle = paste("Log2 Fold Change cutoff:", 1),
       fill = "Gene Count") +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Save the plot
ggsave("Venn_diagram_significant_genes.png", venn_plot, width = 12, height = 8, dpi = 300, bg = "white")

# Print summary statistics
cat("Number of genes with significant log2fold change in each comparison:\n")
for(comp in names(sig_genes_list)) {
  cat(comp, ": ", length(sig_genes_list[[comp]]), "\n")
}

# Get genes in all comparisons
genes_in_all <- Reduce(intersect, sig_genes_list)
cat("\nNumber of genes with significant log2fold change in all comparisons:", length(genes_in_all), "\n")
if(length(genes_in_all) > 0) {
  cat("Genes with significant log2fold change in all comparisons:\n")
  print(genes_in_all)
}

# Get genes unique to each comparison
unique_genes <- lapply(names(sig_genes_list), function(comp) {
  setdiff(sig_genes_list[[comp]], unlist(sig_genes_list[names(sig_genes_list) != comp]))
})
names(unique_genes) <- names(sig_genes_list)

cat("\nNumber of genes with significant log2fold change unique to each comparison:\n")
for(comp in names(unique_genes)) {
  cat(comp, ": ", length(unique_genes[[comp]]), "\n")
}

# Save the lists of genes to CSV files
for(comp in names(sig_genes_list)) {
  write.csv(data.frame(Gene = sig_genes_list[[comp]]), 
            file = paste0("Significant_genes_", gsub(" ", "_", comp), ".csv"), 
            row.names = FALSE)
}

write.csv(data.frame(Gene = genes_in_all), 
          file = "Genes_significant_in_all_comparisons.csv", 
          row.names = FALSE)

for(comp in names(unique_genes)) {
  write.csv(data.frame(Gene = unique_genes[[comp]]), 
            file = paste0("Unique_genes_", gsub(" ", "_", comp), ".csv"), 
            row.names = FALSE)
}

####################### reactome

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReactomePA")
BiocManager::install("org.Mm.eg.db")  # For mouse genes. Use org.Hs.eg.db for human genes


library(ReactomePA)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(dplyr)

# Assuming 'all_results' is your dataframe with DESeq2 results

# Function to perform Reactome analysis for a single comparison
perform_reactome_analysis <- function(data, comparison, lfc_cutoff = 1) {
  # Filter significant genes
  sig_genes <- data %>%
    filter(abs(!!sym(paste0(comparison, "_log2FoldChange"))) > lfc_cutoff) %>%
    pull(entrezgene_id)
  
  # Perform enrichment analysis
  enrich_result <- enrichPathway(gene = sig_genes,
                                 organism = "mouse",
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 readable = TRUE)
  
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    cat("No significant enrichment found for", comparison, "\n")
    return(NULL)
  }
  
  # Create dotplot
  dotplot <- dotplot(enrich_result, showCategory = 20) +
    ggtitle(paste("Reactome Pathways -", comparison)) +
    theme(axis.text.y = element_text(size = 8))
  
  # Save the dotplot
  ggsave(paste0("Reactome_dotplot_", comparison, ".png"), dotplot, width = 12, height = 10, dpi = 300)
  
  # Save results to CSV
  write.csv(as.data.frame(enrich_result), file = paste0("Reactome_results_", comparison, ".csv"), row.names = FALSE)
  
  return(enrich_result)
}

# List of comparisons
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

# Perform Reactome pathway analysis for each comparison
results_list <- list()
for (comp in comparisons) {
  cat("\nProcessing", comp, "...\n")
  
  # Perform analysis and store results
  results_list[[comp]] <- perform_reactome_analysis(all_results, comp)
  
  if (!is.null(results_list[[comp]])) {
    cat("Top Reactome Pathways for", comp, ":\n")
    print(head(results_list[[comp]]))
  }
}

# Create a combined plot for all comparisons
if (length(results_list) > 0) {
  valid_results <- results_list[!sapply(results_list, is.null)]
  if (length(valid_results) > 0) {
    combined_plot <- compareCluster(valid_results, fun = "enrichPathway", organism = "mouse")
    dotplot <- dotplot(combined_plot, showCategory = 15) +
      ggtitle("Comparison of Enriched Reactome Pathways") +
      theme(axis.text.y = element_text(size = 8))
    ggsave("Combined_Reactome_dotplot.png", dotplot, width = 15, height = 12, dpi = 300)
  } else {
    cat("No valid results for combined plot.\n")
  }
} else {
  cat("No results available for combined plot.\n")
}



#################### noah code for heatmap


# Read the data
# Replace 'your_file.txt' with the actual name of your file
# Load count matrix
Data <- read.csv("C:/Users/epige/OneDrive/Documents/Synovial-rnaseq/Syn_ gene_count.csv", header=T, row.names=1)
# Load required libraries
# Ensure Data is a matrix
Data <- as.matrix(Data)

# Log2 transform the data (add a small value to avoid log(0))
data_log <- log2(Data + 1)

# Calculate Z-scores for each row
data_z <- t(scale(t(data_log)))

# Create annotation data frame
annotation_df <- data.frame(
  Sample = c("F83", "F85", "F87", "F101", "F105", "F109", "F121", "F125", "F127", "F145", "F143", "F147"),
  Group = c("HU_No_Ex", "HU_No_Ex", "HU_No_Ex",
            "NL_No_Ex", "NL_No_Ex", "NL_No_Ex",
            "HU_Ex", "HU_Ex", "HU_Ex",
            "NL_Ex", "NL_Ex", "NL_Ex")
)
rownames(annotation_df) <- annotation_df$Sample

# Split the Group into separate columns
annotation_df <- annotation_df %>%
  separate(Group, c("Condition", "Exercise"), sep = "_")

# Define colors for annotation
ann_colors <- list(
  Condition = c(HU = "#E41A1C", NL = "#4DAF4A"),
  Exercise = c(Ex = "#377EB8", No = "#FF7F00")
)

# Create the heatmap
heatmap <- pheatmap(data_z,
                    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                    show_rownames = FALSE,  # Hide row names due to potentially large number of genes
                    cluster_cols = TRUE,    # Cluster columns (samples)
                    cluster_rows = TRUE,    # Cluster rows (genes)
                    main = "Gene Expression Heatmap",
                    fontsize = 10,
                    cellwidth = 15,
                    cellheight = 0.5,       # Adjust this value based on the number of genes
                    filename = "gene_expression_heatmap.png",
                    width = 10,
                    height = 12)

# If you want to display the heatmap in R as well
print(heatmap)

# Optional: Create a heatmap with top variable genes
# Optional: Create a heatmap with top variable genes
# Calculate variance of each gene
gene_variance <- apply(data_z, 1, var)
top_variable_genes <- names(sort(gene_variance, decreasing = TRUE)[1:50])  # Top 100 most variable genes

# Create heatmap with top variable genes
heatmap_top <- pheatmap(data_z[top_variable_genes,],
                        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                        show_rownames = TRUE,  # Show row names for top genes
                        cluster_cols = TRUE,
                        cluster_rows = TRUE,
                        annotation_col = annotation_df,
                        annotation_colors = ann_colors,
                        main = "Top 100 Variable Genes Heatmap",
                        fontsize = 8,
                        cellwidth = 15,
                        cellheight = 10,
                        filename = "top_variable_genes_heatmap.png",
                        width = 16,
                        height = 18)

print(heatmap_top)

# Save the data used for the heatmap
write.csv(data_z, "heatmap_data.csv")



####### calculating group averages and plotting heatmap

# Load count matrix
data <- read.csv("C:/Users/epige/OneDrive/Documents/Synovial-rnaseq/Syn_ gene_count.csv", header=T, row.names=1)

# Create a named vector for the mapping
old_to_new_names <- c(
  "F83" = "HU_No_Ex", "F85" = "HU_No_Ex", "F87" = "HU_No_Ex",
  "F101" = "NL_No_Ex", "F105" = "NL_No_Ex", "F109" = "NL_No_Ex",
  "F121" = "HU_Ex", "F125" = "HU_Ex", "F127" = "HU_Ex",
  "F145" = "NL_Ex", "F143" = "NL_Ex", "F147" = "NL_Ex"
)

# Create a new dataframe for means
Data_mean <- data.frame(row.names = rownames(data))

# Calculate means for each group
for (group in unique(old_to_new_names)) {
  cols <- names(old_to_new_names)[old_to_new_names == group]
  Data_mean[[paste0(group, "_Mean")]] <- rowMeans(data[, cols, drop = FALSE])
}

# Add the mean columns to the original dataframe
data <- cbind(data, Data_mean)

# Check the new structure
str(data)

# Print the names of the new mean columns
print(names(Data_mean))

# Create a named vector for the mapping
old_to_new_names <- c(
  "F83" = "HU_No_Ex", "F85" = "HU_No_Ex", "F87" = "HU_No_Ex",
  "F101" = "NL_No_Ex", "F105" = "NL_No_Ex", "F109" = "NL_No_Ex",
  "F121" = "HU_Ex", "F125" = "HU_Ex", "F127" = "HU_Ex",
  "F145" = "NL_Ex", "F143" = "NL_Ex", "F147" = "NL_Ex"
)

# Rename the columns
names(data) <- ifelse(names(data) %in% names(old_to_new_names), 
                      old_to_new_names[names(data)], 
                      names(data))
# Using base R
Data_mean <- data[, (ncol(data) - 3):ncol(data)]

# use the Data_means frame to make heatmap

# Log2 transform the data (add a small value to avoid log(0))
data_log <- log2(Data_mean + 1)

# Calculate Z-scores for each row
data_z <- t(scale(t(data_log)))

# Select top 50 most variable genes
gene_variance <- apply(data_z, 1, var)
top_50_genes <- names(sort(gene_variance, decreasing = TRUE)[1:50])
data_z_top50 <- data_z[top_50_genes,]

# Create annotation for the columns
annotation_col <- data.frame(
  Group = sub("_Mean", "", colnames(Data_mean)),
  Condition = substr(colnames(Data_mean), 1, 2),
  Exercise = ifelse(grepl("Ex", colnames(Data_mean)), "Ex", "No_Ex")
)
rownames(annotation_col) <- colnames(Data_mean)

# Define colors for annotation
ann_colors <- list(
  Condition = c(HU = "#E41A1C", NL = "#4DAF4A"),
  Exercise = c(Ex = "#377EB8", No_Ex = "#FF7F00")
)

# Create the heatmap
heatmap <- pheatmap(data_z_top50,
                    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                    show_rownames = TRUE,  # Show row names for top 50 genes
                    cluster_cols = FALSE,  # Don't cluster columns as they are already means
                    cluster_rows = TRUE,   # Cluster rows (genes)
                    annotation_col = annotation_col,
                    annotation_colors = ann_colors,
                    main = "Top 50 Variable Genes Heatmap (Group Means)",
                    fontsize = 8,
                    cellwidth = 25,
                    cellheight = 12,
                    filename = "top_50_variable_genes_heatmap_group_means.png",
                    width = 10,
                    height = 16)

# If you want to display the heatmap in R as well
print(heatmap)

# Save the data used for the heatmap
write.csv(data_z_top50, "heatmap_data_top50_group_means.csv")

# Save the data used for the heatmap
write.csv(Data_mean, "synovial_group_means.csv")





########################### trying new heatmap strategies
# Load count matrix

# Remove duplicate rows based on the first column
heat <- heat[!duplicated(heat[,1]),]

# Replace NA values in the first column with a placeholder
heat[,1] <- make.names(as.character(heat[,1]), unique=TRUE)

# Set row names and remove the first column
rownames(heat) <- heat[,1]
heat <- heat[,-1]

# Log2 transform the data (add a small value to avoid log(0))
heat_log <- log2(heat + 1)

# Calculate Z-scores for each row
heat_z <- t(scale(t(heat_log)))

# Select top 50 most variable genes
gene_variance <- apply(heat_z, 1, var)
top_50_genes <- names(sort(gene_variance, decreasing = TRUE)[1:50])
heat_z_top50 <- heat_z[top_50_genes,]

# Create annotation for the columns
annotation_col <- data.frame(
  Group = colnames(heat),
  Condition = substr(colnames(heat), 1, 2),
  Exercise = ifelse(grepl("Ex", colnames(heat)), "Ex", "No_Ex")
)
rownames(annotation_col) <- colnames(heat)

# Define colors for annotation
ann_colors <- list(
  Condition = c(HU = "#E41A1C", NL = "#4DAF4A"),
  Exercise = c(Ex = "#377EB8", No_Ex = "#FF7F00")
)

# Create the heatmap
heatmap <- pheatmap(heat_z_top50,
                    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                    show_rownames = TRUE,  # Show row names for top 50 genes
                    cluster_cols = FALSE,  # Don't cluster columns as they are already means
                    cluster_rows = TRUE,   # Cluster rows (genes)
                    annotation_col = annotation_col,
                    annotation_colors = ann_colors,
                    main = "Top 50 Variable Genes Heatmap for the (Group Means)",
                    fontsize = 8,
                    cellwidth = 25,
                    cellheight = 12,
                    filename = "top_50_variable_genes_heatmap_for_the_group_means.png",
                    width = 10,
                    height = 16)

# If you want to display the heatmap in R as well
print(heatmap)

# Save the data used for the heatmap
write.csv(heat_z_top50, "heatmap_data_top50_group_means_genes.csv")

# Save the original mean data
write.csv(heat, "synovial_group_means_genes.csv") 



######################################################################################
######### analysis with genes of interest





