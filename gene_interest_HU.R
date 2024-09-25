# Load count matrix
syn1 <- read.csv("C:/Users/epige/OneDrive/Documents/Synovial-rnaseq/DESeq2_full_data_HUgoi.csv", header=T)
library(pheatmap)
library(RColorBrewer)

# Function to make names unique and handle NAs
make_unique <- function(x) {
  x[is.na(x)] <- "Unknown"
  if (length(x) > 0) {
    duplicates <- duplicated(x) | duplicated(x, fromLast = TRUE)
    x[duplicates] <- make.unique(x[duplicates])
  }
  x
}

# Select log2FoldChange columns
log2FC_columns <- c("HU_No_Ex_vs_NL_No_Ex_log2FoldChange",
                    "HU_Ex_vs_NL_Ex_log2FoldChange",
                    "HU_No_Ex_vs_HU_Ex_log2FoldChange",
                    "NL_No_Ex_vs_NL_Ex_log2FoldChange")

heatmap1_data <- syn1[, log2FC_columns]

# Make row names unique and handle NAs
rownames(heatmap1_data) <- make_unique(syn1$external_gene_name)

# Remove any rows with NA values in the log2FC columns
heatmap1_data <- na.omit(heatmap1_data)

# Log transform the data (add a small value to avoid log(0))
heatmap1_data_log <- log2(abs(heatmap1_data) + 1) * sign(heatmap1_data)

# Calculate z-statistic
z_statistic <- function(x) {
  (x - mean(x)) / sd(x)
}

heatmap1_data_z <- t(apply(heatmap1_data_log, 1, z_statistic))

# Select top 50 most variable genes
gene_variance <- apply(heatmap1_data_z, 1, var)
top_50_genes <- names(sort(gene_variance, decreasing = TRUE)[1:50])
heatmap1_data_top50 <- heatmap1_data_z[top_50_genes, ]

# Create annotation for the columns
annotation_col <- data.frame(
  Comparison = c("HU_vs_NL_No_Ex", "HU_vs_NL_Ex", "HU_No_Ex_vs_Ex", "NL_No_Ex_vs_Ex"),
  Condition = c("HU_vs_NL", "HU_vs_NL", "HU", "NL"),
  Exercise = c("No_Ex", "Ex", "No_Ex_vs_Ex", "No_Ex_vs_Ex")
)
rownames(annotation_col) <- colnames(heatmap1_data_top50)

# Define colors for annotation
ann_colors <- list(
  Comparison = c("HU_vs_NL_No_Ex" = "#E41A1C", "HU_vs_NL_Ex" = "#377EB8", 
                 "HU_No_Ex_vs_Ex" = "#4DAF4A", "NL_No_Ex_vs_Ex" = "#FF7F00"),
  Condition = c("HU_vs_NL" = "#A65628", "HU" = "#F781BF", "NL" = "#999999"),
  Exercise = c("No_Ex" = "#66C2A5", "Ex" = "#FC8D62", "No_Ex_vs_Ex" = "#8DA0CB")
)

# Create the heatmap
heatmap <- pheatmap(heatmap1_data_top50,
                    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    cluster_cols = FALSE,
                    cluster_rows = TRUE,
                    scale = "none",  # We've already z-scored the data
                    annotation_col = annotation_col,
                    annotation_colors = ann_colors,
                    main = "Top 50 Variable Genes HU_Log2 Fold Changes",
                    fontsize = 8,
                    cellwidth = 25,
                    cellheight = 12,
                    filename = "top_50_variable_genes_log2FC_HU_heatmap.png",
                    width = 12,
                    height = 16)

# Display the heatmap in R
print(heatmap)

# Save the data used for the heatmap
write.csv(heatmap_data_top50, "heatmap_data_top50genes_HU_condition_log2FC.csv")

# volcano plot 
library(ggplot2)
library(ggrepel)
library(dplyr)

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

# Create and save volcano plots
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

for (comparison in comparisons) {
  plot <- create_volcano_plot(syn1, 
                              paste0(comparison, "_log2FoldChange"), 
                              paste0(comparison, "_pvalue"), 
                              gsub("_", " ", comparison))
  
  # Save the plot as a PNG file
  ggsave(paste0("volcano_HU_plot_", comparison, ".png"), 
         plot = plot, 
         width = 10, 
         height = 8, 
         dpi = 300)
  
  # Print summary statistics
  cat("\nSummary Statistics for", comparison, ":\n")
  
  total_de <- sum(syn1[[paste0(comparison, "_pvalue")]] < 0.05 & abs(syn1[[paste0(comparison, "_log2FoldChange")]]) > 1, na.rm = TRUE)
  up_regulated <- sum(syn1[[paste0(comparison, "_pvalue")]] < 0.05 & syn1[[paste0(comparison, "_log2FoldChange")]] > 1, na.rm = TRUE)
  down_regulated <- sum(syn1[[paste0(comparison, "_pvalue")]] < 0.05 & syn1[[paste0(comparison, "_log2FoldChange")]] < -1, na.rm = TRUE)
  
  cat("Total DE genes:", total_de, "\n")
  cat("Up-regulated:", up_regulated, "\n")
  cat("Down-regulated:", down_regulated, "\n")
  
  top_genes <- syn1 %>%
    filter(!!sym(paste0(comparison, "_pvalue")) < 0.05) %>%
    arrange(desc(abs(!!sym(paste0(comparison, "_log2FoldChange"))))) %>%
    head(10)
  
  cat("Top 10 DE genes:\n")
  print(top_genes$external_gene_name)
}

# Save the data used for the volcano plots
write.csv(syn1, "volcano_plot_HU_data.csv", row.names = FALSE)

################# GSEA

library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

# Function to perform GO analysis for each comparison
perform_go_analysis <- function(data, comparison, log2FC_col, pvalue_col) {
  # Convert gene IDs to Entrez IDs
  entrez_ids <- mapIds(org.Mm.eg.db, 
                       keys = data$external_gene_name, 
                       keytype = "SYMBOL", 
                       column = "ENTREZID")
  
  # Create a named vector of log2 fold changes
  gene_list <- data[[log2FC_col]]
  names(gene_list) <- entrez_ids
  
  # Remove NA values
  gene_list <- gene_list[!is.na(names(gene_list))]
  
  # Sort the list in decreasing order (required for GSEA)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Perform GSEA for BP, MF, and CC
  go_results <- list()
  for (ont in c("BP", "MF", "CC")) {
    go_gsea <- gseGO(geneList = gene_list,
                     OrgDb = org.Mm.eg.db,
                     ont = ont,
                     keyType = "ENTREZID",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose = FALSE)
    
    go_results[[ont]] <- go_gsea
    
    # Plot and save results
    if (nrow(go_gsea) > 0) {
      dotplot(go_gsea, showCategory = 30, title = paste("GO", ont, "GSEA -", comparison))
      ggsave(paste0("go_", tolower(ont), "_gsea_dotplot_", comparison, ".png"), width = 12, height = 20)
      
      emapplot(pairwise_termsim(go_gsea))
      ggsave(paste0("go_", tolower(ont), "_gsea_emap_", comparison, ".png"), width = 12, height = 20)
      
      gseaplot2(go_gsea, geneSetID = 1:min(5, nrow(go_gsea)), title = paste("GO", ont, "GSEA -", comparison))
      ggsave(paste0("go_", tolower(ont), "_gsea_plot_", comparison, ".png"), width = 12, height = 20)
    }
    
    # Save results to CSV files
    write.csv(as.data.frame(go_gsea), paste0("go_", tolower(ont), "_gsea_results_", comparison, ".csv"))
  }
  
  return(go_results)
}

# Perform GO analysis for each comparison
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

results <- list()

for (comparison in comparisons) {
  log2FC_col <- paste0(comparison, "_log2FoldChange")
  pvalue_col <- paste0(comparison, "_pvalue")
  
  results[[comparison]] <- perform_go_analysis(syn1, comparison, log2FC_col, pvalue_col)
  
  # Print top 10 enriched GO terms for each ontology
  for (ont in c("BP", "MF", "CC")) {
    cat(paste("\nTop 30 GO", ont, "GSEA terms for", comparison, ":\n"))
    print(head(results[[comparison]][[ont]], 30))
  }
}

# Save the complete results
saveRDS(results, "go_gsea_HU_complete_results.rds")

############################# REACTOME
library(ReactomePA)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Load the data if not already loaded
# syn1 <- read.csv("C:/Users/epige/OneDrive/Documents/Synovial-rnaseq/DESeq2_full_data_HUgoi.csv", header=T)

# Function to perform Reactome pathway analysis for each comparison
perform_reactome_analysis <- function(data, comparison, log2FC_col, pvalue_col) {
  # Convert gene IDs to Entrez IDs
  entrez_ids <- mapIds(org.Mm.eg.db, 
                       keys = data$external_gene_name, 
                       keytype = "SYMBOL", 
                       column = "ENTREZID")
  
  # Create a named vector of log2 fold changes
  gene_list <- data[[log2FC_col]]
  names(gene_list) <- entrez_ids
  
  # Remove NA values
  gene_list <- gene_list[!is.na(names(gene_list))]
  
  # Sort the list in decreasing order (required for GSEA)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Perform Reactome pathway enrichment analysis
  reactome_enrichment <- enrichPathway(gene = names(gene_list)[abs(gene_list) > 1],
                                       universe = names(gene_list),
                                       organism = "mouse",
                                       pvalueCutoff = 0.05,
                                       pAdjustMethod = "BH")
  
  # Perform GSEA with Reactome pathways
  reactome_gsea <- gsePathway(geneList = gene_list,
                              organism = "mouse",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH")
  
  # Plot and save results
  if (nrow(reactome_enrichment) > 0) {
    dotplot(reactome_enrichment, showCategory = 30, title = paste("Reactome Pathway Enrichment -", comparison))
    ggsave(paste0("reactome_HU_enrichment_dotplot_", comparison, ".png"), width = 12, height = 20)
    
    emapplot(pairwise_termsim(reactome_enrichment))
    ggsave(paste0("reactome_HU_enrichment_map_", comparison, ".png"), width = 12, height = 20)
  }
  
  if (nrow(reactome_gsea) > 0) {
    gseaplot2(reactome_gsea, geneSetID = 1:min(5, nrow(reactome_gsea)), title = paste("Reactome Pathway GSEA -", comparison))
    ggsave(paste0("reactome_HU_gsea_plot_", comparison, ".png"), width = 12, height = 20)
  }
  
  # Save results to CSV files
  write.csv(as.data.frame(reactome_enrichment), paste0("reactome_HU_enrichment_results_", comparison, ".csv"))
  write.csv(as.data.frame(reactome_gsea), paste0("reactome_HU_gsea_results_", comparison, ".csv"))
  
  return(list(enrichment = reactome_enrichment, gsea = reactome_gsea))
}

# Perform Reactome analysis for each comparison
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

results <- list()

for (comparison in comparisons) {
  log2FC_col <- paste0(comparison, "_log2FoldChange")
  pvalue_col <- paste0(comparison, "_pvalue")
  
  results[[comparison]] <- perform_reactome_analysis(syn1, comparison, log2FC_col, pvalue_col)
  
  # Print top 10 enriched pathways
  cat("\nTop 10 enriched Reactome pathways for", comparison, ":\n")
  print(head(results[[comparison]]$enrichment, 30))
  
  cat("\nTop 10 GSEA Reactome pathways for", comparison, ":\n")
  print(head(results[[comparison]]$gsea, 30))
}

# Save the complete results
saveRDS(results, "reactome1_analysis_HU_complete_results.rds")





##################### VENN DIAGRAM
library(ggVennDiagram)
library(ggplot2)
library(dplyr)

# Function to get significant genes
get_sig_genes <- function(data, log2FC_col, pvalue_col, log2FC_threshold = 1, pvalue_threshold = 0.05) {
  data %>%
    filter(abs(!!sym(log2FC_col)) > log2FC_threshold, !!sym(pvalue_col) < pvalue_threshold) %>%
    pull(external_gene_name)
}

# Get significant genes for each comparison
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")
sig_genes <- list()

for (comparison in comparisons) {
  log2FC_col <- paste0(comparison, "_log2FoldChange")
  pvalue_col <- paste0(comparison, "_pvalue")
  sig_genes[[comparison]] <- get_sig_genes(syn1, log2FC_col, pvalue_col)
}

# Create Venn diagram
venn_plot <- ggVennDiagram(
  sig_genes,
  category.names = comparisons,
  label = "count",
  edge_size = 0.5,
  set_size = 4
) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF", "#95D840FF")) +
  theme(legend.position = "right")

# Save the Venn diagram
ggsave("venn_diagram_HU_DEGs.png", venn_plot, width = 10, height = 8, dpi = 300)

# Print overlap information
cat("Number of significant genes in each comparison:\n")
for (comparison in comparisons) {
  cat(comparison, ": ", length(sig_genes[[comparison]]), "\n")
}

cat("\nNumber of genes overlapping in all comparisons: ", 
    length(Reduce(intersect, sig_genes)), "\n")

# Function to get pairwise overlaps
get_pairwise_overlap <- function(list1, list2) {
  length(intersect(list1, list2))
}

cat("\nPairwise overlaps:\n")
for (i in 1:(length(comparisons) - 1)) {
  for (j in (i + 1):length(comparisons)) {
    overlap <- get_pairwise_overlap(sig_genes[[comparisons[i]]], sig_genes[[comparisons[j]]])
    cat(comparisons[i], "and", comparisons[j], ": ", overlap, "\n")
  }
}

# Save overlapping genes to a file
overlapping_genes <- Reduce(intersect, sig_genes)
write.csv(data.frame(gene = overlapping_genes), "genes_overlapping_all_comparisons_HU.csv", row.names = FALSE)

# Save all significant genes for each comparison
all_sig_genes <- unique(unlist(sig_genes))
sig_gene_matrix <- matrix(0, nrow = length(all_sig_genes), ncol = length(comparisons))
colnames(sig_gene_matrix) <- comparisons
rownames(sig_gene_matrix) <- all_sig_genes

for (i in 1:length(comparisons)) {
  sig_gene_matrix[sig_genes[[comparisons[i]]], i] <- 1
}

write.csv(as.data.frame(sig_gene_matrix), "all_significant_genes_HU.csv")



##################### UPSET PLOT

library(UpSetR)
library(dplyr)

# Function to get significant genes
get_sig_genes <- function(data, log2FC_col, pvalue_col, log2FC_threshold = 1, pvalue_threshold = 0.05) {
  data %>%
    filter(abs(!!sym(log2FC_col)) > log2FC_threshold, !!sym(pvalue_col) < pvalue_threshold) %>%
    pull(external_gene_name)
}

# Get significant genes for each comparison
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")
sig_genes <- list()

for (comparison in comparisons) {
  log2FC_col <- paste0(comparison, "_log2FoldChange")
  pvalue_col <- paste0(comparison, "_pvalue")
  sig_genes[[comparison]] <- get_sig_genes(syn1, log2FC_col, pvalue_col)
}

# Create a data frame for UpSetR
all_genes <- unique(unlist(sig_genes))
upset_data <- data.frame(gene = all_genes)

for (comparison in comparisons) {
  upset_data[[comparison]] <- as.integer(upset_data$gene %in% sig_genes[[comparison]])
}

# Create the UpSet plot
upset_plot <- upset(upset_data, 
                    sets = comparisons,
                    order.by = "freq",
                    nsets = length(comparisons),
                    nintersects = NA,
                    point.size = 3,
                    line.size = 1,
                    mainbar.y.label = "Intersection Size",
                    sets.x.label = "Set Size",
                    text.scale = c(1.3, 1.3, 1, 1, 1.3, 1))

# Save the UpSet plot
png("upset_plot_HU_DEGs.png", width = 12, height = 8, units = "in", res = 300)
print(upset_plot)
dev.off()

# Print summary information
cat("Number of significant genes in each comparison:\n")
for (comparison in comparisons) {
  cat(comparison, ": ", length(sig_genes[[comparison]]), "\n")
}

cat("\nNumber of genes overlapping in all comparisons: ", 
    length(Reduce(intersect, sig_genes)), "\n")

# Save overlapping genes to a file
overlapping_genes <- Reduce(intersect, sig_genes)
write.csv(data.frame(gene = overlapping_genes), "genes_overlapping_all_HU_comparisons.csv", row.names = FALSE)

# Save all significant genes with their comparison information
sig_gene_info <- upset_data %>%
  gather(key = "comparison", value = "is_significant", -gene) %>%
  filter(is_significant == 1) %>%
  select(-is_significant)

write.csv(sig_gene_info, "all_significant_genes_HU_with_comparisons.csv", row.names = FALSE)
