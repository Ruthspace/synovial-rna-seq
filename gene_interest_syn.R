# Load count matrix
syn <- read.csv("C:/Users/epige/OneDrive/Documents/Synovial-rnaseq/DESeq2_EXgoi1.csv", header=T)
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

heatmap_data <- syn[, log2FC_columns]

# Make row names unique and handle NAs
rownames(heatmap_data) <- make_unique(syn$external_gene_name)

# Remove any rows with NA values in the log2FC columns
heatmap_data <- na.omit(heatmap_data)

# Log transform the data (add a small value to avoid log(0))
heatmap_data_log <- log2(abs(heatmap_data) + 1) * sign(heatmap_data)

# Calculate z-statistic
z_statistic <- function(x) {
  (x - mean(x)) / sd(x)
}

heatmap_data_z <- t(apply(heatmap_data_log, 1, z_statistic))

# Select top 50 most variable genes
gene_variance <- apply(heatmap_data_z, 1, var)
top_50_genes <- names(sort(gene_variance, decreasing = TRUE)[1:50])
heatmap_data_top50 <- heatmap_data_z[top_50_genes, ]

# Create annotation for the columns
annotation_col <- data.frame(
  Comparison = c("HU_vs_NL_No_Ex", "HU_vs_NL_Ex", "HU_No_Ex_vs_Ex", "NL_No_Ex_vs_Ex"),
  Condition = c("HU_vs_NL", "HU_vs_NL", "HU", "NL"),
  Exercise = c("No_Ex", "Ex", "No_Ex_vs_Ex", "No_Ex_vs_Ex")
)
rownames(annotation_col) <- colnames(heatmap_data_top50)

# Define colors for annotation
ann_colors <- list(
  Comparison = c("HU_vs_NL_No_Ex" = "#E41A1C", "HU_vs_NL_Ex" = "#377EB8", 
                 "HU_No_Ex_vs_Ex" = "#4DAF4A", "NL_No_Ex_vs_Ex" = "#FF7F00"),
  Condition = c("HU_vs_NL" = "#A65628", "HU" = "#F781BF", "NL" = "#999999"),
  Exercise = c("No_Ex" = "#66C2A5", "Ex" = "#FC8D62", "No_Ex_vs_Ex" = "#8DA0CB")
)

# Create the heatmap
heatmap <- pheatmap(heatmap_data_top50,
                    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    cluster_cols = FALSE,
                    cluster_rows = TRUE,
                    scale = "none",  # We've already z-scored the data
                    annotation_col = annotation_col,
                    annotation_colors = ann_colors,
                    main = "Top 50 Variable Genes Log2 Fold Changes",
                    fontsize = 8,
                    cellwidth = 25,
                    cellheight = 12,
                    filename = "top_50_variable_genes_log2FC_heatmap.png",
                    width = 12,
                    height = 16)

# Display the heatmap in R
print(heatmap)

# Save the data used for the heatmap
write.csv(heatmap_data_top50, "heatmap_data_top50genes_firstcondition_log2FC.csv")

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
  plot <- create_volcano_plot(syn, 
                              paste0(comparison, "_log2FoldChange"), 
                              paste0(comparison, "_pvalue"), 
                              gsub("_", " ", comparison))
  
  # Save the plot as a PNG file
  ggsave(paste0("volcano_plot_", comparison, ".png"), 
         plot = plot, 
         width = 10, 
         height = 8, 
         dpi = 300)
  
  # Print summary statistics
  cat("\nSummary Statistics for", comparison, ":\n")
  
  total_de <- sum(syn[[paste0(comparison, "_pvalue")]] < 0.05 & abs(syn[[paste0(comparison, "_log2FoldChange")]]) > 1, na.rm = TRUE)
  up_regulated <- sum(syn[[paste0(comparison, "_pvalue")]] < 0.05 & syn[[paste0(comparison, "_log2FoldChange")]] > 1, na.rm = TRUE)
  down_regulated <- sum(syn[[paste0(comparison, "_pvalue")]] < 0.05 & syn[[paste0(comparison, "_log2FoldChange")]] < -1, na.rm = TRUE)
  
  cat("Total DE genes:", total_de, "\n")
  cat("Up-regulated:", up_regulated, "\n")
  cat("Down-regulated:", down_regulated, "\n")
  
  top_genes <- syn %>%
    filter(!!sym(paste0(comparison, "_pvalue")) < 0.05) %>%
    arrange(desc(abs(!!sym(paste0(comparison, "_log2FoldChange"))))) %>%
    head(10)
  
  cat("Top 10 DE genes:\n")
  print(top_genes$external_gene_name)
}



################# GSEA


library(clusterProfiler)
library(org.Mmu.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

# Assuming 'syn' is your data frame
# Remove the 'X' column if it's empty or not needed
syn <- syn %>% select(-X)

# Function to perform KEGG analysis for each comparison
perform_kegg_analysis <- function(data, comparison, log2FC_col, pvalue_col) {
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
  
  # Perform KEGG pathway enrichment analysis
  kegg_enrichment <- enrichKEGG(gene = names(gene_list)[abs(gene_list) > 1],
                                universe = names(gene_list),
                                organism = 'mmu',
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH")
  
  # Perform GSEA with KEGG pathways
  kegg_gsea <- gseKEGG(geneList = gene_list,
                       organism = 'mmu',
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
  
  # Plot and save results
  if (nrow(kegg_enrichment) > 0) {
    dotplot(kegg_enrichment, showCategory = 20, title = paste("KEGG Pathway Enrichment -", comparison))
    ggsave(paste0("kegg_enrichment_dotplot_", comparison, ".png"), width = 10, height = 8)
    
    emapplot(pairwise_termsim(kegg_enrichment))
    ggsave(paste0("kegg_enrichment_map_", comparison, ".png"), width = 12, height = 10)
  }
  
  if (nrow(kegg_gsea) > 0) {
    gseaplot2(kegg_gsea, geneSetID = 1:min(5, nrow(kegg_gsea)), title = paste("KEGG Pathway GSEA -", comparison))
    ggsave(paste0("kegg_gsea_plot_", comparison, ".png"), width = 12, height = 8)
  }
  
  # Save results to CSV files
  write.csv(as.data.frame(kegg_enrichment), paste0("kegg_enrichment_results_", comparison, ".csv"))
  write.csv(as.data.frame(kegg_gsea), paste0("kegg_gsea_results_", comparison, ".csv"))
  
  return(list(enrichment = kegg_enrichment, gsea = kegg_gsea))
}

# Perform KEGG analysis for each comparison
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

results <- list()

for (comparison in comparisons) {
  log2FC_col <- paste0(comparison, "_log2FoldChange")
  pvalue_col <- paste0(comparison, "_pvalue")
  
  results[[comparison]] <- perform_kegg_analysis(syn, comparison, log2FC_col, pvalue_col)
  
  # Print top 10 enriched pathways
  cat("\nTop 10 enriched KEGG pathways for", comparison, ":\n")
  print(head(results[[comparison]]$enrichment, 10))
  
  cat("\nTop 10 GSEA KEGG pathways for", comparison, ":\n")
  print(head(results[[comparison]]$gsea, 10))
}


############################# REACTOME

library(ReactomePA)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Assuming 'syn' is your data frame
# Remove the 'X' column if it's empty or not needed
syn <- syn %>% select(-X)

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
    dotplot(reactome_enrichment, showCategory = 20, title = paste("Reactome Pathway Enrichment -", comparison))
    ggsave(paste0("reactome_enrichment_dotplot_", comparison, ".png"), width = 12, height = 10)
    
    emapplot(pairwise_termsim(reactome_enrichment))
    ggsave(paste0("reactome_enrichment_map_", comparison, ".png"), width = 12, height = 10)
  }
  
  if (nrow(reactome_gsea) > 0) {
    gseaplot2(reactome_gsea, geneSetID = 1:min(5, nrow(reactome_gsea)), title = paste("Reactome Pathway GSEA -", comparison))
    ggsave(paste0("reactome_gsea_plot_", comparison, ".png"), width = 12, height = 8)
  }
  
  # Save results to CSV files
  write.csv(as.data.frame(reactome_enrichment), paste0("reactome_enrichment_results_", comparison, ".csv"))
  write.csv(as.data.frame(reactome_gsea), paste0("reactome_gsea_results_", comparison, ".csv"))
  
  return(list(enrichment = reactome_enrichment, gsea = reactome_gsea))
}

# Perform Reactome analysis for each comparison
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

results <- list()

for (comparison in comparisons) {
  log2FC_col <- paste0(comparison, "_log2FoldChange")
  pvalue_col <- paste0(comparison, "_pvalue")
  
  results[[comparison]] <- perform_reactome_analysis(syn, comparison, log2FC_col, pvalue_col)
  
  # Print top 10 enriched pathways
  cat("\nTop 10 enriched Reactome pathways for", comparison, ":\n")
  print(head(results[[comparison]]$enrichment, 10))
  
  cat("\nTop 10 GSEA Reactome pathways for", comparison, ":\n")
  print(head(results[[comparison]]$gsea, 10))
}


#################GSEA

library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Assuming 'syn' is your data frame
# Remove the 'X' column if it's empty or not needed
syn <- syn %>% select(-X)

# Function to perform GSEA for GO terms (MF and CC)
perform_go_gsea <- function(data, comparison, log2FC_col, pvalue_col) {
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
  
  # Perform GSEA for Molecular Function (MF)
  gsea_mf <- gseGO(geneList = gene_list,
                   OrgDb = org.Mm.eg.db,
                   ont = "MF",
                   keyType = "ENTREZID",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH")
  
  # Perform GSEA for Cellular Component (CC)
  gsea_cc <- gseGO(geneList = gene_list,
                   OrgDb = org.Mm.eg.db,
                   ont = "CC",
                   keyType = "ENTREZID",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH")
  
  # Plot and save results for MF
  if (nrow(gsea_mf) > 0) {
    dotplot(gsea_mf, showCategory = 20, title = paste("GO Molecular Function GSEA -", comparison))
    ggsave(paste0("go_mf_gsea_dotplot_", comparison, ".png"), width = 12, height = 10)
    
    emapplot(pairwise_termsim(gsea_mf))
    ggsave(paste0("go_mf_gsea_emap_", comparison, ".png"), width = 12, height = 10)
    
    gseaplot2(gsea_mf, geneSetID = 1:min(5, nrow(gsea_mf)), title = paste("GO MF GSEA -", comparison))
    ggsave(paste0("go_mf_gsea_plot_", comparison, ".png"), width = 12, height = 8)
  }
  
  # Plot and save results for CC
  if (nrow(gsea_cc) > 0) {
    dotplot(gsea_cc, showCategory = 20, title = paste("GO Cellular Component GSEA -", comparison))
    ggsave(paste0("go_cc_gsea_dotplot_", comparison, ".png"), width = 12, height = 10)
    
    emapplot(pairwise_termsim(gsea_cc))
    ggsave(paste0("go_cc_gsea_emap_", comparison, ".png"), width = 12, height = 10)
    
    gseaplot2(gsea_cc, geneSetID = 1:min(5, nrow(gsea_cc)), title = paste("GO CC GSEA -", comparison))
    ggsave(paste0("go_cc_gsea_plot_", comparison, ".png"), width = 12, height = 8)
  }
  
  # Save results to CSV files
  write.csv(as.data.frame(gsea_mf), paste0("go_mf_gsea_results_", comparison, ".csv"))
  write.csv(as.data.frame(gsea_cc), paste0("go_cc_gsea_results_", comparison, ".csv"))
  
  return(list(mf = gsea_mf, cc = gsea_cc))
}

# Perform GO GSEA analysis for each comparison
comparisons <- c("HU_No_Ex_vs_NL_No_Ex", "HU_Ex_vs_NL_Ex", "HU_No_Ex_vs_HU_Ex", "NL_No_Ex_vs_NL_Ex")

results <- list()

for (comparison in comparisons) {
  log2FC_col <- paste0(comparison, "_log2FoldChange")
  pvalue_col <- paste0(comparison, "_pvalue")
  
  results[[comparison]] <- perform_go_gsea(syn, comparison, log2FC_col, pvalue_col)
  
  # Print top 10 enriched GO MF terms
  cat("\nTop 10 enriched GO Molecular Function terms for", comparison, ":\n")
  print(head(results[[comparison]]$mf, 10))
  
  # Print top 10 enriched GO CC terms
  cat("\nTop 10 enriched GO Cellular Component terms for", comparison, ":\n")
  print(head(results[[comparison]]$cc, 10))
}


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
  sig_genes[[comparison]] <- get_sig_genes(syn, log2FC_col, pvalue_col)
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
ggsave("venn_diagram_DEGs.png", venn_plot, width = 10, height = 8, dpi = 300)

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
write.csv(Reduce(intersect, sig_genes), "genes_overlapping_all_comparisons.csv", row.names = FALSE)






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
  sig_genes[[comparison]] <- get_sig_genes(syn, log2FC_col, pvalue_col)
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
png("upset_plot_DEGs.png", width = 12, height = 8, units = "in", res = 300)
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
write.csv(Reduce(intersect, sig_genes), "genes_overlapping_all_comparisons.csv", row.names = FALSE)