# Plot Heatmap for SLE samples by treatment
# Author: Sofia Salazar
# December 2023
# -----

# Usage: plot_heatmap(dds_object, DEresultsFile, metadata, figdir, plotName)

# Arguments

# dds_object = DESeq dataset object
# DEresultsFile = DESeq results in csv
# metadata = metadata of the samples contained in the dds_object
# figdir = output directory for figure
# plotName = name of the heatmap plot to be created
# columnContrast = name of the column with the sample metadata to contrast in heatmap

# ---- Libraries ----

library(ggplot2)
library(ComplexHeatmap)
library(DESeq2)
library(RColorBrewer)
library(circlize)

# ---- Define plotting function ----

plot_heatmap <- function(dds_object, DEresultsFile, metadata, figdir, plotName, columnContrast, clusterRows = F){
  DGE <- read.csv(file = paste0(DEresultsFile), row.names = 'X')
  DGE_ordered <- DGE[order(DGE$log2FoldChange, decreasing = T),]
  
  # ---Compute z-score transformation---
  
  counts <- assay(dds_object, 'counts')
  
  var_non_zero <- apply(counts, 1, var) !=0 #  filter out genes that have low variance across samples
  filtered_counts <- counts[var_non_zero, ]
  
  zscores <- t(scale(t(filtered_counts)))
  dim(zscores)
  
  
  mat <- as.matrix(zscores)
  
  # ---getting log2FoldChange values---
  
  row_names_dge <- rownames(DGE_ordered)
  row_names_mat <- rownames(mat)
  matching_rows <- row_names_dge %in% row_names_mat
  
  filtered_DGE <- DGE_ordered[matching_rows, ]
  l2_val <- as.matrix(filtered_DGE$log2FoldChange)
  
  # ---order samples by dosage
  ordered_samples <- metadata[order(metadata$Dosis), ]
  
  # ---order samples by condition---
  
  ordered_samples <- ordered_samples[order(ordered_samples[[columnContrast]]), ]
  
  
  
  # order columns of count matrix
  ordered_mat <- mat[, ordered_samples$sample_ID]
  colnames(ordered_mat) <- NULL
  rownames(ordered_mat) <- NULL
  
  # Expression color
  col_exp <- colorRamp2(c(min(ordered_mat),0, max(ordered_mat)), c('blue','white', 'red'))
  
  split = data.frame(Sample = ordered_samples[[columnContrast]])
  # samples group color
  print(levels(as.factor(ordered_samples[[columnContrast]]))[1])
  print(levels(as.factor(ordered_samples[[columnContrast]]))[2])
  ha_samples <- HeatmapAnnotation(Group = c(ordered_samples[,columnContrast]),
                                  col = list(Group = c('monocyte' = '#40E0D0', 'tolDCIMQ' = '#a62d6b')))
  # logFC color
  
  # col_logFC <- colorRamp2(c(min(l2_val, na.rm = TRUE), 0, max(l2_val, na.rm = TRUE)), c('#34a624','white','#a709eb'))
  #row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC = col_logFC))
  
  # GC dosage annotation
  gc_val = ordered_samples$Dosis
  levels(as.factor(gc_val))
  
  gc_col =c("0" = "#7b807a","2.5" = "#a1db95", "3" = "#a195db", "5"= "#8d49c4", "8" = "#423485")
  gc_samples <- HeatmapAnnotation(GC_Dosage = gc_val, col = list(GC_Dosage = gc_col))
  
  # Age annotation
  age_val = ordered_samples$Edad
  col_age <- colorRamp2(c(min(age_val), max(age_val)), c('#e6d1be','#bd5a00'))
  age_samples <- HeatmapAnnotation(Age = age_val, col = list(Age = col_age))
  
  heatmap_figure <- Heatmap(ordered_mat, cluster_rows = clusterRows, cluster_columns = F, name = 'Z-score',
                            # left_annotation = row_ha, 
                            top_annotation = c(ha_samples ,age_samples, gc_samples),
                            column_split = split, col = col_exp)
  
  
  png(filename =paste0(figdir, "heatmapZscore_", plotName,".png"), width = 15, height = 10, units = 'cm', res = 300)
  print(heatmap_figure)
  dev.off()
  
  save(figdir, file = paste0(figdir, 'heatmapZcore_', plotName, '.RData'))
  
}

# ----PLOT-----

workdir = '/mnt/Citosina/amedina/lupus/RNA_lupus/'
load(file = paste0(workdir, 'counts/txi.RData'))

dir(paste0(workdir, 'DE_Type_Treatment_P/DE_files'))
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = metadata,
                                design = ~ Group)
SLE_df <- metadata[metadata$Group == "SLE",]
SLE_df$group <-as.factor(SLE_df$Group)
dds_SLE <- dds[, SLE_df$sample_ID] # select columns
dds_SLE$Internal <- factor(paste(dds_SLE$Group, dds_SLE$Cell_type, sep="_"))
dds_SLE$Internal <-  relevel(dds_SLE$Internal, ref = 'SLE_monocyte')

comparisons <- list(c('monocyte', 'moDC'),c('monocyte', 'moDCIMQ'), c('monocyte', 'tolDC'), c('monocyte', 'tolDCIMQ'))
de_files <- c('moVSmoDC', 'moVSmoDCIMQ', 'moVStolDC','moVStolDCIMQ')
figdir <- paste0(workdir, 'DE_Type_Treatment_P/figures/')

# ----monocyte vs moDC---

meta_sub <- SLE_df[SLE_df$Cell_type %in% c(comparisons[[1]]),]
meta_sub <- SLE_df[SLE_df$Cell_type %in% c(comparisons[[1]]),]
meta_sub$Group <- as.character(meta_sub$Group)
dds_sub <- dds_SLE[, meta_sub$sample_ID]
DEresultsFile = paste0(workdir, 'DE_Type_Treatment_P/DE_files/DE_SLE_',de_files[1],'.csv')
plotName = de_files[1]
columnContrast = 'Cell_type'

plot_heatmap(dds_sub, DEresultsFile, meta_sub, figdir, plotName, columnContrast, clusterRows = F)

# ----monocyte vs moDC + IMQ----

meta_sub <- SLE_df[SLE_df$Cell_type %in% c(comparisons[[2]]),]
meta_sub <- SLE_df[SLE_df$Cell_type %in% c(comparisons[[2]]),]
meta_sub$Group <- as.character(meta_sub$Group)
dds_sub <- dds_SLE[, meta_sub$sample_ID]
DEresultsFile = paste0(workdir, 'DE_Type_Treatment_P/DE_files/DE_SLE_',de_files[2],'.csv')
plotName = de_files[2]
columnContrast = 'Cell_type'

plot_heatmap(dds_sub, DEresultsFile, meta_sub, figdir, plotName, columnContrast, clusterRows = F)

# ----monocyte vs tolDC----

meta_sub <- SLE_df[SLE_df$Cell_type %in% c(comparisons[[3]]),]
meta_sub <- SLE_df[SLE_df$Cell_type %in% c(comparisons[[3]]),]
meta_sub$Group <- as.character(meta_sub$Group)
dds_sub <- dds_SLE[, meta_sub$sample_ID]
DEresultsFile = paste0(workdir, 'DE_Type_Treatment_P/DE_files/DE_SLE_',de_files[3],'.csv')
plotName = de_files[3]
columnContrast = 'Cell_type'

plot_heatmap(dds_sub, DEresultsFile, meta_sub, figdir, plotName, columnContrast, clusterRows = F)

# ----monocyte vs tolDC + IMQ----
meta_sub <- SLE_df[SLE_df$Cell_type %in% c(comparisons[[4]]),]
meta_sub <- SLE_df[SLE_df$Cell_type %in% c(comparisons[[4]]),]
meta_sub$Group <- as.character(meta_sub$Group)
dds_sub <- dds_SLE[, meta_sub$sample_ID]
DEresultsFile = paste0(workdir, 'DE_Type_Treatment_P/DE_files/DE_SLE_',de_files[4],'.csv')
plotName = de_files[4]
columnContrast = 'Cell_type'

plot_heatmap(dds_sub, DEresultsFile, meta_sub, figdir, plotName, columnContrast, clusterRows = F)


# -----
# for (i in length(comparisons)){
#   meta_sub <- control_df[control_df$Cell_type %in% c(comparisons[[i]]),]
#   meta_sub$Group <- as.character(meta_sub$Group)
#   dds_sub <- dds_control[, meta_sub$sample_ID]
#   DEresultsFile = paste0(workdir, 'DE_Type_Treatment_C/DE_files/DE_control_',de_files[i],'.csv')
#   plotName = de_files[i]
#   columnContrast = 'Cell_type'
#   
#   plot_heatmap(dds_sub, DEresultsFile, meta_sub, figdir, plotName, columnContrast, clusterRows = F)
#   
#   }

