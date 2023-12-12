# Plot Heatmap for cell types DE between SLE patients and controls
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
  

  print(identical(colnames(ordered_mat), ordered_samples$sample_ID))
  
  colnames(ordered_mat) <- NULL
  rownames(ordered_mat) <- NULL
  
  # Expression color
  col_exp <- colorRamp2(c(min(ordered_mat),0, max(ordered_mat)), c('blue','white', 'red'))
  
  split = data.frame(Sample = ordered_samples[[columnContrast]])
  # samples group color
  print(levels(as.factor(ordered_samples[[columnContrast]]))[1])
  print(levels(as.factor(ordered_samples[[columnContrast]]))[2])
  ha_samples <- HeatmapAnnotation(Group = c(ordered_samples[,columnContrast]),
                                  col = list(Group = c('Ctrl' = '#40E0D0', 'SLE' = '#a62d6b')))
  # logFC color
  
  # col_logFC <- colorRamp2(c(min(l2_val, na.rm = TRUE), 0, max(l2_val, na.rm = TRUE)), c('#34a624','white','#a709eb'))
  #row_ha <- rowAnnotation(log2FC = l2_val, col = list(log2FC = col_logFC))
  
  # GC dosage annotation
  gc_val = ordered_samples$Dosis
  print(levels(as.factor(gc_val)))
  levels(as.factor(gc_val))
  
  gc_col =c("0" = "#7b807a","2.5" = "#a1db95", "3" = "#a195db", "5"= "#8d49c4", "8" = "#423485")
  gc_samples <- HeatmapAnnotation(GC_Dosage = gc_val, col = list(GC_Dosage = gc_col))
  
  # Age annotation
  age_val = ordered_samples$Edad
  col_age <- colorRamp2(c(min(age_val), max(age_val)), c('#e6d1be','#bd5a00'))
  age_samples <- HeatmapAnnotation(Age = age_val, col = list(Age = col_age))
  
  heatmap_figure <- Heatmap(ordered_mat, cluster_rows = clusterRows, cluster_columns = F, name = 'Z-score',
                            # left_annotation = row_ha, 
                            top_annotation = c(ha_samples,age_samples,gc_samples),
                            column_split = split, col = col_exp)
  
  
  png(filename =paste0(figdir, "heatmapZscore_", plotName,".png"), width = 15, height = 10, units = 'cm', res = 300)
  print(heatmap_figure)
  dev.off()
  
  save(figdir, file = paste0(figdir, 'heatmapZcore_', plotName, '.RData'))
  
}

# ----PLOT-----

workdir = '/mnt/Citosina/amedina/lupus/RNA_lupus/'
load(file = paste0(workdir, 'counts/txi.RData'))


# -----Monocytes-----

DEresultsFile = paste0(workdir, 'DE_P_C_celltypes/DE_files/DE_monocytes.csv')

load(paste0(workdir, 'DE_P_C_celltypes/dds_objects/dds_mo.RData')) # dds_mo

dds_mo
figdir <- paste0(workdir, 'DE_P_C_celltypes/figures/')
plotName = 'monocytes'
mo_meta <- metadata[metadata$Cell_type == 'monocyte',]
mo_meta$Group <- as.character(mo_meta$Group)
columnContrast <- 'Group'


plot_heatmap(dds_mo, DEresultsFile, mo_meta, figdir, plotName = 'monocytes', columnContrast, clusterRows = F)

# -----moDCs-----

DEresultsFile = paste0(workdir, 'DE_P_C_celltypes/DE_files/DE_moDCs.csv')

load(paste0(workdir, 'DE_P_C_celltypes/dds_objects/dds_moDCs.RData')) # dds_mo

dds_moDC
figdir <- paste0(workdir, 'DE_P_C_celltypes/figures/')
plotName = 'moDCs'
moDC_meta <- metadata[metadata$Cell_type == 'moDC',]
moDC_meta$Group <- as.character(moDC_meta$Group)
columnContrast <- 'Group'

plot_heatmap(dds_moDC, DEresultsFile, moDC_meta, figdir, plotName, columnContrast, clusterRows = F)

# -----tolDCs-----

DEresultsFile = paste0(workdir, 'DE_P_C_celltypes/DE_files/DE_tolDCs.csv')

load(paste0(workdir, 'DE_P_C_celltypes/dds_objects/dds_tolDCs.RData')) # dds_mo

dds_tolDC
figdir <- paste0(workdir, 'DE_P_C_celltypes/figures/')
plotName = 'tolDCs'
tolDC_meta <- metadata[metadata$Cell_type == 'tolDC',]
tolDC_meta$Group <- as.character(tolDC_meta$Group)
columnContrast <- 'Group'

plot_heatmap(dds_tolDC, DEresultsFile, tolDC_meta, figdir, plotName = 'tolDCs' , columnContrast, clusterRows = F)

# -----tolDCs + IMQ-----

DEresultsFile = paste0(workdir, 'DE_P_C_celltypes/DE_files/DE_tolDCs_IMQ.csv')

load(paste0(workdir, 'DE_P_C_celltypes/dds_objects/dds_tolDCs_IMQ.RData')) # dds_mo

dds_tolDC_IMQ
figdir <- paste0(workdir, 'DE_P_C_celltypes/figures/')
plotName = 'tolDCs_IMQ'
tolDC_IMQ_meta <- metadata[metadata$Cell_type == 'tolDCIMQ',]
tolDC_IMQ_meta$Group <- as.character(tolDC_IMQ_meta$Group)
columnContrast <- 'Group'

plot_heatmap(dds_tolDC_IMQ, DEresultsFile, tolDC_IMQ_meta, figdir, plotName = 'tolDCs_IMQ' , columnContrast, clusterRows = F)

# -----moDCs + IMQ-----

DEresultsFile = paste0(workdir, 'DE_P_C_celltypes/DE_files/DE_moDCs_IMQ.csv')

load(paste0(workdir, 'DE_P_C_celltypes/dds_objects/dds_moDCs_IMQ.RData')) # dds_mo

dds_moDC_IMQ
figdir <- paste0(workdir, 'DE_P_C_celltypes/figures/')
plotName = 'moDCs_IMQ'
moDC_IMQ_meta <- metadata[metadata$Cell_type == 'moDCIMQ',]
moDC_IMQ_meta$Group <- as.character(moDC_IMQ_meta$Group)
columnContrast <- 'Group'
plot_heatmap(dds_moDC_IMQ, DEresultsFile, moDC_IMQ_meta, figdir, plotName = 'moDCs_IMQ' , columnContrast, clusterRows = F)

