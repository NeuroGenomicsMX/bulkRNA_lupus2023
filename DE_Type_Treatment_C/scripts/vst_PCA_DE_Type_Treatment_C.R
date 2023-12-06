# Transform counts and plot PCA
# Author: Evelia Coss
# 5 December 2023
# ----

library(DESeq2)
library(tidyverse)
workdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/'
metadata <- read.csv(file = paste0(workdir, 'metadata/metadata.csv'), header = T)
figdir <- paste0(workdir, 'DE_P_C_celltypes/figures/')
load(paste0(workdir, 'counts/dds_txi.RData'))

# Count normalization
counts <- assays(dds)$counts
vsd <- vst(counts, blind=FALSE)
assays(dds)$vsd <- vsd

# ----- Plot PCA --------

mat <- as.matrix(vsd) # col = samples, rows  = genes
pc <- prcomp(t(mat)) # col = genes , rows = samples

pca_df <- as.data.frame(pc$x)
pca_df$sample <- rownames(pca_df)
meta_subset <- metadata[, c("sample_ID", "Group", "Cell_type")]

# add column with cell name

meta_subset <- meta_subset %>%
  mutate(cell_name = case_when(
    Cell_type == 'monocyte' ~ "monocyte",
    Cell_type == 'moDC' ~ "moDC",
    Cell_type == 'moDCIMQ' ~ "moDC + IMQ",
    Cell_type == 'tolDC' ~ "tolDC",
    Cell_type == 'tolDCIMQ' ~ "tolDC + IMQ",
    TRUE ~ as.character(Cell_type)
  ))

meta_subset$cell_name <- as.factor(meta_subset$cell_name)

pca_df <- merge(pca_df, meta_subset, by.x = "sample", by.y = "sample_ID", all.x = T)
pca_df$Group <- as.factor(pca_df$Group)

levels(pca_df$Group)
# [1] "Ctrl" "SLE" 

variance_explained <- pc$sdev^2 / sum(pc$sdev^2) * 100

# ---
# PCA
# Color by cell_type and shape by group

# COLOR CODE

# GROUP: 0        1         2         3         4
# COLOR: #ff3d54  #7e8f61   #68999c   #9379b5   #E69F00
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_name, shape = Group)) +
  geom_point() + 
  scale_shape_manual(values=c(0, 16)) + 
  scale_color_manual(values = c('#ff3d54','forestgreen','#56B4E9','#9379b5','#E69F00')) +
  labs(title = 'PCA Plot: All samples', 
       x = sprintf('Principal Component 1 (%.2f%% Variance Explained)', variance_explained[1]),
       y = sprintf('Principal Component 2 (%.2f%% Variance Explained)',variance_explained[2]),
       color = "Cell type",
       shape = "Group") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white"))

ggsave(p, file = paste0(figdir,'PCA_grouped.png'), dpi = 300)
ggsave(p, file = '/mnt/Citosina/amedina/ssalazar/RNA_lupus/figures/PCA_grouped.png', dpi = 300)
save(p, file = paste0(figdir, 'PCA_grouped.RData'))
# ---
# Remove monocyte samples
no_mo_df <- metadata[metadata$Cell_type != 'monocyte', ]
no_mo_samples <- no_mo_df$sample_ID
no_mo_mat <- as.matrix(vsd[, no_mo_samples])

pc_no_mo <- prcomp(t(no_mo_mat)) # col = genes , rows = samples
pca_no_mo_df <- as.data.frame(pc_no_mo$x)
pca_no_mo_df$sample <- rownames(pca_no_mo_df)
meta_no_mo_subset <- no_mo_df[, c("sample_ID", "Group", "Cell_type")]

meta_no_mo_subset <- meta_no_mo_subset %>%
  mutate(cell_name = case_when(
    Cell_type == 'moDC' ~ "moDC",
    Cell_type == 'moDCIMQ' ~ "moDC + IMQ",
    Cell_type == 'tolDC' ~ "tolDC",
    Cell_type == 'tolDCIMQ' ~ "tolDC + IMQ",
    TRUE ~ as.character(Cell_type)
  ))

pca_no_mo_df <- merge(pca_no_mo_df, meta_no_mo_subset, by.x = "sample", by.y = "sample_ID", all.x = T)
pca_no_mo_df$Group <- as.factor(pca_no_mo_df$Group)
pca_no_mo_df$cell_name <- as.factor(pca_no_mo_df$cell_name)

levels(pca_no_mo_df$Group)
# [1] "Ctrl" "SLE" 

variance_explained_no_mo <- pc_no_mo$sdev^2 / sum(pc_no_mo$sdev^2) * 100

p2 <- ggplot(pca_no_mo_df, aes(x = PC1, y = PC2, color = cell_name, shape = Group)) +
  geom_point() + 
  scale_shape_manual(values=c(0, 16)) + 
  scale_color_manual(values = c('forestgreen','#56B4E9','#9379b5','#E69F00')) +
  labs(title = 'PCA Plot: Without monocytes', 
       x = sprintf('Principal Component 1 (%.2f%% Variance Explained)', variance_explained_no_mo[1]),
       y = sprintf('Principal Component 2 (%.2f%% Variance Explained)',variance_explained_no_mo[2]),
       color = "Cell type",
       shape = "Group") +
  theme_minimal() + 
  theme(plot.background = element_rect(fill = "white"))

ggsave(p2, file = paste0(figdir,'PCA_no_mo_grouped.png'), dpi = 300)
ggsave(p2, file = '/mnt/Citosina/amedina/ssalazar/RNA_lupus/figures/PCA_no_mo_grouped.png', dpi = 300)
save(p2, file = paste0(figdir, 'PCA_no_mo_grouped.RData'))