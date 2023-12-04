# DE cell types
# Author: Sofia Salazar
# Date: December 2023
# Differential expression analysis of the 5 cell types in controls and samples
# ---

# Cluster setup
# qlogin
# module load r/4.0.2
# R

# ---
# Libraries

library(DESeq2)

# ---
# Load data

workdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/'

#load dds
load(paste0(workdir, 'counts/dds_txi.RData'))

dim(assays(dds)$counts) # [1] 29744   165
counts <- assays(dds)$counts
colnames(counts)
metadata <- read.csv(file = '/mnt/Citosina/amedina/lupus/RNA_lupus/metadata/metadata.csv', header = T)

#########
# DE of monoctyes 

# ---
print('Monocytes')
mo_df <- metadata[metadata$Cell_type == 'monocyte',]
mo_df$group <-as.factor(mo_df$Group)

mo_counts <- counts[, mo_df$sample_ID]
mo_counts <- as.matrix(mo_counts)
print('Dimension of counts:')
dim(mo_counts)
print('Samples are the same in the metadata and in the count matrix')
all((mo_df$sample_ID) %in% colnames(mo_counts))


identical(mo_df$sample_ID, colnames(mo_counts))
dds_mo <- DESeqDataSetFromMatrix(countData = round(mo_counts),
                                 colData = mo_df,
                                 design = ~ Group)
print('Finished dds')
dds_mo
print('dds dimension:')
dim(dds_mo)
dds_mo$Group <- relevel(dds_mo$Group, ref = 'Ctrl')
# ---
# Differential expression analysis

dds_mo <- DESeq(dds_mo)
res_mo <- results(dds_mo)

res_mo_Ordered <- res_mo[order(res_mo$padj),]
print('DE results summary')
summary(res_mo_Ordered)

write.csv(res_mo_Ordered, file=paste0(workdir, 'DE_P_C_celltypes/DE_files/DE_monocytes.csv'))
print('Saved results')

#########
# DE of moDC (group 1)
# ---
print('moDCs')
moDC_df <- metadata[metadata$Cell_type == "moDC",]
moDC_df$Group <-as.factor(moDC_df$Group)

moDC_counts <- counts[, moDC_df$sample_ID]
moDC_counts <- as.matrix(moDC_counts)
print('Dimension of counts:')
dim(moDC_counts)

print('Samples are the same in the metadata and in the count matrix')
all((moDC_df$sample_ID) %in% colnames(moDC_counts))


identical(moDC_df$sample_ID, colnames(moDC_counts))


dds_moDC <- DESeqDataSetFromMatrix(countData = round(moDC_counts),
                                 colData = moDC_df,
                                 design = ~ Group)
print('Finished dds')

dds_moDC

dds_moDC$Group <- relevel(dds_moDC$Group, ref = 'Ctrl')
print('dds dimension:')
dim(dds_moDC)
# ---
# Differential expression analysis

dds_moDC <- DESeq(dds_moDC)
res_moDC <- results(dds_moDC)
print('DE results summary')
res_moDC_Ordered <- res_moDC[order(res_moDC$padj),]
summary(res_moDC_Ordered)

write.csv(res_moDC_Ordered, file=paste0(workdir, 'DE_P_C_celltypes/DE_files/DE_moDCs.csv'))
print('Saved results')

#######
# DE moDC + IMQ (group 2)
print('moDCs + IMQ')

# ---
# quality control
moDC_IMQ_df <- metadata[metadata$Cell_type == 'moDCIMQ',]
moDC_IMQ_df$Group <-as.factor(moDC_IMQ_df$Group)


moDC_IMQ_counts <- counts[, moDC_IMQ_df$sample_ID]
moDC_IMQ_counts <- as.matrix(moDC_IMQ_counts)
print('Dimension of counts:')

dim(moDC_IMQ_counts)
print('Samples are the same in the metadata and in the count matrix')

all((moDC_IMQ_df$sample) %in% colnames(moDC_IMQ_counts))

identical(moDC_IMQ_df$sample, colnames(moDC_IMQ_counts))

dds_moDC_IMQ <- DESeqDataSetFromMatrix(countData = round(moDC_IMQ_counts),
                                   colData = moDC_IMQ_df,
                                   design = ~ Group)
print('Finished dds')

dds_moDC_IMQ

dds_moDC_IMQ$Group <- relevel(dds_moDC_IMQ$Group, ref = 'Ctrl')

# ---
# Differential expression analysis

dds_moDC_IMQ <- DESeq(dds_moDC_IMQ)
res_moDC_IMQ <- results(dds_moDC_IMQ)


res_moDC_IMQ_Ordered <- res_moDC_IMQ[order(res_moDC_IMQ$padj),]
print('DE results summary')
summary(res_moDC_IMQ_Ordered)

print('Saved results')
write.csv(res_moDC_IMQ_Ordered, file=paste0(workdir, 'DE_P_C_celltypes/DE_files/DE_moDCs_IMQ.csv'))

#######
# DE tolDC (Group 3)
print('tolDCs')

tolDC_df <- metadata[metadata$Cell_type == 'tolDC',]
tolDC_df$Group <-as.factor(tolDC_df$Group)


tolDC_counts <- counts[, tolDC_df$sample_ID]
tolDC_counts <- as.matrix(tolDC_counts)
print('Dimension of counts:')
dim(tolDC_counts)
print('Samples are the same in the metadata and in the count matrix')
all((tolDC_df$sample) %in% colnames(tolDC_counts))


identical(tolDC_df$sample, colnames(tolDC_counts))


dds_tolDC <- DESeqDataSetFromMatrix(countData = round(tolDC_counts),
                                   colData = tolDC_df,
                                   design = ~ Group)
print('Finished dds')
dds_tolDC



dds_tolDC$Group <- relevel(dds_tolDC$Group, ref = 'Ctrl')

# ---
# Differential expression analysis

dds_tolDC <- DESeq(dds_tolDC)
res_tolDC <- results(dds_tolDC)
res_tolDC



res_tolDC_Ordered <- res_tolDC[order(res_tolDC$padj),]
print('DE results summary')
summary(res_tolDC_Ordered)

print('Saved results')
write.csv(res_tolDC_Ordered, file=paste0(workdir, 'DE_P_C_celltypes/DE_files/DE_tolDCs.csv'))

#########
# DE tolDCs + IMQ (group 4)
# ---
print('tolDCs + IMQ')
tolDC_IMQ_df <- metadata[metadata$Cell_type == 'tolDCIMQ',]
tolDC_IMQ_df$Group <-as.factor(tolDC_IMQ_df$Group)
dim(tolDC_IMQ_df)

tolDC_IMQ_counts <- counts[, tolDC_IMQ_df$sample_ID]
tolDC_IMQ_counts <- as.matrix(tolDC_IMQ_counts)
print('Dimension of counts:')

dim(tolDC_IMQ_counts)
print('Samples are the same in the metadata and in the count matrix')

all((tolDC_IMQ_df$sample) %in% colnames(tolDC_IMQ_counts))


identical(tolDC_IMQ_df$sample, colnames(tolDC_IMQ_counts))


dds_tolDC_IMQ <- DESeqDataSetFromMatrix(countData = round(tolDC_IMQ_counts),
                                       colData = tolDC_IMQ_df,
                                       design = ~ Group)
print('Finished dds')
dds_tolDC_IMQ



dds_tolDC_IMQ$Group <- relevel(dds_tolDC_IMQ$Group, ref = 'Ctrl')

# ---
# Differential expression analysis

dds_tolDC_IMQ <- DESeq(dds_tolDC_IMQ)
res_tolDC_IMQ <- results(dds_tolDC_IMQ)

res_tolDC_IMQ_Ordered <- res_tolDC_IMQ[order(res_tolDC_IMQ$padj),]
print('DE results summary')
summary(res_tolDC_IMQ_Ordered)
print('Saved results')
write.csv(res_tolDC_IMQ_Ordered, file=paste0(workdir, 'DE_P_C_celltypes/DE_files/DE_tolDCs_IMQ.csv'))

# ----
sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /cm/shared/apps/r/4.0.2-studio/lib64/R/lib/libRblas.so
# LAPACK: /cm/shared/apps/r/4.0.2-studio/lib64/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     
# 
# other attached packages:
#   [1] DESeq2_1.30.1               SummarizedExperiment_1.20.0
# [3] Biobase_2.50.0              MatrixGenerics_1.2.1       
# [5] matrixStats_1.0.0           GenomicRanges_1.42.0       
# [7] GenomeInfoDb_1.26.7         IRanges_2.24.1             
# [9] S4Vectors_0.28.1            BiocGenerics_0.36.1        
# 
# loaded via a namespace (and not attached):
#   [1] genefilter_1.72.1      locfit_1.5-9.4         tidyselect_1.2.0      
# [4] splines_4.0.2          lattice_0.20-41        generics_0.1.0        
# [7] colorspace_2.0-2       vctrs_0.5.1            utf8_1.2.2            
# [10] blob_1.2.1             XML_3.99-0.6           survival_3.5-5        
# [13] rlang_1.0.6            pillar_1.6.2           glue_1.4.2            
# [16] DBI_1.1.1              BiocParallel_1.24.1    bit64_4.0.5           
# [19] RColorBrewer_1.1-2     GenomeInfoDbData_1.2.4 lifecycle_1.0.3       
# [22] zlibbioc_1.36.0        munsell_0.5.0          gtable_0.3.0          
# [25] memoise_2.0.0          geneplotter_1.68.0     fastmap_1.1.0         
# [28] AnnotationDbi_1.52.0   fansi_0.5.0            Rcpp_1.0.7            
# [31] xtable_1.8-4           scales_1.1.1           cachem_1.0.5          
# [34] DelayedArray_0.16.3    annotate_1.68.0        XVector_0.30.0        
# [37] bit_4.0.4              ggplot2_3.3.5          dplyr_1.0.10          
# [40] grid_4.0.2             cli_3.6.0              tools_4.0.2           
# [43] bitops_1.0-7           magrittr_2.0.1         RCurl_1.98-1.3        
# [46] RSQLite_2.2.7          tibble_3.1.3           pkgconfig_2.0.3       
# [49] crayon_1.4.1           ellipsis_0.3.2         Matrix_1.3-4          
# [52] assertthat_0.2.1       httr_1.4.2             R6_2.5.0              
# [55] compiler_4.0.2  
