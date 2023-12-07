# DE cell types
# Author: Sofia Salazar
# Date: 26 Nov 2023
# Differential expression analysis of the 5 cell types in controls and samples
# ---

# Cluster setup
# qlogin
# module load r/4.0.2
# R

# ---
# Libraries

library(DESeq2)
library(tximport)
# ---
# Load data
workdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/'
outdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_P_C_celltypes/DE_files/'
load(file = paste0(workdir, 'counts/txi.RData')) # metadata, txi, tx2gene
#########
# Construct the DESEQDataset Object (dds)

dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = metadata,
                                design = ~Group)
#########
# DE of monoctyes control VS patients
# ---
mo_df <- metadata[metadata$Cell_type == 'monocyte',]
mo_df$Group <-as.factor(mo_df$Group)
dim(mo_df)
# [1] 33  6
# ---

# select dds only monocytes
dds_mo <- dds[,mo_df$sample_ID]
dds_mo$Group <- relevel(dds_mo$Group, ref = 'Ctrl')
unique(dds_mo$Group)
design(dds_mo) <- ~Group
# [1] Ctrl SLE 
# ---
# Add comparison
# dds_mo$Internal <- factor(dds_mo$Group)
# dds_mo$Internal <-  relevel(dds_mo$Internal, ref = 'Ctrl')
# unique(dds_mo$Internal)
# [1] Ctrl SLE 
# Levels: Ctrl SLE
#design(dds_mo) <- ~Internal
# design(dds_mo) # ~Internal
# ---
# Differential expression analysis
dds_mo <- DESeq(dds_mo)
resultsNames(dds_mo) # lists the coefficients
# [1] "Intercept"            "Group_SLE_vs_Ctrl"
res_mo <- results(dds_mo, contrast=c("Group", "SLE", "Ctrl"))
res_mo_Ordered <- res_mo[order(res_mo$padj),]
summary(res_mo_Ordered)
# LFC > 0 (up)       : 121, 0.55%
# LFC < 0 (down)     : 23, 0.11%
# outliers [1]       : 0, 0%
#low counts [2]     : 7903, 36%
write.csv(res_mo_Ordered, file=paste0(outdir, 'DE_monocytes.csv'))

#########
# DE of moDC (group 1)
# ---
moDC_df <- metadata[metadata$Cell_type == "moDC",]
moDC_df$Group <-as.factor(moDC_df$Group)
dim(moDC_df)
# [1] 33  6
dds_moDC <- dds[,moDC_df$sample_ID]
dds_moDC$Group <- relevel(dds_moDC$Group, ref = 'Ctrl')
unique(dds_moDC$Group)
design(dds_moDC) <- ~Group

# ---
# Differential expression analysis

dds_moDC <- DESeq(dds_moDC)
res_moDC <- results(dds_moDC, contrast=c("Group", "SLE", "Ctrl"))
# res_moDC <- results(dds_moDC)


res_moDC_Ordered <- res_moDC[order(res_moDC$padj),]
summary(res_moDC_Ordered)
# out of 21576 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 27, 0.13%
# LFC < 0 (down)     : 154, 0.71%
# outliers [1]       : 0, 0%
# low counts [2]     : 3699, 17%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

write.csv(res_moDC_Ordered, file=paste0(outdir, 'DE_moDCs.csv'))

#######
# DE moDC + IMQ (group 2)

# ---
moDC_IMQ_df <- metadata[metadata$Cell_type == 'moDCIMQ',]
moDC_IMQ_df$Group <-as.factor(moDC_IMQ_df$Group)
dim(moDC_IMQ_df)
# [1] 33  6

dds_moDC_IMQ <- dds[,moDC_IMQ_df$sample_ID]
dds_moDC_IMQ$Group <- relevel(dds_moDC_IMQ$Group, ref = 'Ctrl')
unique(dds_moDC_IMQ$Group)
design(dds_moDC_IMQ) <- ~Group
# ---
# Differential expression analysis

dds_moDC_IMQ <- DESeq(dds_moDC_IMQ)
res_moDC_IMQ <- results(dds_moDC_IMQ, contrast=c("Group", "SLE", "Ctrl"))
res_moDC_IMQ



res_moDC_IMQ_Ordered <- res_moDC_IMQ[order(res_moDC_IMQ$padj),]
summary(res_moDC_IMQ_Ordered)
# out of 21577 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2, 0.0093%
# LFC < 0 (down)     : 10, 0.046%
# outliers [1]       : 0, 0%
# low counts [2]     : 4108, 19%

write.csv(res_moDC_IMQ_Ordered, file=paste0(outdir, 'DE_moDCs_IMQ.csv'))

#######
# DE tolDC (Group 3)
# ---
# quality control
tolDC_df <- metadata[metadata$Cell_type == 'tolDC',]
tolDC_df$Group <-as.factor(tolDC_df$Group)
dim(tolDC_df)
# [1] 33  6


dds_tolDC <- dds[,tolDC_df$sample_ID]
dds_tolDC$Group <- relevel(dds_tolDC$Group, ref = 'Ctrl')
unique(dds_tolDC$Group)
design(dds_tolDC) <- ~Group

# ---
# Differential expression analysis

dds_tolDC <- DESeq(dds_tolDC)
res_tolDC <- results(dds_tolDC, contrast=c("Group", "SLE", "Ctrl"))
res_tolDC

res_tolDC_Ordered <- res_tolDC[order(res_tolDC$padj),]
summary(res_tolDC_Ordered)
# out of 23126 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 13, 0.056%
# LFC < 0 (down)     : 72, 0.31%
# outliers [1]       : 0, 0%
# low counts [2]     : 8020, 35%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


write.csv(res_tolDC_Ordered, file=paste0(outdir, 'DE_tolDCs.csv'))

#########
# DE tolDCs + IMQ (group 4)
# ---
# quality control
tolDC_IMQ_df <- metadata[metadata$Cell_type == 'tolDCIMQ',]
tolDC_IMQ_df$Group <-as.factor(tolDC_IMQ_df$Group)
dim(tolDC_IMQ_df)
# [1] 33  6

dds_tolDC_IMQ <- dds[,tolDC_IMQ_df$sample_ID]
dds_tolDC_IMQ$Group <- relevel(dds_tolDC_IMQ$Group, ref = 'Ctrl')
unique(dds_tolDC_IMQ$Group)
design(dds_tolDC_IMQ) <- ~Group
# ---
# Differential expression analysis

dds_tolDC_IMQ <- DESeq(dds_tolDC_IMQ)
res_tolDC_IMQ <- results(dds_tolDC_IMQ, contrast=c("Group", "SLE", "Ctrl"))
res_tolDC_IMQ


res_tolDC_IMQ_Ordered <- res_tolDC_IMQ[order(res_tolDC_IMQ$padj),]
summary(res_tolDC_IMQ_Ordered)
# out of 21512 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1, 0.0046%
# LFC < 0 (down)     : 2, 0.0093%
# outliers [1]       : 0, 0%
# low counts [2]     : 9, 0.042%
write.csv(res_tolDC_IMQ_Ordered, file=paste0(outdir, 'DE_tolDCs_IMQ.csv'))

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
