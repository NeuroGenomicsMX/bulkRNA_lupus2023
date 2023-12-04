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

# ---
# Load data

indir <- '/mnt/Citosina/amedina/alhernandez/Lupus/DGE_SS/ASE/resultados/star_salmon/'
outdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_P_C_celltypes/DE_files/'
counts <- read.table(paste0(indir, 'salmon.merged.gene_counts.tsv'), sep = '\t', header = T)
dim(counts) # [1] 29744   189
rownames(counts) <- counts$gene_id
colnames(counts)
metadata <- read.csv(file = '/mnt/Citosina/amedina/lupus/RNA_lupus/metadata/metadata.csv', header = T)

#########
# DE of monoctyes 

# ---
# quality control
mo_df <- metadata[metadata$Cell_type == 'monocyte',]
mo_df$group <-as.factor(mo_df$Group)
dim(mo_df)
# [1] 33  7

mo_counts <- counts[, mo_df$sample_ID]
mo_counts <- as.matrix(mo_counts)
dim(mo_counts)
# [1] 29744    33 (Antes 39)

all((mo_df$sample_ID) %in% colnames(mo_counts))
# [1] TRUE

identical(mo_df$sample_ID, colnames(mo_counts))
# [1] TRUE

dds_mo <- DESeqDataSetFromMatrix(countData = round(mo_counts),
                                 colData = mo_df,
                                 design = ~ Group)
dds_mo
# class: DESeqDataSet 
# dim: 29744 33 
# metadata(1): version
# assays(1): counts
# rownames(29744): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
# rowData names(0):
#   colnames(33): QR011_0 QR013_0 ... QR107_0 QR108_0
# colData names(7): sample_ID Group ... Cell_type group

# filter genes with low reads
# keep <- rowSums(counts(dds_mo) >= 10)
# dds_mo <- dds_mo[keep,]
dim(dds_mo)
# [1] 16464    33

dds_mo$Group <- relevel(dds_mo$Group, ref = 'Ctrl')

# ---
# Differential expression analysis

dds_mo <- DESeq(dds_mo)
res_mo <- results(dds_mo)
res_mo

# log2 fold change (MLE): Group SLE vs Ctrl 
# Wald test p-value: Group SLE vs Ctrl 
# DataFrame with 29744 rows and 6 columns
# baseMean log2FoldChange     lfcSE       stat    pvalue      padj
# <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
#   A1BG     195.8398853     -0.1136430  0.182327 -0.6232938  0.533091  0.911528
# A1BG-AS1  76.2669470     -0.0804208  0.182261 -0.4412397  0.659039  0.941983
# A1CF       0.0762707      0.1030405  3.232919  0.0318723  0.974574        NA
# A2M       26.3523048     -0.5190098  0.320718 -1.6182745  0.105603  0.639842
# A2M-AS1   26.1716991     -0.2194198  0.287866 -0.7622300  0.445923  0.873786
# ...              ...            ...       ...        ...       ...       ...
# ZYG11A       3.49082     0.12209481 0.6951456  0.1756392 0.8605774  0.983021
# ZYG11B    1046.59604     0.00600105 0.0900490  0.0666421 0.9468666  0.991971
# ZYX      25387.22971     0.31308386 0.1536245  2.0379817 0.0415518  0.503499
# ZZEF1     2979.29555     0.07260051 0.1135338  0.6394614 0.5225228  0.908964
# ZZZ3       642.03098    -0.02092301 0.0866567 -0.2414472 0.8092086  0.974203


res_mo_Ordered <- res_mo[order(res_mo$padj),]
summary(res_mo_Ordered)
# out of 21857 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 116, 0.53%
# LFC < 0 (down)     : 20, 0.092%
# outliers [1]       : 0, 0%
# low counts [2]     : 4994, 23%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

write.csv(res_mo_Ordered, file=paste0(outdir, 'DE_monocytes.csv'))

#########
# DE of moDC (group 1)
# ---
# quality control
moDC_df <- metadata[metadata$Cell_type == "moDC",]
moDC_df$Group <-as.factor(moDC_df$Group)
dim(moDC_df)
# [1] 33  6

moDC_counts <- counts[, moDC_df$sample_ID]
moDC_counts <- as.matrix(moDC_counts)
dim(moDC_counts)
# [1] 29744    33

all((moDC_df$sample_ID) %in% colnames(moDC_counts))
# [1] TRUE

identical(moDC_df$sample_ID, colnames(moDC_counts))
# [1] TRUE

dds_moDC <- DESeqDataSetFromMatrix(countData = round(moDC_counts),
                                   colData = moDC_df,
                                   design = ~ Group)
dds_moDC
# class: DESeqDataSet 
# dim: 29744 33 
# metadata(1): version
# assays(1): counts
# rownames(29744): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
# rowData names(0):
#   colnames(33): QR011_1 QR013_1 ... QR107_1 QR108_1
# colData names(6): sample_ID Group ... Edad Cell_type


dds_moDC$Group <- relevel(dds_moDC$Group, ref = 'Ctrl')

# ---
# Differential expression analysis

dds_moDC <- DESeq(dds_moDC)
res_moDC <- results(dds_moDC)
res_moDC

# log2 fold change (MLE): Group SLE vs Ctrl 
# Wald test p-value: Group SLE vs Ctrl 
# DataFrame with 29744 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat    pvalue      padj
# <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
#   A1BG     1.59920e+02       0.211916  0.222953  0.950496  0.341860  0.955629
# A1BG-AS1 1.40651e+02       0.348923  0.214270  1.628425  0.103435  0.764383
# A1CF     3.38072e-02      -0.358145  3.232919 -0.110781  0.911790        NA
# A2M      2.69539e+04      -0.278239  0.275066 -1.011533  0.311762  0.942457


res_moDC_Ordered <- res_moDC[order(res_moDC$padj),]
summary(res_moDC_Ordered)
# out of 21576 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 27, 0.13%
# LFC < 0 (down)     : 154, 0.71%
# outliers [1]       : 0, 0%
# low counts [2]     : 4108, 19%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

write.csv(res_moDC_Ordered, file=paste0(outdir, 'DE_moDCs.csv'))

#######
# DE moDC + IMQ (group 2)

# ---
# quality control
moDC_IMQ_df <- metadata[metadata$Cell_type == 'moDCIMQ',]
moDC_IMQ_df$Group <-as.factor(moDC_IMQ_df$Group)
dim(moDC_IMQ_df)
# [1] 33  6

moDC_IMQ_counts <- counts[, moDC_IMQ_df$sample_ID]
moDC_IMQ_counts <- as.matrix(moDC_IMQ_counts)
dim(moDC_IMQ_counts)
# [1] 29744    33

all((moDC_IMQ_df$sample) %in% colnames(moDC_IMQ_counts))
# [1] TRUE

identical(moDC_IMQ_df$sample, colnames(moDC_IMQ_counts))
# [1] TRUE

dds_moDC_IMQ <- DESeqDataSetFromMatrix(countData = round(moDC_IMQ_counts),
                                       colData = moDC_IMQ_df,
                                       design = ~ Group)
dds_moDC_IMQ
# class: DESeqDataSet 
# dim: 29744 33 
# metadata(1): version
# assays(1): counts
# rownames(29744): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
# rowData names(0):
#   colnames(33): QR011_2 QR013_2 ... QR107_2 QR108_2
# colData names(6): sample_ID Group ... Edad Cell_type


dds_moDC_IMQ$Group <- relevel(dds_moDC_IMQ$Group, ref = 'Ctrl')

# ---
# Differential expression analysis

dds_moDC_IMQ <- DESeq(dds_moDC_IMQ)
res_moDC_IMQ <- results(dds_moDC_IMQ)
res_moDC_IMQ



res_moDC_IMQ_Ordered <- res_moDC_IMQ[order(res_moDC_IMQ$padj),]
summary(res_moDC_IMQ_Ordered)


write.csv(res_moDC_IMQ_Ordered, file=paste0(outdir, 'DE_moDCs_IMQ.csv'))

#######
# DE tolDC (Group 3)
# ---
# quality control
tolDC_df <- metadata[metadata$Cell_type == 'tolDC',]
tolDC_df$Group <-as.factor(tolDC_df$Group)
dim(tolDC_df)
# [1] 33  6

tolDC_counts <- counts[, tolDC_df$sample_ID]
tolDC_counts <- as.matrix(tolDC_counts)
dim(tolDC_counts)
# [1] 29744    33

all((tolDC_df$sample) %in% colnames(tolDC_counts))
# [1] TRUE

identical(tolDC_df$sample, colnames(tolDC_counts))
# [1] TRUE

dds_tolDC <- DESeqDataSetFromMatrix(countData = round(tolDC_counts),
                                    colData = tolDC_df,
                                    design = ~ Group)
dds_tolDC



dds_tolDC$Group <- relevel(dds_tolDC$Group, ref = 'Ctrl')

# ---
# Differential expression analysis

dds_tolDC <- DESeq(dds_tolDC)
res_tolDC <- results(dds_tolDC)
res_tolDC



res_tolDC_Ordered <- res_tolDC[order(res_tolDC$padj),]
summary(res_tolDC_Ordered)
# out of 23126 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 17, 0.074%
# LFC < 0 (down)     : 96, 0.42%
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

tolDC_IMQ_counts <- counts[, tolDC_IMQ_df$sample_ID]
tolDC_IMQ_counts <- as.matrix(tolDC_IMQ_counts)
dim(tolDC_IMQ_counts)
# [1] 29744    33

all((tolDC_IMQ_df$sample) %in% colnames(tolDC_IMQ_counts))
# [1] TRUE

identical(tolDC_IMQ_df$sample, colnames(tolDC_IMQ_counts))
# [1] TRUE

dds_tolDC_IMQ <- DESeqDataSetFromMatrix(countData = round(tolDC_IMQ_counts),
                                        colData = tolDC_IMQ_df,
                                        design = ~ Group)
dds_tolDC_IMQ



dds_tolDC_IMQ$Group <- relevel(dds_tolDC_IMQ$Group, ref = 'Ctrl')

# ---
# Differential expression analysis

dds_tolDC_IMQ <- DESeq(dds_tolDC_IMQ)
res_tolDC_IMQ <- results(dds_tolDC_IMQ)
res_tolDC_IMQ

# log2 fold change (MLE): group SLE vs Ctrl 
# Wald test p-value: group SLE vs Ctrl 
# DataFrame with 29744 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat    pvalue      padj
# <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
#   A1BG         171.160      0.0210617  0.166300  0.126649  0.899218  0.999931
# A1BG-AS1     160.711     -0.0846240  0.175739 -0.481533  0.630138  0.999931
# A1CF           0.000             NA        NA        NA        NA        NA
# A2M         9540.525     -0.0733994  0.486095 -0.150998  0.879977  0.999931
# A2M-AS1       22.453      0.4600238  0.350518  1.312412  0.189381  0.999931


res_tolDC_IMQ_Ordered <- res_tolDC_IMQ[order(res_tolDC_IMQ$padj),]
summary(res_tolDC_IMQ_Ordered)

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
