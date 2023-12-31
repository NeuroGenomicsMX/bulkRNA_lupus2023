# DE types - Treatment P
# Author: Evelia Coss
# Date: 5 Dec 2023
# Differential expression between 5 cell types in the Lupus condition
# ---

# Cluster setup
# qlogin
# module load r/4.0.2
# R

# --- Libraries ----
library(tximport)
library(DESeq2)

# --- Load data -----
workdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/'
outdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_P/DE_files/'
load(file = paste0(workdir, 'counts/txi.RData')) # metadata, txi, tx2gene
samples <- metadata$sample_ID #165, colnames(txi$counts)

# ---- DE of monoctyes vs moDC (only SLE) ----

# Add comparison
#metadata$Internal <- factor(paste(metadata$Group, metadata$Cell_type, sep="_"))

## --- quality control -----
# Only control from all cell types
SLE_df <- metadata[metadata$Group == "SLE",]
SLE_df$group <-as.factor(SLE_df$Group)
dim(SLE_df)
# 23 *5 = 115 samples, 8 columnas

## ---- Construct DESEQDataSet Object (dds object) -----
# If you have performed transcript quantification (with Salmon, kallisto, RSEM, etc.) 
# you could import the data with tximport, which produces a list, and then you can 
# use DESeqDataSetFromTximport().

# dds <- DESeqDataSetFromTximport(txi = txi,
#                                 colData = metadata,
#                                 design = ~ Group)
# dds

#class: DESeqDataSet
#dim: 29744 165
#metadata(1): version
#assays(2): counts avgTxLength
#rownames(29744): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
#rowData names(0):
#  colnames(165): QR011_0 QR011_1 ... QR108_3 QR108_4
#colData names(7): sample_ID Group ... Cell_type Internal

dim(dds)
# [1] 29744   165

# Select only SLE
dds_SLE <- dds[, SLE_df$sample_ID] # select columns
dim(dds_SLE)
#[1] 29744   115

# to compare Ctrl_moDC vs Ctrl_monocyte, make Ctrl_monocyte the reference level,
# and select the last coefficient

# Add comparison
dds_SLE$Internal <- factor(paste(dds_SLE$Group, dds_SLE$Cell_type, sep="_"))
dds_SLE$Internal <-  relevel(dds_SLE$Internal, ref = 'SLE_monocyte')

# check (only SLE)
unique(dds_SLE$Internal)

# [1] SLE_monocyte SLE_moDC     SLE_moDCIMQ  SLE_tolDC    SLE_tolDCIMQ
#Levels: SLE_monocyte SLE_moDC SLE_moDCIMQ SLE_tolDC SLE_tolDCIMQ

design(dds_SLE) <- ~Internal
design(dds_SLE)
# ~Internal

## --- Differential expression analysis ----
dds_SLE <- DESeq(dds_SLE)
resultsNames(dds_SLE) # lists the coefficients

#[1] "Intercept"
#[2] "Internal_SLE_moDC_vs_SLE_monocyte"
#[3] "Internal_SLE_moDCIMQ_vs_SLE_monocyte"
#[4] "Internal_SLE_tolDC_vs_SLE_monocyte"
#[5] "Internal_SLE_tolDCIMQ_vs_SLE_monocyte"

# --- Contrast 1 (monocytes vs moDC) -----
# results(dds, contrast=c("condition","treated","untreated"))
res_SLE_moVSmoDC <- results(dds_SLE, contrast=c("Internal", "SLE_moDC", "SLE_monocyte"))
res_SLE_moVSmoDC

res_SLE_moVSmoDC_Ordered <- res_SLE_moVSmoDC[order(res_SLE_moVSmoDC$padj),]
summary(res_SLE_moVSmoDC_Ordered)
# out of 23214 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 6957, 30%
# LFC < 0 (down)     : 8088, 35%
# outliers [1]       : 0, 0%
# low counts [2]     : 3105, 13%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Save results
write.csv(res_SLE_moVSmoDC_Ordered, file=paste0(outdir, 'DE_SLE_moVSmoDC.csv'))

# --- Contrast 2 (monocytes vs moDC+IMQ) -----
res_SLE_moVSmoDCIMQ <- results(dds_SLE, contrast=c("Internal", "SLE_moDCIMQ", "SLE_monocyte"))
res_SLE_moVSmoDCIMQ

res_SLE_moVSmoDCIMQ_Ordered <- res_SLE_moVSmoDCIMQ[order(res_SLE_moVSmoDCIMQ$padj),]
summary(res_SLE_moVSmoDCIMQ_Ordered)

# out of 23214 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 6823, 29%
# LFC < 0 (down)     : 8046, 35%
# outliers [1]       : 0, 0%
# low counts [2]     : 2218, 9.6%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Save results
write.csv(res_SLE_moVSmoDCIMQ_Ordered, file=paste0(outdir, 'DE_SLE_moVSmoDCIMQ.csv'))

# --- Contrast 3 (monocytes vs tolDC) -----
res_SLE_moVStolDC <- results(dds_SLE, contrast=c("Internal", "SLE_tolDC", "SLE_monocyte"))
res_SLE_moVStolDC

res_SLE_moVStolDC_Ordered <- res_SLE_moVStolDC[order(res_SLE_moVStolDC$padj),]
summary(res_SLE_moVStolDC_Ordered)

# out of 23214 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 6829, 29%
# LFC < 0 (down)     : 7995, 34%
# outliers [1]       : 0, 0%
# low counts [2]     : 3105, 13%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Save results
write.csv(res_SLE_moVStolDC_Ordered, file=paste0(outdir, 'DE_SLE_moVStolDC.csv'))

# --- Contrast 4 (monocytes vs tolDC+IMQ) -----
res_SLE_moVStolDCIMQ <- results(dds_SLE, contrast=c("Internal", "SLE_tolDCIMQ", "SLE_monocyte"))
res_SLE_moVStolDCIMQ

res_SLE_moVStolDCIMQ_Ordered <- res_SLE_moVStolDCIMQ[order(res_SLE_moVStolDCIMQ$padj),]
summary(res_SLE_moVStolDCIMQ_Ordered)

# out of 23214 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 6853, 30%
# LFC < 0 (down)     : 8083, 35%
# outliers [1]       : 0, 0%
# low counts [2]     : 3548, 15%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Save results
write.csv(res_SLE_moVStolDCIMQ_Ordered, file=paste0(outdir, 'DE_SLE_moVStolDCIMQ.csv'))

# --- Save workspace ----
save(res_SLE_moVSmoDC_Ordered, res_SLE_moVSmoDCIMQ_Ordered, res_SLE_moVStolDC_Ordered, 
     res_SLE_moVStolDCIMQ_Ordered, file = paste0(outdir, 'DE_SLE_results.RData'))

# ---- Information of session ----
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
