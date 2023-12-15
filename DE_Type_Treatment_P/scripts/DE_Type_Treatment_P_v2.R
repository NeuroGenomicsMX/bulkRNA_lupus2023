# DE types - Treatment P
# Author: Evelia Coss
# Date: 14 Dec 2023
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

# as.factor
metadata$Age <- as.factor(metadata$Age)
metadata$Dose <- as.factor(metadata$Dose)
metadata$GC <- as.factor(metadata$GC)
metadata$Group <- as.factor(metadata$Group)
metadata$Group_age <- as.factor(metadata$Group_age)
metadata$Nephritis <- as.factor(metadata$Nephritis)
metadata$Other_tx <- as.factor(metadata$Other_tx)
metadata$Cell_type <- as.factor(metadata$Cell_type)

## --- quality control -----
# Only control from all cell types
SLE_df <- metadata[metadata$Group == "SLE",]
dim(SLE_df)
# 23 *5 = 115 samples, 9 columnas

## ---- Construct DESEQDataSet Object (dds object) -----
# If you have performed transcript quantification (with Salmon, kallisto, RSEM, etc.) 
# you could import the data with tximport, which produces a list, and then you can 
# use DESeqDataSetFromTximport().

dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = metadata,
                                design = ~ Group + Cell_type + Group_age)
dds

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

# corregir levels
dds_SLE$Cell_type <-  relevel(dds_SLE$Cell_type, ref = 'monocyte')

# corregir levels
dds_SLE$Dose <-  relevel(dds_SLE$Dose, ref = '0')

# to compare Ctrl_moDC vs Ctrl_monocyte, make Ctrl_monocyte the reference level,
# and select the last coefficient

# Design
design(dds_SLE) <- ~ Group_age + Cell_type + Dose
design(dds_SLE) # ~ Group_age + Cell_type + Dose

## --- Differential expression analysis ----
dds_SLE <- DESeq(dds_SLE)
resultsNames(dds_SLE) # lists the coefficients

# [1] "Intercept"                      "Group_age_more30_vs_less30"
# [3] "Cell_type_moDC_vs_monocyte"     "Cell_type_moDCIMQ_vs_monocyte"
# [5] "Cell_type_tolDC_vs_monocyte"    "Cell_type_tolDCIMQ_vs_monocyte"
# [7] "Dose_2.5_vs_0"                  "Dose_3_vs_0"
# [9] "Dose_5_vs_0"                    "Dose_8_vs_0"

# save
save(metadata, dds_SLE, file = paste0(outdir, 'dds_SLE.RData'))

# --- Contrast 0 (more_30 vs less_30) -----
res_SLE_Group_age <- results(dds_SLE, name="Group_age_more30_vs_less30")
res_SLE_Group_age

res_SLE_Group_age_Ordered <- res_SLE_Group_age[order(res_SLE_Group_age$padj),]
summary(res_SLE_Group_age_Ordered)

#out of 23215 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 300, 1.3%
# LFC < 0 (down)     : 344, 1.5%
# outliers [1]       : 7, 0.03%
# low counts [2]     : 5764, 25%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Save results
write.csv(res_SLE_Group_age_Ordered, file=paste0(outdir, 'DE_SLE_GroupAge_more30_vs_less30.csv'))

# --- Contrast 1 (monocytes vs moDC) -----
# results(dds, contrast=c("condition","treated","untreated"))
res_SLE_moVSmoDC <- results(dds_SLE, contrast=c("Cell_type", "moDC", "monocyte"))
res_SLE_moVSmoDC

res_SLE_moVSmoDC_Ordered <- res_SLE_moVSmoDC[order(res_SLE_moVSmoDC$padj),]
summary(res_SLE_moVSmoDC_Ordered)

# 
# out of 23215 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 6904, 30%
# LFC < 0 (down)     : 8105, 35%
# outliers [1]       : 7, 0.03%
# low counts [2]     : 3104, 13%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_SLE_moVSmoDC_Ordered, file=paste0(outdir, 'DE_SLE_moVSmoDC.csv'))

# --- Contrast 2 (monocytes vs moDC+IMQ) -----
res_SLE_moVSmoDCIMQ <- results(dds_SLE, contrast=c("Cell_type", "moDCIMQ", "monocyte"))
res_SLE_moVSmoDCIMQ

res_SLE_moVSmoDCIMQ_Ordered <- res_SLE_moVSmoDCIMQ[order(res_SLE_moVSmoDCIMQ$padj),]
summary(res_SLE_moVSmoDCIMQ_Ordered)

#
# out of 23215 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 6808, 29%
# LFC < 0 (down)     : 8120, 35%
# outliers [1]       : 7, 0.03%
# low counts [2]     : 2661, 11%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_SLE_moVSmoDCIMQ_Ordered, file=paste0(outdir, 'DE_SLE_moVSmoDCIMQ.csv'))

# --- Contrast 3 (monocytes vs tolDC) -----
res_SLE_moVStolDC <- results(dds_SLE, contrast=c("Cell_type", "tolDC", "monocyte"))
res_SLE_moVStolDC

res_SLE_moVStolDC_Ordered <- res_SLE_moVStolDC[order(res_SLE_moVStolDC$padj),]
summary(res_SLE_moVStolDC_Ordered)

# out of 23215 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 6781, 29%
# LFC < 0 (down)     : 8082, 35%
# outliers [1]       : 7, 0.03%
# low counts [2]     : 3548, 15%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_SLE_moVStolDC_Ordered, file=paste0(outdir, 'DE_SLE_moVStolDC.csv'))

# --- Contrast 4 (monocytes vs tolDC+IMQ) -----
res_SLE_moVStolDCIMQ <- results(dds_SLE, contrast=c("Cell_type", "tolDCIMQ", "monocyte"))
res_SLE_moVStolDCIMQ

res_SLE_moVStolDCIMQ_Ordered <- res_SLE_moVStolDCIMQ[order(res_SLE_moVStolDCIMQ$padj),]
summary(res_SLE_moVStolDCIMQ_Ordered)

# out of 23215 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 6819, 29%
# LFC < 0 (down)     : 8112, 35%
# outliers [1]       : 7, 0.03%
# low counts [2]     : 3548, 15%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Save results
write.csv(res_SLE_moVStolDCIMQ_Ordered, file=paste0(outdir, 'DE_SLE_moVStolDCIMQ.csv'))

# --- Contrast 5 (Dosis 2.5  vs 0) (todos los tipos celulares) -----
res_SLE_Dose2_5VS0 <- results(dds_SLE, name="Dose_2.5_vs_0")
res_SLE_Dose2_5VS0

res_SLE_Dose2_5VS0_Ordered <- res_SLE_Dose2_5VS0[order(res_SLE_Dose2_5VS0$padj),]
summary(res_SLE_Dose2_5VS0_Ordered)

# out of 23215 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1774, 7.6%
# LFC < 0 (down)     : 2350, 10%
# outliers [1]       : 7, 0.03%
# low counts [2]     : 5764, 25%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_SLE_Dose2_5VS0_Ordered, file=paste0(outdir, 'DE_SLE_Dose2_5VS0.csv'))

# --- Contrast 6 (Dosis 3 vs 0) (todos los tipos celulares) -----
res_SLE_Dose3VS0 <- results(dds_SLE, name="Dose_3_vs_0")
res_SLE_Dose3VS0

res_SLE_Dose3VS0_Ordered <- res_SLE_Dose3VS0[order(res_SLE_Dose3VS0$padj),]
summary(res_SLE_Dose3VS0_Ordered)

# out of 23215 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 194, 0.84%
# LFC < 0 (down)     : 486, 2.1%
# outliers [1]       : 7, 0.03%
# low counts [2]     : 7534, 32%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_SLE_Dose3VS0_Ordered, file=paste0(outdir, 'DE_SLE_Dose3VS0.csv'))

# --- Contrast 7 (Dosis 5  vs 0) (todos los tipos celulares) -----
res_SLE_Dose5VS0 <- results(dds_SLE, name="Dose_5_vs_0")
res_SLE_Dose5VS0

res_SLE_Dose5VS0_Ordered <- res_SLE_Dose5VS0[order(res_SLE_Dose5VS0$padj),]
summary(res_SLE_Dose5VS0_Ordered)

# out of 23215 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1234, 5.3%
# LFC < 0 (down)     : 1877, 8.1%
# outliers [1]       : 7, 0.03%
# low counts [2]     : 6205, 27%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_SLE_Dose5VS0_Ordered, file=paste0(outdir, 'DE_SLE_Dose5VS0.csv'))

# --- Contrast 7 (Dosis 8  vs 0) (todos los tipos celulares) -----
res_SLE_Dose8VS0 <- results(dds_SLE, name="Dose_8_vs_0")
res_SLE_Dose8VS0

res_SLE_Dose8VS0_Ordered <- res_SLE_Dose8VS0[order(res_SLE_Dose8VS0$padj),]
summary(res_SLE_Dose8VS0_Ordered)

# out of 23215 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 544, 2.3%
# LFC < 0 (down)     : 1149, 4.9%
# outliers [1]       : 7, 0.03%
# low counts [2]     : 6648, 29%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_SLE_Dose8VS0_Ordered, file=paste0(outdir, 'DE_SLE_Dose8VS0.csv'))

# --- Save workspace ----
save(res_SLE_moVSmoDC_Ordered, res_SLE_moVSmoDCIMQ_Ordered, res_SLE_moVStolDC_Ordered, 
     res_SLE_moVStolDCIMQ_Ordered, res_SLE_Dose8VS0_Ordered, res_SLE_Dose5VS0_Ordered,
     res_SLE_Dose3VS0_Ordered, res_SLE_Dose2_5VS0_Ordered,
     file = paste0(outdir, 'DE_SLE_results.RData'))

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
