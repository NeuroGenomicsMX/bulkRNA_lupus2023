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
outdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_moDC_tolDC_P/DE_files/'
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
SLE_df <- SLE_df[SLE_df$Cell_type %in% c("moDC", "tolDC"),]
dim(SLE_df)
# 46 samples, 9 columnas

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
#[1] 29744   46

# corregir levels
unique(dds_SLE$Cell_type)
dds_SLE$Cell_type <- factor(dds_SLE$Cell_type, levels = unique(dds_SLE$Cell_type))

# corregir levels
unique(dds_SLE$Dose)
dds_SLE$Dose <- factor(dds_SLE$Dose, levels = unique(dds_SLE$Dose))
dds_SLE$Dose <-  relevel(dds_SLE$Dose, ref = '0')

# to compare Ctrl_moDC vs Ctrl_monocyte, make Ctrl_monocyte the reference level,
# and select the last coefficient

# Design
design(dds_SLE) <- ~ Group_age + Cell_type + Dose
design(dds_SLE)

# Number of participants
table(dds_SLE$Cell_type)

# moDC tolDC
# 23    23

table(dds_SLE$Dose)
# 0   5   8 2.5   3
# 18  20   2   4   2

## --- Differential expression analysis ----
dds_SLE <- DESeq(dds_SLE)
resultsNames(dds_SLE) # lists the coefficients

# Cuando hay este error, significa que debes corregir los levels que hay en las variables
# Error in designAndArgChecker(object, betaPrior) :
#   full model matrix is less than full rank

# [1] "Intercept"                  "Group_age_more30_vs_less30"
# [3] "Cell_type_tolDC_vs_moDC"    "Dose_5_vs_0"
# [5] "Dose_8_vs_0"                "Dose_2.5_vs_0"
# [7] "Dose_3_vs_0"

# --- Contrast 0 (more_30 vs less_30) -----
res_SLE_Group_age_compare <- results(dds_SLE, name="Group_age_more30_vs_less30")
res_SLE_Group_age_compare

res_SLE_Group_age_compare_Ordered <- res_SLE_Group_age_compare[order(res_SLE_Group_age_compare$padj),]
summary(res_SLE_Group_age_compare_Ordered)

# out of 21874 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 8, 0.037%
# LFC < 0 (down)     : 11, 0.05%
# outliers [1]       : 13, 0.059%
# low counts [2]     : 1, 0.0046%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_SLE_Group_age_compare_Ordered, file=paste0(outdir, 'DE_SLE_GroupAge_more30_vs_less30_moDCVStolDC.csv'))

# --- Contrast 1 (moDC vs tolDC) -----
# results(dds, contrast=c("condition","treated","untreated"))
res_SLE_moDCVStolDC_P <- results(dds_SLE, contrast=c("Cell_type", "moDC", "tolDC"))
res_SLE_moDCVStolDC_P

res_SLE_moDCVStolDC_P_Ordered <- res_SLE_moDCVStolDC_P[order(res_SLE_moDCVStolDC_P$padj),]
summary(res_SLE_moDCVStolDC_P_Ordered)

# out of 21874 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 4005, 18%
# LFC < 0 (down)     : 3930, 18%
# outliers [1]       : 13, 0.059%
# low counts [2]     : 3746, 17%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Save results
write.csv(res_SLE_moDCVStolDC_P_Ordered, file=paste0(outdir, 'DE_SLE_moDCs_tolDCs.csv'))

# --- Contrast 2 (Dosis 2.5  vs 0) (solo en moDC y tolDC) -----
res_SLE_Dose2_5VS0_moDCVStolDC <- results(dds_SLE, name="Dose_2.5_vs_0")
res_SLE_Dose2_5VS0_moDCVStolDC

res_SLE_Dose2_5VS0_moDCVStolDC_Ordered <- res_SLE_Dose2_5VS0_moDCVStolDC[order(res_SLE_Dose2_5VS0_moDCVStolDC$padj),]
summary(res_SLE_Dose2_5VS0_moDCVStolDC_Ordered)

# out of 21874 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 356, 1.6%
# LFC < 0 (down)     : 637, 2.9%
# outliers [1]       : 13, 0.059%
# low counts [2]     : 6653, 30%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_SLE_Dose2_5VS0_moDCVStolDC_Ordered, file=paste0(outdir, 'DE_Dose2_5VS0_moDCVStolDC.csv'))

# --- Contrast 3 (Dosis 3  vs 0) (solo en moDC y tolDC) -----
res_SLE_Dose3VS0_moDCVStolDC <- results(dds_SLE, name="Dose_3_vs_0")
res_SLE_Dose3VS0_moDCVStolDC

res_SLE_Dose3VS0_moDCVStolDC_Ordered <- res_SLE_Dose3VS0_moDCVStolDC[order(res_SLE_Dose3VS0_moDCVStolDC$padj),]
summary(res_SLE_Dose3VS0_moDCVStolDC_Ordered)

# out of 21874 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2, 0.0091%
# LFC < 0 (down)     : 31, 0.14%
# outliers [1]       : 13, 0.059%
# low counts [2]     : 1249, 5.7%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Save results
write.csv(res_SLE_Dose3VS0_moDCVStolDC_Ordered, file=paste0(outdir, 'DE_Dose3VS0_moDCVStolDC.csv'))

# --- Contrast 4 (Dosis 5  vs 0) (solo en moDC y tolDC) -----
res_SLE_Dose5VS0_moDCVStolDC <- results(dds_SLE, name="Dose_5_vs_0")
res_SLE_Dose5VS0_moDCVStolDC

res_SLE_Dose5VS0_moDCVStolDC_Ordered <- res_SLE_Dose5VS0_moDCVStolDC[order(res_SLE_Dose5VS0_moDCVStolDC$padj),]
summary(res_SLE_Dose5VS0_moDCVStolDC_Ordered)

# out of 21874 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 163, 0.75%
# LFC < 0 (down)     : 442, 2%
# outliers [1]       : 13, 0.059%
# low counts [2]     : 5408, 25%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_SLE_Dose5VS0_moDCVStolDC_Ordered, file=paste0(outdir, 'DE_Dose5VS0_moDCVStolDC.csv'))

# --- Contrast 5 (Dosis 8  vs 0) (solo en moDC y tolDC) -----
res_SLE_Dose8VS0_moDCVStolDC <- results(dds_SLE, name="Dose_8_vs_0")
res_SLE_Dose8VS0_moDCVStolDC

res_SLE_Dose8VS0_moDCVStolDC_Ordered <- res_SLE_Dose8VS0_moDCVStolDC[order(res_SLE_Dose8VS0_moDCVStolDC$padj),]
summary(res_SLE_Dose8VS0_moDCVStolDC_Ordered)


# out of 21874 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 25, 0.11%
# LFC < 0 (down)     : 79, 0.36%
# outliers [1]       : 13, 0.059%
# low counts [2]     : 6237, 29%
# (mean count < 3)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_SLE_Dose8VS0_moDCVStolDC_Ordered, file=paste0(outdir, 'DE_Dose8VS0_moDCVStolDC.csv'))

# --- Save workspace ----
save(res_SLE_Group_age_compare_Ordered, res_SLE_moDCVStolDC_P_Ordered, 
     res_SLE_Dose2_5VS0_moDCVStolDC_Ordered, res_SLE_Dose3VS0_moDCVStolDC_Ordered, 
     res_SLE_Dose5VS0_moDCVStolDC_Ordered, res_SLE_Dose8VS0_moDCVStolDC_Ordered, 
     file = paste0(outdir, 'DE_SLE_Dose_moDCVStolDC.RData'))

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
