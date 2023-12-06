# DE types - Treatment C
# Author: Evelia Coss
# Date: 5 Dec 2023
# Differential expression between 5 cell types in the control condition
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
outdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/DE_files/'
load(file = paste0(workdir, 'counts/txi.RData')) # metadata, txi, tx2gene
samples <- metadata$sample_ID #165, colnames(txi$counts)

# ---- DE of monoctyes vs moDC (only control) ----

# Add comparison
# metadata$Internal <- factor(paste(metadata$Group, metadata$Cell_type, sep="_"))

## --- quality control -----
# Only control from all cell types
control_df <- metadata[metadata$Group == "Ctrl",]
control_df$group <-as.factor(control_df$Group)
dim(control_df)
# 50 samples, 8  columns

## ---- Construct DESEQDataSet Object (dds object) -----
# If you have performed transcript quantification (with Salmon, kallisto, RSEM, etc.) 
# you could import the data with tximport, which produces a list, and then you can 
# use DESeqDataSetFromTximport().

names(txi) # columnas en Txi
## [1] "abundance"           "counts"              "length"             
## [4] "countsFromAbundance"

dds <- DESeqDataSetFromTximport(txi = txi,
                                 colData = metadata,
                                 design = ~ Group)
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

# Select only control
dds_control <- dds[, control_df$sample_ID] # select columns
dim(dds_control)
# [1] 29744    50

# to compare Ctrl_moDC vs Ctrl_monocyte, make Ctrl_monocyte the reference level,
# and select the last coefficient

# Add comparison
dds_control$Internal <- factor(paste(dds_control$Group, dds_control$Cell_type, sep="_"))
dds_control$Internal <-  relevel(dds_control$Internal, ref = 'Ctrl_monocyte')

# check (only Control)
unique(dds_control$Internal)
# [1] Ctrl_monocyte Ctrl_moDC     Ctrl_moDCIMQ  Ctrl_tolDC    Ctrl_tolDCIMQ
#Levels: Ctrl_monocyte Ctrl_moDC Ctrl_moDCIMQ Ctrl_tolDC Ctrl_tolDCIMQ

design(dds_control) <- ~Internal
design(dds_control)
# ~Internal

## --- Differential expression analysis ----
dds_control <- DESeq(dds_control)
resultsNames(dds_control) # lists the coefficients

# [1] "Intercept"
# [2] "Internal_Ctrl_moDC_vs_Ctrl_monocyte"
# [3] "Internal_Ctrl_moDCIMQ_vs_Ctrl_monocyte"
# [4] "Internal_Ctrl_tolDC_vs_Ctrl_monocyte"
# [5] "Internal_Ctrl_tolDCIMQ_vs_Ctrl_monocyte"

# --- Contrast 1 (monocytes vs moDC) -----
# results(dds, contrast=c("condition","treated","untreated"))
res_control_moVSmoDC <- results(dds_control, contrast=c("Internal", "Ctrl_moDC", "Ctrl_monocyte"))
res_control_moVSmoDC

# log2 fold change (MLE): Internal Ctrl_moDC vs Ctrl_monocyte
# Wald test p-value: Internal Ctrl moDC vs Ctrl monocyte
# DataFrame with 29744 rows and 6 columns
# baseMean log2FoldChange     lfcSE       stat      pvalue
# <numeric>      <numeric> <numeric>  <numeric>   <numeric>
#  A1BG     1.96421e+02     -0.6632559  0.257652 -2.5742344 1.00462e-02
#A1BG-AS1 1.37586e+02      0.4339117  0.288150  1.5058518 1.32105e-01
#A1CF     7.55939e-02      0.0916079  4.652350  0.0196907 9.84290e-01
#A2M      1.54537e+04      9.7377721  0.464425 20.9673758 1.30253e-97
#A2M-AS1  2.12225e+01     -1.0755376  0.351643 -3.0586044 2.22371e-03
#...              ...            ...       ...        ...         ...
#ZYG11A        1.5605      -1.247414  1.266295   -0.98509 3.24580e-01
#ZYG11B     1249.7482       0.279113  0.126985    2.19800 2.79493e-02
#ZYX       14173.6144      -0.866997  0.217763   -3.98139 6.85146e-05
#ZZEF1      2882.9569      -0.229063  0.138719   -1.65127 9.86838e-02
#ZZZ3        684.8323      -0.178049  0.128343   -1.38728 1.65356e-01
#padj
#<numeric>
#  A1BG     1.86034e-02
#A1BG-AS1 1.86585e-01
#A1CF              NA
#A2M      3.92923e-95
#A2M-AS1  4.66517e-03
#...              ...
#ZYG11A   0.403250968
#ZYG11B   0.047113903
#ZYX      0.000183035
#ZZEF1    0.144975533
#ZZZ3     0.226205204

res_control_moVSmoDC_Ordered <- res_control_moVSmoDC[order(res_control_moVSmoDC$padj),]
summary(res_control_moVSmoDC_Ordered)
#out of 23687 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 5818, 25%
#LFC < 0 (down)     : 6300, 27%
#outliers [1]       : 0, 0%
#low counts [2]     : 4997, 21%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_control_moVSmoDC_Ordered, file=paste0(outdir, 'DE_control_moVSmoDC.csv'))

# --- Contrast 2 (monocytes vs moDC+IMQ) -----
res_control_moVSmoDCIMQ <- results(dds_control, contrast=c("Internal", "Ctrl_moDCIMQ", "Ctrl_monocyte"))
res_control_moVSmoDCIMQ

res_control_moVSmoDCIMQ_Ordered <- res_control_moVSmoDCIMQ[order(res_control_moVSmoDCIMQ$padj),]
summary(res_control_moVSmoDCIMQ_Ordered)

# out of 23687 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 5827, 25%
#LFC < 0 (down)     : 6360, 27%
#outliers [1]       : 0, 0%
#low counts [2]     : 4997, 21%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results


# Save results
write.csv(res_control_moVSmoDCIMQ_Ordered, file=paste0(outdir, 'DE_control_moVSmoDCIMQ.csv'))

# --- Contrast 3 (monocytes vs tolDC) -----
res_control_moVStolDC <- results(dds_control, contrast=c("Internal", "Ctrl_tolDC", "Ctrl_monocyte"))
res_control_moVStolDC

res_control_moVStolDC_Ordered <- res_control_moVStolDC[order(res_control_moVStolDC$padj),]
summary(res_control_moVStolDC_Ordered)

# out of 23687 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 5968, 25% , 5968/23687
#LFC < 0 (down)     : 6100, 26%
#outliers [1]       : 0, 0%
#low counts [2]     : 5450, 23%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

# Save results
write.csv(res_control_moVStolDC_Ordered, file=paste0(outdir, 'DE_control_moVStolDC.csv'))

# --- Contrast 4 (monocytes vs tolDC+IMQ) -----
res_control_moVStolDCIMQ <- results(dds_control, contrast=c("Internal", "Ctrl_tolDCIMQ", "Ctrl_monocyte"))
res_control_moVStolDCIMQ

res_control_moVStolDCIMQ_Ordered <- res_control_moVStolDCIMQ[order(res_control_moVStolDCIMQ$padj),]
summary(res_control_moVStolDCIMQ_Ordered)

#out of 23687 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 5690, 24%
#LFC < 0 (down)     : 6286, 27%
#outliers [1]       : 0, 0%
#low counts [2]     : 4997, 21%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

# Save results
write.csv(res_control_moVStolDCIMQ_Ordered, file=paste0(outdir, 'DE_control_moVStolDCIMQ.csv'))

# --- Save workspace ----
save(res_control_moVSmoDC_Ordered, res_control_moVSmoDCIMQ_Ordered, res_control_moVStolDC_Ordered, 
     res_control_moVStolDCIMQ_Ordered, file = paste0(outdir, 'DE_control_results.RData'))

save(metadata, dds, file = paste0(outdir, 'dds_TreatmentVSCelltype.RData'))


# ---- Information of session ----
sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)

# Matrix products: default
# BLAS:   /cm/shared/apps/r/4.0.2-studio/lib64/R/lib/libRblas.so
# LAPACK: /cm/shared/apps/r/4.0.2-studio/lib64/R/lib/libRlapack.so

# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

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
# [11] tximport_1.18.0
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
# [31] readr_1.4.0            xtable_1.8-4           scales_1.1.1
# [34] cachem_1.0.5           DelayedArray_0.16.3    jsonlite_1.7.2
# [37] annotate_1.68.0        XVector_0.30.0         bit_4.0.4
# [40] hms_1.1.0              ggplot2_3.3.5          dplyr_1.0.10
# [43] grid_4.0.2             cli_3.6.0              tools_4.0.2
# [46] bitops_1.0-7           magrittr_2.0.1         RCurl_1.98-1.3
# [49] RSQLite_2.2.7          tibble_3.1.3           pkgconfig_2.0.3
# [52] crayon_1.4.1           ellipsis_0.3.2         Matrix_1.3-4
# [55] assertthat_0.2.1       httr_1.4.2             R6_2.5.0
# [58] compiler_4.0.2

