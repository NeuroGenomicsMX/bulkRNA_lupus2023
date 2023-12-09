# DE in moDCs and tolDCs in patients
# Author: Sofia Salazar
# December 2023
# Differential expression analysis between moDCs and tolDCs in SLE patients
# ----

# qlogin
# module load r/4.0.2
# R

# ----Load data and libraries----

library(DESeq2)
library(tximport)

workdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/'
outdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_moDC_tolDC_P/DE_files/'
ddsDir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_moDC_tolDC_P/dds_objects/'
load(file = paste0(workdir, 'counts/txi.RData')) # metadata, txi, tx2gene

# ----Filter metadata----

meta_sub <- metadata[(metadata$Group == 'SLE') & (metadata$Cell_type %in% c('moDC', 'tolDC')),]
dim(meta_sub) # [1] 46  6

# ----Construct the DESeqDataset Object (dds)----

dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = metadata,
                                design = ~Group)

# Select only moDCs and tolDCs in patients

dds_sub <- dds[, meta_sub$sample_ID] # select columns
dim(dds_sub) # [1] 29744    46

# Add comparison to compare SLE moDCs VS SLE tolDCs

dds_sub$Internal <- factor(paste(dds_sub$Group, dds_sub$Cell_type, sep="_"))
dds_sub$Internal <-  relevel(dds_sub$Internal, ref = 'SLE_moDC')
unique(dds_sub$Internal)
# [1] SLE_moDC  SLE_tolDC
# Levels: SLE_moDC SLE_tolDC

design(dds_sub) <- ~Internal

# ----Differential expression analysis----

dds_sub <- DESeq(dds_sub)
resultsNames(dds_sub) # lists the coefficients
# [1] "Intercept"                      "Internal_SLE_tolDC_vs_SLE_moDC"
res_DE <- results(dds_sub, contrast = c('Internal', 'SLE_tolDC', 'SLE_moDC'))
res_DE <- res_DE[order(res_DE$padj),]
summary(res_DE)

# out of 21870 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 3793, 17%
# LFC < 0 (down)     : 3856, 18%
# outliers [1]       : 0, 0%
# low counts [2]     : 3333, 15%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# ----Save results----

write.csv(res_DE, file = paste0(outdir, 'DE_moDCs_tolDCs_P.csv'))

# ----Session info----

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
# [31] xtable_1.8-4           scales_1.1.1           cachem_1.0.5          
# [34] DelayedArray_0.16.3    annotate_1.68.0        XVector_0.30.0        
# [37] bit_4.0.4              ggplot2_3.3.5          dplyr_1.0.10          
# [40] grid_4.0.2             cli_3.6.0              tools_4.0.2           
# [43] bitops_1.0-7           magrittr_2.0.1         RCurl_1.98-1.3        
# [46] RSQLite_2.2.7          tibble_3.1.3           pkgconfig_2.0.3       
# [49] crayon_1.4.1           ellipsis_0.3.2         Matrix_1.3-4          
# [52] assertthat_0.2.1       httr_1.4.2             R6_2.5.0              
# [55] compiler_4.0.2  