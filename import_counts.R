# Import salmon counts
# Author: Sofia Salazar
# Date: December 2023
# Import transcript abundance with tximport
# ---

library(tximport)
library(DESeq2)

workdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/'
metadata <- read.csv(file = paste0(workdir, 'metadata/metadata.csv'), header = T)
samples <- metadata$sample_ID
dir <- '/mnt/Citosina/amedina/alhernandez/Lupus/DGE_SS/ASE/resultados'

files <- file.path(dir, "star_salmon", samples, "quant.sf")
names(files) <- samples
all(file.exists(files))
# [1] TRUE
tx2gene <- read.table(paste0(dir, '/star_salmon/salmon_tx2gene.tsv'), sep = '\t')
colnames(tx2gene) <- c('TXNAME', 'GENEID', 'GENEID2')
head(tx2gene, 20)

# gene-level summarization
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)

head(txi$counts)

metadata$Group <- as.factor(metadata$Group)
rownames(metadata) <- colnames(txi$counts)
all(rownames(metadata) == metadata$sample_ID)


dds <- DESeqDataSetFromTximport(txi, metadata, ~Group)
head(assays(dds)$counts)

save(dds, file = paste0(workdir, 'counts/dds_txi.RData'))

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
