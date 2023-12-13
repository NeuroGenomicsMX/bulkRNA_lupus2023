# GO enrichment -control condition
# Author: Evelia Coss
# Date: 12 Dec 2023
# Make your own UpSetR
# ---
# Cluster setup
# qlogin
# module load r/4.0.2
# R

# --- Libraries ----
library(gprofiler2)
library(enrichplot) # BiocManager::install("enrichplot")
library(DOSE)
library(clusterProfiler) # BiocManager::install("clusterProfiler")
library(ggplot2)
library(tidyverse)
library(dplyr)

# --- load data ----
indir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/DE_files/' 
outdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/DE_files/' 
figdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/figures/'

#load function
functiondir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/'
source(paste0(functiondir, "GOTerm_function.R")) # ConvertData and plot_GOTerm function

# ---- Extraer informacion ----
# Lista de genes diferencialmente expresados
Expressed_list <- ConvertData(indir, abslogFC = 1)
str(Expressed_list)
length(Expressed_list) #4

plot_GOTerm(Expressed_list, outdir= outdir, figdir = figdir)


#####
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
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] forcats_0.5.1          stringr_1.4.0          dplyr_1.0.10          
# [4] purrr_0.3.4            readr_1.4.0            tidyr_1.2.1           
# [7] tibble_3.1.3           tidyverse_1.3.0        ggplot2_3.3.5         
# [10] clusterProfiler_3.18.1 DOSE_3.16.0            enrichplot_1.10.2     
# [13] gprofiler2_0.2.0      

# loaded via a namespace (and not attached):
#   [1] fs_1.5.0             lubridate_1.7.9.2    bit64_4.0.5         
# [4] RColorBrewer_1.1-2   httr_1.4.2           tools_4.0.2         
# [7] backports_1.2.1      utf8_1.2.2           R6_2.5.0            
# [10] DBI_1.1.1            lazyeval_0.2.2       BiocGenerics_0.36.1 
# [13] colorspace_2.0-2     withr_2.4.2          tidyselect_1.2.0    
# [16] gridExtra_2.3        bit_4.0.4            compiler_4.0.2      
# [19] rvest_0.3.6          cli_3.6.0            Biobase_2.50.0      
# [22] xml2_1.3.2           scatterpie_0.1.6     plotly_4.9.3        
# [25] shadowtext_0.0.8     scales_1.1.1         digest_0.6.27       
# [28] pkgconfig_2.0.3      htmltools_0.5.1.1    dbplyr_2.2.1        
# [31] fastmap_1.1.0        readxl_1.3.1         htmlwidgets_1.5.3   
# [34] rlang_1.0.6          rstudioapi_0.13      RSQLite_2.2.7       
# [37] farver_2.1.0         generics_0.1.0       jsonlite_1.7.2      
# [40] BiocParallel_1.24.1  GOSemSim_2.16.1      magrittr_2.0.1      
# [43] GO.db_3.12.1         Matrix_1.3-4         Rcpp_1.0.7          
# [46] munsell_0.5.0        S4Vectors_0.28.1     fansi_0.5.0         
# [49] viridis_0.5.1        lifecycle_1.0.3      stringi_1.6.2       
# [52] ggraph_2.0.5         MASS_7.3-53          plyr_1.8.6          
# [55] qvalue_2.22.0        grid_4.0.2           blob_1.2.1          
# [58] parallel_4.0.2       ggrepel_0.9.1        DO.db_2.9           
# [61] crayon_1.4.1         lattice_0.20-41      graphlayouts_0.7.1  
# [64] haven_2.3.1          cowplot_1.1.1        splines_4.0.2       
# [67] hms_1.1.0            pillar_1.6.2         fgsea_1.16.0        
# [70] igraph_1.2.8         reshape2_1.4.4       stats4_4.0.2        
# [73] fastmatch_1.1-0      reprex_1.0.0         glue_1.4.2          
# [76] downloader_0.4       modelr_0.1.8         data.table_1.14.0   
# [79] BiocManager_1.30.21  vctrs_0.5.1          tweenr_1.0.2        
# [82] cellranger_1.1.0     gtable_0.3.0         polyclip_1.10-0     
# [85] assertthat_0.2.1     cachem_1.0.5         ggforce_0.3.3       
# [88] broom_0.7.9          tidygraph_1.2.0      viridisLite_0.4.0   
# [91] rvcheck_0.1.8        AnnotationDbi_1.52.0 memoise_2.0.0       
# [94] IRanges_2.24.1       ellipsis_0.3.2      