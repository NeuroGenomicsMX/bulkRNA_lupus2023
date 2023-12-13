# GO enrichment -function
# Author: Evelia Coss y Sofia salazar
# Date: 12 Dec 2023
# Procesa la salida de DESEq (res), extrae los genes diferencialmente expresados y
# Rdetermina los terminos GO, obteniendo graficas de los terminos GO y un manhatan plot de los procesos biologicos.

# ----- ConvertData (function) -----
# Extracts the information and converts it into a list, and also classifies the differentially expressed 
# genes according to the pvalue of 0.05 and the LogFoldChange (abslogFC).

# Usage: ConvertData(indir, abslogFC = 1)

# Arguments
# indir = input directory, csv files must start with `DE_` and end with cn `.csv`, example: "DE_control_moVSmoDC.csv"
# abslogFC = log2FoldChange cutoff value, default = 1

# output = contrasts_list

# ----- ConvertData (function) -----

# Usage: plot_GOTerm(contrasts_list, figdir = figdir, outdir = outdir)

# Arguments
# contrasts_list = output of ConvertData (function)
# figdir =  output directory for figures
# outdir =  output directory for files

# output = Manthatan plot, 

# --- Libraries ----
library(gprofiler2)
library(enrichplot) # BiocManager::install("enrichplot")
library(DOSE)
library(clusterProfiler) # BiocManager::install("clusterProfiler")
library(ggplot2)
library(tidyverse)
library(dplyr)

# ---- Example ------------------
#indir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/DE_files/' 
#outdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/DE_files/' 
#figdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/figures/'

# Expressed_list <- ConvertData(indir, abslogFC = 1)
# str(Expressed_list)
# length(Expressed_list)

# plot_GOTerm(Expressed_list, outdir= outdir, figdir = figdir)

# ---- Extraer informacion ----

ConvertData <- function(indir, abslogFC = 1){ 
  contrasts_list <- list() # create a empty list
  # --- Load data -----
  # change directory
  setwd(indir)
  # change to place where "DE_plotName.csv" files are:
  files <- dir(getwd(), pattern = "^DE_(.+)\\.csv$") #only CSV
  
  # --- multiple files ----
  for (f in files){
    print(f)
    plot_name <- gsub("^DE_(.+)\\.csv$", "\\1", f) #name
    
    # Load data
    df <- read.csv(file = f, row.names = 'X')
    df <- na.omit(df)
    
    # --
    # 1. Add expression column
    df$genes <- row.names(df) # Gene names
    df <- df %>% 
      dplyr::mutate(Expression = case_when(log2FoldChange >= abslogFC & padj < 0.05 ~ "Up-regulated",
                                           log2FoldChange <= -(abslogFC) & padj < 0.05 ~ "Down-regulated",
                                           TRUE ~ "Unchanged")) 
    
    # Up names
    up <- df %>% filter(Expression == 'Up-regulated') %>% 
      arrange(padj, desc(abs(log2FoldChange)))
    
    down <- df %>% filter(Expression == 'Down-regulated') %>% 
      arrange(padj, desc(abs(log2FoldChange)))
    
    # Results
    tmp <- list(UP = up, DOWN = down); 
    contrasts_list[[plot_name]] <- tmp
    
  } #end loop
  return(contrasts_list)
  
} # end function


# outfile, example
# $ control_moVStolDCIMQ:List of 2
# ..$ UP  :'data.frame':	3099 obs. of  8 variables:
#   .. ..$ baseMean      : num [1:3099] 89209 5727 87429 4123 3685 ...
# .. ..$ log2FoldChange: num [1:3099] 9.58 8.18 10.23 8.76 9.2 ...
# .. ..$ lfcSE         : num [1:3099] 0.233 0.231 0.292 0.25 0.269 ...
# .. ..$ stat          : num [1:3099] 41 35.4 35.1 35 34.2 ...
# .. ..$ pvalue        : num [1:3099] 0.00 2.17e-274 3.78e-269 2.47e-268 3.33e-256 ...
# .. ..$ padj          : num [1:3099] 0.00 2.03e-270 2.36e-265 1.16e-264 1.24e-252 ...
# .. ..$ genes         : chr [1:3099] "MRC1" "SLC7A8" "GPNMB" "PPAP2B" ...
# .. ..$ Expression    : chr [1:3099] "Up-regulated" "Up-regulated" "Up-regulated" "Up-regulated" ...
# .. ..- attr(*, "na.action")= 'omit' Named int [1:11041] 18704 18705 18706 18707 18708 18709 18710 18711 18712 18713 ...
# .. .. ..- attr(*, "names")= chr [1:11041] "A1CF" "A2ML1" "A2MP1" "A3GALT2" ...
# ..$ DOWN:'data.frame':	3589 obs. of  8 variables:
#   .. ..$ baseMean      : num [1:3589] 2616 3316 5508 1080 3268 ...
# .. ..$ log2FoldChange: num [1:3589] -3.9 -2.23 -1.74 -5.68 -2.88 ...
# .. ..$ lfcSE         : num [1:3589] 0.1299 0.0811 0.0691 0.2475 0.131 ...
# .. ..$ stat          : num [1:3589] -30 -27.5 -25.2 -23 -22 ...
# .. ..$ pvalue        : num [1:3589] 3.87e-198 2.02e-166 1.56e-140 1.26e-116 4.43e-107 ...
# .. ..$ padj          : num [1:3589] 4.83e-195 1.50e-163 8.59e-138 4.79e-114 1.31e-104 ...
# .. ..$ genes         : chr [1:3589] "RBM3" "BIN2" "EIF2S3" "VNN2" ...
# .. ..$ Expression    : chr [1:3589] "Down-regulated" "Down-regulated" "Down-regulated" "Down-regulated" ...
# .. ..- attr(*, "na.action")= 'omit' Named int [1:11041] 18704 18705 18706 18707 18708 18709 18710 18711 18712 18713 ...
# .. .. ..- attr(*, "names")= chr [1:11041] "A1CF" "A2ML1" "A2MP1" "A3GALT2" ...


# ---- Extraer informacion ----

plot_GOTerm <- function(contrasts_list, figdir = figdir, outdir = outdir){ 
  
  for(c in 1:length(contrasts_list)) {
  plot_name <- names(contrasts_list)[c] # name
  
  print(plot_name)
  # Constrats - list
  multi_gp <- gost(list("Upregulated" = Expressed_list[[c]]$UP$genes, "Downregulated" = Expressed_list[[c]]$DOWN$genes), 
                   correction_method = "fdr", multi_query = F, ordered_query = T, 
                   organism = 'hsapiens')
  
  gost_query <- as.data.frame(multi_gp$result)
  
  ## ---- colors ---
  # paleta de colores
  Category_colors <- data.frame(
    category = c("GO:BP", "GO:CC", "GO:MP", "KEGG",
                 'REAC', 'TF', 'CORUM', 'MIRNA', 'HPA', 'HP', 'WP'), 
    label = c('Biological Process', 'Cellular Component', 'Molecular Function',  "KEGG",
              'REAC', 'TF', 'CORUM', 'MIRNA', 'HPA', 'HP', 'WP'),
    colors =  c('#3C6997', '#DD7230','#B4DC7F', '#081d58',
                '#d9c621','#25ced1', 'purple', '#b30039', '#dd3497', '#fb6a4a', '#4daf4a'))
  
  ## ----manhattan plot--------
  print("Create Manhattan plot")
  gostp1 <- gostplot(multi_gp, interactive = FALSE)
  
  # save Graph
  print("Save Manhattan plot")
  ggsave(paste0(figdir, "ManhattanGO_", plot_name, ".png"),
         plot = gostp1, dpi = 300)
  
  # barplots
  print("All GO terms")
  bar_data <- data.frame("term" = as.factor(gost_query$term_name), "condition" = gost_query$query, 
                         "count" = gost_query$term_size, "p.adjust" = gost_query$p_value, 'category' = as.factor(gost_query$source))
  
  top_terms <- head(bar_data[order(bar_data$p.adjust),],40)
  top_terms <- subset(top_terms, p.adjust < 1e-20)
  
  # save dataset
  print("Save TOP Terms dataset in CSV")
  write.csv(top_terms, file = paste0(outdir, "Gprofiler_TopTerms_", plot_name, ".csv"))
  
  bar_data_up <- subset(bar_data, condition == 'Upregulated')
  bar_data_up <- head(bar_data_up[order(bar_data_up$p.adjust),],15)
  
  bar_data_down <- subset(bar_data, condition == 'Downregulated')
  bar_data_down <- head(bar_data_down[order(bar_data_down$p.adjust),],15)
  
  bar_data_reduced <- rbind(bar_data_up, bar_data_down)
  
  bar_data_ordered <- bar_data_reduced[order(bar_data_reduced$p.adjust),] # order by count
  bar_data_ordered<- bar_data_ordered[order(bar_data_ordered$category),] # order by category
  bar_data_ordered$num <- seq(1:nrow(bar_data_ordered)) # num category for plot
  
  # save dataset
  print("Save dataset in RData")
  save(bar_data_ordered, file = paste0(outdir, "all_GO_", plot_name, ".RData"))

  # join colors
  bar_data_ordered_mod <- left_join(bar_data_ordered, Category_colors, by= "category")
  
  print("All GO terms -  graph")
  g <- ggplot(bar_data_ordered_mod, aes(count, reorder(term, -num), fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = count),
      color = "black",
      hjust = -0.1,
      size = 2,
      position = position_dodge(0.9)
    ) +
    labs(x = "Gene counts" , y = NULL) +
    scale_fill_manual(name='Category', 
                      labels = unique(bar_data_ordered_mod$label), 
                      values = unique(bar_data_ordered_mod$colors)) +
    theme(
      legend.position = "right",
      # panel.grid = element_blank(),
      # axis.text.x = element_blank(),
      # axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      strip.text.x = element_text(size = 5, face = "bold"),
      strip.background = element_blank()
    )
  
  print("Save graph")
  ggsave(paste0(figdir, "barplotGO_", plot_name, ".png"),
         #plot = g, dpi = 300, width = 1000, height = 800, units = 'px')
         plot = g + theme_classic(), dpi = 600, width = 10, height = 5)
  
  # save graph
  print("Save graph in RData")
  save(g, file = paste0(figdir, "barplotGO_", plot_name, ".RData"))
  
  
  ## barplot only upregulated
  print("UP GO terms")
  bar_data_up <- subset(bar_data, condition == 'Upregulated')
  bar_data_up <-head(bar_data_up[order(bar_data_up$p.adjust),],40) # order by pvalue
  bar_data_up_ordered <- bar_data_up[order(bar_data_up$p.adjust),] # order by pvalue
  bar_data_up_ordered<- bar_data_up_ordered[order(bar_data_up_ordered$category),] # order by category
  bar_data_up_ordered$p.val <- round(-log10(bar_data_up_ordered$p.adjust), 2)
  bar_data_up_ordered$num <- seq(1:nrow(bar_data_up_ordered)) # num category for plot
  
  # save graph
  print("Save dataset in RData")
  save(bar_data_up_ordered, file = paste0(outdir, "UP_GO_", plot_name, ".RData"))
  
  # join color
  bar_data_up_ordered_mod <- left_join(bar_data_up_ordered, Category_colors, by= "category")
  
  print("UP GO terms -  graph")
  g.up <- ggplot(bar_data_up_ordered_mod, aes(p.val, reorder(term, -num), fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = p.val),
      color = "black",
      hjust = 0,
      size = 2.2,
      position = position_dodge(0)
    ) +
    labs(x = "-log10(p-value)" , y = NULL) +
    scale_fill_manual(name='Category', 
                      labels = unique(bar_data_up_ordered_mod$label), 
                      values = unique(bar_data_up_ordered_mod$colors)) +
    theme(
      legend.position = "right",
      # panel.grid = element_blank(),
      # axis.text.x = element_blank(),
      # axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      strip.text.x = element_text(size = 11, face = "bold"),
      strip.background = element_blank() 
    ) + theme_classic()
  
  print("Save graph")
  ggsave(paste0(figdir, "barplotUP_GO_", plot_name, ".png"),
         plot = g.up + theme_classic(), dpi = 600, width = 10, height = 5)
  
  # save graph
  print("Save graph in RData")
  save(g.up, file = paste0(figdir, "barplotUP_GO", plot_name, ".RData"))
  
  
  ## barplot only downregulated
  print("DOWN GO terms")
  bar_data_down <- subset(bar_data, condition == 'Downregulated')
  bar_data_down <-head(bar_data_down[order(bar_data_down$p.adjust),],40) # order by pvalue
  bar_data_down_ordered <- bar_data_down[order(bar_data_down$p.adjust),] # order by pvalue
  bar_data_down_ordered<- bar_data_down_ordered[order(bar_data_down_ordered$category),] # order by category
  bar_data_down_ordered$p.val <- round(-log10(bar_data_down_ordered$p.adjust), 2)
  bar_data_down_ordered$num <- seq(1:nrow(bar_data_down_ordered)) # num category for plot
  
  # save dataset
  print("Save dataset in RData")
  save(bar_data_down_ordered, file = paste0(outdir, "DOWN_GO_", plot_name, ".RData"))
  
  # join color
  bar_data_down_ordered_mod <- left_join(bar_data_down_ordered, Category_colors, by= "category")
  
  print("DOWN GO terms -  graph")
  g.down <- ggplot(bar_data_down_ordered_mod, aes(p.val, reorder(term, -num), fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = p.val),
      color = "black",
      hjust = 0,
      size = 2.2,
      position = position_dodge(0)
    ) +
    labs(x = "-log10(p-value)" , y = NULL) +
    scale_fill_manual(name='Category', 
                      labels = unique(bar_data_down_ordered_mod$label), 
                      values = unique(bar_data_down_ordered_mod$colors)) +
    theme(
      legend.position = "right",
      axis.title.y = element_blank(),
      strip.text.x = element_text(size = 11, face = "bold"),
      strip.background = element_blank()
    )+ theme_classic()
  
  print("Save graph")
  ggsave(paste0(figdir,"barplotDOWN_GO_", plot_name, ".png"),
         plot = g.down + theme_classic(), dpi = 600, width = 10, height = 5)
  # dev.off()
  
  # save graph
  print("Save graph in RData")
  save(g.down, file = paste0(figdir, "barplotDOWN_GO_", plot_name, ".RData"))
  
  } #end for
}#end function


# sessionInfo()
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
#   [1] forcats_0.5.1          stringr_1.4.0          dplyr_1.0.10
# [4] purrr_0.3.4            readr_1.4.0            tidyr_1.2.1
# [7] tibble_3.1.3           tidyverse_1.3.0        ggplot2_3.3.5
# [10] clusterProfiler_3.18.1 DOSE_3.16.0            enrichplot_1.10.2
# [13] gprofiler2_0.2.0
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-7         fs_1.5.0             lubridate_1.7.9.2
# [4] bit64_4.0.5          RColorBrewer_1.1-2   httr_1.4.2
# [7] tools_4.0.2          backports_1.2.1      utf8_1.2.2
# [10] R6_2.5.0             DBI_1.1.1            lazyeval_0.2.2
# [13] BiocGenerics_0.36.1  colorspace_2.0-2     withr_2.4.2
# [16] tidyselect_1.2.0     gridExtra_2.3        bit_4.0.4
# [19] compiler_4.0.2       rvest_0.3.6          cli_3.6.0
# [22] Biobase_2.50.0       xml2_1.3.2           scatterpie_0.1.6
# [25] plotly_4.9.3         labeling_0.4.2       shadowtext_0.0.8
# [28] scales_1.1.1         digest_0.6.27        pkgconfig_2.0.3
# [31] htmltools_0.5.1.1    dbplyr_2.2.1         fastmap_1.1.0
# [34] readxl_1.3.1         htmlwidgets_1.5.3    rlang_1.0.6
# [37] rstudioapi_0.13      RSQLite_2.2.7        farver_2.1.0
# [40] generics_0.1.0       jsonlite_1.7.2       BiocParallel_1.24.1
# [43] GOSemSim_2.16.1      RCurl_1.98-1.3       magrittr_2.0.1
# [46] GO.db_3.12.1         Matrix_1.3-4         Rcpp_1.0.7
# [49] munsell_0.5.0        S4Vectors_0.28.1     fansi_0.5.0
# [52] viridis_0.5.1        lifecycle_1.0.3      stringi_1.6.2
# [55] ggraph_2.0.5         MASS_7.3-53          plyr_1.8.6
# [58] qvalue_2.22.0        grid_4.0.2           blob_1.2.1
# [61] parallel_4.0.2       ggrepel_0.9.1        DO.db_2.9
# [64] crayon_1.4.1         lattice_0.20-41      graphlayouts_0.7.1
# [67] haven_2.3.1          cowplot_1.1.1        splines_4.0.2
# [70] hms_1.1.0            pillar_1.6.2         fgsea_1.16.0
# [73] igraph_1.2.8         reshape2_1.4.4       stats4_4.0.2
# [76] fastmatch_1.1-0      reprex_1.0.0         glue_1.4.2
# [79] downloader_0.4       modelr_0.1.8         data.table_1.14.0
# [82] BiocManager_1.30.21  vctrs_0.5.1          tweenr_1.0.2
# [85] cellranger_1.1.0     gtable_0.3.0         polyclip_1.10-0
# [88] assertthat_0.2.1     cachem_1.0.5         ggforce_0.3.3
# [91] broom_0.7.9          tidygraph_1.2.0      viridisLite_0.4.0
# [94] rvcheck_0.1.8        AnnotationDbi_1.52.0 memoise_2.0.0
# [97] IRanges_2.24.1       ellipsis_0.3.2