# UpSetR - plots
# Author: Evelia Coss
# Date: 7 Dec 2023
# Make your own UpSetR
# ---

# Cluster setup
# qlogin
# module load r/4.0.2
# R

# ---- Libraries ----
library(pacman)

pacman::p_load(
  dplyr,  # Manipulacion de datos
  ggplot, # Graficas
  UpSetR  # UpsetR
)

# Or get the latest version directly from GitHub
#devtools::install_github("const-ae/ggupset")
#library(ggupset) # UpsetR

library(ComplexHeatmap)

# LAVIS sever
#indir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/DE_files/' 
# My computer
indir <- "C:/Users/ecoss/OneDrive - CINVESTAV/Documentos/Posdoc_LIIGH/LUPUS_Proyecto_2023/bulkRNA_lupus2023/DE_Type_Treatment_C/DE_files/"
load(paste0(outdir, 'DE_control_results.RData'))

# figdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_P_C_celltypes/figures/'
figdir <- 'C:/Users/ecoss/OneDrive - CINVESTAV/Documentos/Posdoc_LIIGH/LUPUS_Proyecto_2023/bulkRNA_lupus2023/DE_Type_Treatment_C/figures/'

# ---- UpsetR ----

# Funcion 1
ConvertData_UpsetR <- function(indir, abslogFC = 1){ 
  results <- list() # create a empty list
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
  up_names <- df %>%
    filter(Expression == "Up-regulated")
  up_names <- up_names$genes # only names
  
  # Down names
  down_names <- df %>%
    filter(Expression == "Down-regulated")
  down_names <- down_names$genes # only names

  # Results
  name1 = paste(plot_name, "UP", sep="_")
  #print(name1)
  name2 = paste(plot_name, "DOWN", sep="_")
  #print(name2)
  
  # add information
  tmp <- list( name1 = up_names)
  names(tmp) <- name1 #renames
  results <- c(results, tmp)
  tmp <- list( name2 = down_names)
  names(tmp) <- name2 #renames
  results <- c(results, tmp)
  
  } #end loop
  return(results)

} # end function





