# Plot Volcano
# Author: Sofia Salazar
# Date: 26 Nov 2023
# Plot volcano plots of DE results
# ----

# Cluster setup
# qlogin
# module load r/4.0.2
# R

# ---- Libraries ---
library(ggplot2)
library(tidyverse)

# ---- Load data and function -----
# change to place where "DE_plotName.csv" files are:
# 1. Dataset with names
indir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/DE_files/'
figdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/figures/'
  
# function
source(paste0(indir, "plot_Volcano_function.R"))

# ---- Run multiplots -----
# Usage: make_volcano(indir, figdir, abslogFC = 1, noGeneNames = 20)

# Arguments
# indir = input directory
# figdir =  output directory for figures
# abslogFC = log2FoldChange cutoff value, default = 1
# noGeneNames = Number of genes with higher expression to name them, default = 20

make_volcano(indir, figdir, abslogFC = 1, noGeneNames = 20)




