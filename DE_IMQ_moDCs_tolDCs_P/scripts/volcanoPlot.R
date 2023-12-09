# Volcano Plot moDCs VS tolDCs in patients
# Author: Sofia Salazar
# December 2023
# -----

# Setup
# qlogin
# module load r/4.0.2
# R

# ----Libraries----
library(ggplot2)
library(tidyverse)

# ----Load data and function----

workdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/'
indir <- paste0(workdir, 'DE_IMQ_moDC_tolDC_P/DE_files/')
figdir <- paste0(workdir, 'DE_IMQ_moDC_tolDC_P/figures/')
# Function
source(paste0(workdir, 'plot_Volcano_function.R'))

# ----Create volcano plot----

make_volcano(indir, figdir, abslogFC = 1, noGeneNames = 20)
