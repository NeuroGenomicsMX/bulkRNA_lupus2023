# DE types - Treatment C
# Volcano plot
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
library(ggplot2)
library(ggrepel)

# --- Load data -----
indir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/'
figdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/figures/'
outdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_C/'
load(file = paste0(outdir, 'DE_control_results.RData'))

# --- Contrast 1 (monocytes vs moDC) -----
# add a column of NAs
control_moVSmoDC_df <- as.data.frame(res_control_moVSmoDC_Ordered)
control_moVSmoDC_df$diffexpressed <- "NO"

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
control_moVSmoDC_df$diffexpressed[control_moVSmoDC_df$log2FoldChange > 1 & control_moVSmoDC_df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
control_moVSmoDC_df$diffexpressed[control_moVSmoDC_df$log2FoldChange < -1 & control_moVSmoDC_df$pvalue < 0.05] <- "DOWN"

# Create a new column "names" to de, that will contain the name of a 
# subset if genes differentially expressed (NA in case they are not)
control_moVSmoDC_df$names <- NA

# filter for a subset of interesting genes
filter <- which(control_moVSmoDC_df$diffexpressed != "NO" & control_moVSmoDC_df$padj < 0.05 & (control_moVSmoDC_df$log2FoldChange >= 5  | control_moVSmoDC_df$log2FoldChange <= -5))
control_moVSmoDC_df$names[filter] <- rownames(control_moVSmoDC_df)[filter]

# Graph
png(file = paste0(figdir, "volcano_control_moVSmoDC.png"), width = 800, height = 800) # save plot
ggplot(data=control_moVSmoDC_df, aes(x=log2FoldChange, y=-log10(pvalue), 
                                     col=diffexpressed, label=names)) +
  geom_point() +
  scale_color_manual(values=c("blue", "black", "red")) + # colors
  theme_minimal() +
  geom_text_repel() # +
  #xlim(-15,15)

dev.off() #end save plot

# --- Contrast 2 (monocytes vs moDC+IMQ) -----
# add a column of NAs
control_moVSmoDCIMQ_df <- as.data.frame(res_control_moVSmoDCIMQ_Ordered)
control_moVSmoDCIMQ_df$diffexpressed <- "NO"

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
control_moVSmoDCIMQ_df$diffexpressed[control_moVSmoDCIMQ_df$log2FoldChange > 1 & control_moVSmoDCIMQ_df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
control_moVSmoDCIMQ_df$diffexpressed[control_moVSmoDCIMQ_df$log2FoldChange < -1 & control_moVSmoDCIMQ_df$pvalue < 0.05] <- "DOWN"

# Create a new column "names" to de, that will contain the name of a 
# subset if genes differentially expressed (NA in case they are not)
control_moVSmoDCIMQ_df$names <- NA

# filter for a subset of interesting genes
filter <- which(control_moVSmoDCIMQ_df$diffexpressed != "NO" & control_moVSmoDCIMQ_df$padj < 0.05 & (control_moVSmoDCIMQ_df$log2FoldChange >= 5  | control_moVSmoDCIMQ_df$log2FoldChange <= -5))
control_moVSmoDCIMQ_df$names[filter] <- rownames(control_moVSmoDCIMQ_df)[filter]

# Graph
png(file = paste0(figdir, "volcano_control_moVSmoDCIMQ.png"), width = 800, height = 800) # save plot
ggplot(data=control_moVSmoDCIMQ_df, aes(x=log2FoldChange, y=-log10(pvalue), 
        col=diffexpressed, label=names)) +
  geom_point() +
  scale_color_manual(values=c("blue", "black", "red")) + # colors
  theme_minimal() +
  geom_text_repel() 

dev.off() #end save plot

# --- Contrast 3 (monocytes vs tolDC) -----
# add a column of NAs
control_moVStolDC_df <- as.data.frame(res_control_moVStolDC_Ordered)
control_moVStolDC_df$diffexpressed <- "NO"

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
control_moVStolDC_df$diffexpressed[control_moVStolDC_df$log2FoldChange > 1 & control_moVStolDC_df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
control_moVStolDC_df$diffexpressed[control_moVStolDC_df$log2FoldChange < -1 & control_moVStolDC_df$pvalue < 0.05] <- "DOWN"

# Create a new column "names" to de, that will contain the name of a 
# subset if genes differentially expressed (NA in case they are not)
control_moVStolDC_df$names <- NA

# filter for a subset of interesting genes
filter <- which(control_moVStolDC_df$diffexpressed != "NO" & control_moVStolDC_df$padj < 0.05 & (control_moVStolDC_df$log2FoldChange >= 5  | control_moVStolDC_df$log2FoldChange <= -5))
control_moVStolDC_df$names[filter] <- rownames(control_moVStolDC_df)[filter]

# Graph
png(file = paste0(figdir, "volcano_control_moVStolDC.png"), width = 800, height = 800) # save plot
ggplot(data=control_moVStolDC_df, aes(x=log2FoldChange, y=-log10(pvalue), 
                                        col=diffexpressed, label=names)) +
  geom_point() +
  scale_color_manual(values=c("blue", "black", "red")) + # colors
  theme_minimal() +
  geom_text_repel() 

dev.off() #end save plot

# --- Contrast 4 (monocytes vs tolDC+IMQ) -----
# add a column of NAs
control_moVStolDCIMQ_df <- as.data.frame(res_control_moVStolDCIMQ_Ordered)
control_moVStolDCIMQ_df$diffexpressed <- "NO"

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
control_moVStolDCIMQ_df$diffexpressed[control_moVStolDCIMQ_df$log2FoldChange > 1 & control_moVStolDCIMQ_df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
control_moVStolDCIMQ_df$diffexpressed[control_moVStolDCIMQ_df$log2FoldChange < -1 & control_moVStolDCIMQ_df$pvalue < 0.05] <- "DOWN"

# Create a new column "names" to de, that will contain the name of a 
# subset if genes differentially expressed (NA in case they are not)
control_moVStolDCIMQ_df$names <- NA

# filter for a subset of interesting genes
filter <- which(control_moVStolDCIMQ_df$diffexpressed != "NO" & control_moVStolDCIMQ_df$padj < 0.05 & (control_moVStolDCIMQ_df$log2FoldChange >= 5  | control_moVStolDCIMQ_df$log2FoldChange <= -5))
control_moVStolDCIMQ_df$names[filter] <- rownames(control_moVStolDCIMQ_df)[filter]

# Graph
png(file = paste0(figdir, "volcano_control_moVStolDCIMQ.png"), width = 800, height = 800) # save plot
ggplot(data=control_moVStolDCIMQ_df, aes(x=log2FoldChange, y=-log10(pvalue), 
                                      col=diffexpressed, label=names)) +
  geom_point() +
  scale_color_manual(values=c("blue", "black", "red")) + # colors
  theme_minimal() +
  geom_text_repel() 

dev.off() #end save plot

# --- for loop (trabajando en esto)-----
# Usage: volcano_plot(res_dds, LFG1 = 1.5, pv = 0.05, LFG2 = 3, 
# output = "volcano.png", width = 800, height = 800)

# load function
# list.files(indir)
# source(paste0(indir,"volcano_plot_function.R"))
# 
# volcano_plot(res_control_moVSmoDC_Ordered, output = paste0(outdir,"Control_moVSmoDC_vp.png"))
# 
