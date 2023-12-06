# DE types - Treatment P
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
figdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_P/figures/'
outdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_Type_Treatment_P/'
load(file = paste0(outdir, 'DE_SLE_results.RData'))

# --- Contrast 1 (monocytes vs moDC) -----
# add a column of NAs
SLE_moVSmoDC_df <- as.data.frame(res_SLE_moVSmoDC_Ordered)
SLE_moVSmoDC_df$diffexpressed <- "NO"

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
SLE_moVSmoDC_df$diffexpressed[SLE_moVSmoDC_df$log2FoldChange > 1 & SLE_moVSmoDC_df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
SLE_moVSmoDC_df$diffexpressed[SLE_moVSmoDC_df$log2FoldChange < -1 & SLE_moVSmoDC_df$pvalue < 0.05] <- "DOWN"

# Create a new column "names" to de, that will contain the name of a 
# subset if genes differentially expressed (NA in case they are not)
SLE_moVSmoDC_df$names <- NA

# filter for a subset of interesting genes
filter <- which(SLE_moVSmoDC_df$diffexpressed != "NO" & SLE_moVSmoDC_df$padj < 0.05 & (SLE_moVSmoDC_df$log2FoldChange >= 5  | SLE_moVSmoDC_df$log2FoldChange <= -5))
SLE_moVSmoDC_df$names[filter] <- rownames(SLE_moVSmoDC_df)[filter]

# Graph
png(file = paste0(figdir, "volcano_SLE_moVSmoDC.png"), width = 800, height = 800) # save plot
ggplot(data=SLE_moVSmoDC_df, aes(x=log2FoldChange, y=-log10(pvalue), 
                                     col=diffexpressed, label=names)) +
  geom_point() +
  scale_color_manual(values=c("blue", "black", "red")) + # colors
  theme_minimal() +
  geom_text_repel() # +
#xlim(-15,15)

dev.off() #end save plot

# --- Contrast 2 (monocytes vs moDC+IMQ) -----
# add a column of NAs
SLE_moVSmoDCIMQ_df <- as.data.frame(res_SLE_moVSmoDCIMQ_Ordered)
SLE_moVSmoDCIMQ_df$diffexpressed <- "NO"

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
SLE_moVSmoDCIMQ_df$diffexpressed[SLE_moVSmoDCIMQ_df$log2FoldChange > 1 & SLE_moVSmoDCIMQ_df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
SLE_moVSmoDCIMQ_df$diffexpressed[SLE_moVSmoDCIMQ_df$log2FoldChange < -1 & SLE_moVSmoDCIMQ_df$pvalue < 0.05] <- "DOWN"

# Create a new column "names" to de, that will contain the name of a 
# subset if genes differentially expressed (NA in case they are not)
SLE_moVSmoDCIMQ_df$names <- NA

# filter for a subset of interesting genes
filter <- which(SLE_moVSmoDCIMQ_df$diffexpressed != "NO" & SLE_moVSmoDCIMQ_df$padj < 0.05 & (SLE_moVSmoDCIMQ_df$log2FoldChange >= 5  | SLE_moVSmoDCIMQ_df$log2FoldChange <= -5))
SLE_moVSmoDCIMQ_df$names[filter] <- rownames(SLE_moVSmoDCIMQ_df)[filter]

# Graph
png(file = paste0(figdir, "volcano_SLE_moVSmoDCIMQ.png"), width = 800, height = 800) # save plot
ggplot(data=SLE_moVSmoDCIMQ_df, aes(x=log2FoldChange, y=-log10(pvalue), 
                                        col=diffexpressed, label=names)) +
  geom_point() +
  scale_color_manual(values=c("blue", "black", "red")) + # colors
  theme_minimal() +
  geom_text_repel() 

dev.off() #end save plot

# --- Contrast 3 (monocytes vs tolDC) -----
# add a column of NAs
SLE_moVStolDC_df <- as.data.frame(res_SLE_moVStolDC_Ordered)
SLE_moVStolDC_df$diffexpressed <- "NO"

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
SLE_moVStolDC_df$diffexpressed[SLE_moVStolDC_df$log2FoldChange > 1 & SLE_moVStolDC_df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
SLE_moVStolDC_df$diffexpressed[SLE_moVStolDC_df$log2FoldChange < -1 & SLE_moVStolDC_df$pvalue < 0.05] <- "DOWN"

# Create a new column "names" to de, that will contain the name of a 
# subset if genes differentially expressed (NA in case they are not)
SLE_moVStolDC_df$names <- NA

# filter for a subset of interesting genes
filter <- which(SLE_moVStolDC_df$diffexpressed != "NO" & SLE_moVStolDC_df$padj < 0.05 & (SLE_moVStolDC_df$log2FoldChange >= 5  | SLE_moVStolDC_df$log2FoldChange <= -5))
SLE_moVStolDC_df$names[filter] <- rownames(SLE_moVStolDC_df)[filter]

# Graph
png(file = paste0(figdir, "volcano_SLE_moVStolDC.png"), width = 800, height = 800) # save plot
ggplot(data=SLE_moVStolDC_df, aes(x=log2FoldChange, y=-log10(pvalue), 
                                      col=diffexpressed, label=names)) +
  geom_point() +
  scale_color_manual(values=c("blue", "black", "red")) + # colors
  theme_minimal() +
  geom_text_repel() 

dev.off() #end save plot

# --- Contrast 4 (monocytes vs tolDC+IMQ) -----
# add a column of NAs
SLE_moVStolDCIMQ_df <- as.data.frame(res_SLE_moVStolDCIMQ_Ordered)
SLE_moVStolDCIMQ_df$diffexpressed <- "NO"

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
SLE_moVStolDCIMQ_df$diffexpressed[SLE_moVStolDCIMQ_df$log2FoldChange > 1 & SLE_moVStolDCIMQ_df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
SLE_moVStolDCIMQ_df$diffexpressed[SLE_moVStolDCIMQ_df$log2FoldChange < -1 & SLE_moVStolDCIMQ_df$pvalue < 0.05] <- "DOWN"

# Create a new column "names" to de, that will contain the name of a 
# subset if genes differentially expressed (NA in case they are not)
SLE_moVStolDCIMQ_df$names <- NA

# filter for a subset of interesting genes
filter <- which(SLE_moVStolDCIMQ_df$diffexpressed != "NO" & SLE_moVStolDCIMQ_df$padj < 0.05 & (SLE_moVStolDCIMQ_df$log2FoldChange >= 5  | SLE_moVStolDCIMQ_df$log2FoldChange <= -5))
SLE_moVStolDCIMQ_df$names[filter] <- rownames(SLE_moVStolDCIMQ_df)[filter]

# Graph
png(file = paste0(figdir, "volcano_SLE_moVStolDCIMQ.png"), width = 800, height = 800) # save plot
ggplot(data=SLE_moVStolDCIMQ_df, aes(x=log2FoldChange, y=-log10(pvalue), 
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
