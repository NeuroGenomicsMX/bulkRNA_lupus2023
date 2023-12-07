# Plot Volcano
# Author: Sofia Salazar
# Date: 26 Nov 2023
# Plot volcano plots of DE results
# ----

# Cluster setup
# qlogin
# module load r/4.0.2
# R

# ----
# Libraries
library(ggplot2)
library(tidyverse)
# ----
# Load data
# change to place where "DE_plotName.csv" files are:
indir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_P_C_celltypes/DE_files/' 

setwd(indir)
files <- dir(getwd(), pattern = "^DE_*")
# ---
# Define plotting function
plot_volcano <- function(file_name, figdir, abslogFC = 1, noGeneNames = 20){ 
  # ---
  # file_name: "DE_plotName.csv" with DESeq results are, the plotName part will be used as plot title
  # figdir: place where volcano plots will be saved
  # abslogFC: threshold of logFC (absolute value) for coloring the DE genes
  # noGeneNames: maximum number of genes to write their name for each group (downregulated/upregulated)
  # ---
  plot_name <- gsub("^DE_(.+)\\.csv$", "\\1", file_name)
  
  df <- read.csv(file = paste0(file_name), row.names = 'X')
  df <- na.omit(df)
  # --
  # Add expression column
  df <- df %>% 
    mutate(Expression = case_when(log2FoldChange >= abslogFC & padj < 0.05 ~ "Up-regulated",
                                  log2FoldChange <= -(abslogFC) & padj < 0.05 ~ "Down-regulated",
                                  TRUE ~ "Unchanged"))
  # --
  # Plot
  volcanoplot <- ggplot(df, aes(log2FoldChange, -log(padj,10))) +
    geom_point(aes(color = Expression), size = 0.7) +
    labs(title = plot_name) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"p-adj")) +
    scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    geom_vline(xintercept = abslogFC, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_vline(xintercept = -(abslogFC), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5)
  
  
  top <- noGeneNames # no. of highlighted genes in plot
  top_genes <- bind_rows(df %>% 
                           filter(Expression == 'Up-regulated') %>% 
                           arrange(padj, desc(abs(log2FoldChange))) %>% 
                           head(top),
                         df %>% 
                           filter(Expression == 'Down-regulated') %>% 
                           arrange(padj, desc(abs(log2FoldChange))) %>% 
                           head(top)
  )
  dim(df[df$Expression == 'Up-regulated',])
  dim(df[df$Expression == 'Down-regulated',])
  
  volcanoplot_names <-  volcanoplot +
    ggrepel::geom_label_repel(data = top_genes,
                              mapping = aes(log2FoldChange, -log(padj,10), label = rownames(top_genes)),
                              size = 2) + theme_classic() + theme(legend.position = 'bottom')
  ggsave(file = paste0(figdir, "volcano_", plot_name,".png"), plot = volcanoplot_names, dpi = 300)
  # save Rdata object
  save(volcanoplot_names, file =paste0(figdir, "volcano_", plot_name, ".RData"))
}

# ---
# make plots


for (f in files){
  print(f)
  plot_volcano(f, figdir = '/mnt/Citosina/amedina/lupus/RNA_lupus/DE_P_C_celltypes/figures/',
               abslogFC = 1,
               noGeneNames = 20)
}

